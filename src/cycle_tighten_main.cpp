/*
 *  cycle_tighten_main.cpp
 *  mplp
 *
 *  Created by Amir Globerson and David Sontag on 8/10/08.
 *  Updated by David Sontag, Do Kook Choe, and Yitao Li in 2012.
 *  Copyright 2008 MIT, 2012 NYU. All rights reserved.
 *
 *  Notes:
 *  Setting the environmental INF_TIME to the total number of seconds allowed to
 *  run will result in global decoding being called once 1/3 through, and (if turned
 *  on) decimation being called 2/3 through (very helpful for CSP intances).
 *  We did not use this for the UAI 2012 paper (i.e., we did not set INF_TIME).
 */
#include <iostream>
#include <ctime>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>

#include <MPLP/cycle.h> // Most of the logic is in here.

using namespace std;
using namespace mplpLib;

#define MPLP_MAX_TIGHT_ITERS 1000000 // 100 used for testing

int main( int argc, char *argv[] ){
    // Note: MPLP_DEBUG_MODE flag can be found in "cycle.h"
    // DEFAULT values:
    MPLPIndexType niter = 1000, niter_later = 20, nclus_to_add_min = 5, nclus_to_add_max = 20;
    double obj_del_thr = .0002, int_gap_thr = .0002;
    bool UAIsettings = false;  // Settings for UAI inference competition override all others

    // defaults. UAIsettings modulates the value of many of these
    bool addEdgeIntersections = true;
    bool doGlobalDecoding = false;
    bool useDecimation=false;
    bool lookForCSPs = false;

    if(UAIsettings) {
        doGlobalDecoding = true;
        useDecimation=true;
        lookForCSPs = true;
    }

    if (argc < 2 || (argc > 2 && argc < 5) || argc > 6) {
        printf("Syntax: ./solver input-model-file input-evidence-file random-seed MPE log-file\n");
        return -1;
    }

    double time_elapsed, /*time,*/ time_limit;
    bool LOG_MODE=false; // default

    clock_t start = clock();
    // TODO: consider checking INF_TIME and changing parameters as a result
    char *t = getenv("INF_TIME");
    if (!t) {
        // Time limit. Also affects when global decoding & decimation are called.
        time_limit = 99999999;
    } else {
        time_limit = (double)atoi(t);
    }

    if(MPLP_DEBUG_MODE)
        cout << "Time limit = " << time_limit << endl;

    bool decimation_has_started = false;
    bool force_decimation = false;
    bool prevGlobalDecodingWas1 = true;

    // Keep track of triplets added so far
    map<vector<MPLPIndexType>, bool> triplet_set;

    /*      // We probably do not need to worry about memory limit yet?
	char *m = getenv("INF_MEMORY");
	if (!m){
		std::cerr<<"Error: environmental variable INF_MEMORY is missing"<<std::endl;
	}
	int mem_limit = atoi(m);
     */

    std::string input_file(argv[1]), evidence_file;
    if(argc > 2)
        evidence_file = argv[2];
    else
        evidence_file = "no.evid"; // default if not provided as input

    unsigned int seed;
    if(argc > 2)
        sscanf(argv[3],"%u",&seed);
    else seed = 0; // default if not provided as input

    if(MPLP_DEBUG_MODE) cout << "Random seed = " << seed << endl;
    srand(seed);

    // currently not used (since we only perform MPE)
    //	sscanf(argv[4],"%s",task);

    // Log file
    FILE *log_file = 0;
    if(argc == 6) {
        LOG_MODE = true;
        std::string log_filename(argv[5]);
        log_file = fopen(log_filename.c_str(), "w");
    }

    if(LOG_MODE) {
        fprintf(log_file, "I niter=%lu, niter_later=%lu, nclus_to_add_min=%lu, nclus_to_add_max=%lu, obj_del_thr=%lg, int_gap_thr=%lg\n", niter, niter_later, nclus_to_add_min, nclus_to_add_max, obj_del_thr,int_gap_thr);
    }
    if (MPLP_DEBUG_MODE) cout << "niter=" << niter << "\nniter_later=" << niter_later << "\nnclus_to_add=" << nclus_to_add_min << "\nobj_del_thr=" << obj_del_thr << "\nint_gap_thr=" << int_gap_thr << endl;

    // Load in the MRF and initialize GMPLP state
    MPLPAlg mplp(start, time_limit, input_file, evidence_file, log_file, lookForCSPs);

    if (MPLP_DEBUG_MODE) cout << "Initially running MPLP for " << niter << " iterations" << endl;
    mplp.RunMPLP(niter, obj_del_thr, int_gap_thr);

    for(MPLPIndexType iter=1; iter<MPLP_MAX_TIGHT_ITERS; iter++){  // Break when problem is solved
        if(LOG_MODE) fflush(log_file);
        if (MPLP_DEBUG_MODE) cout << "\n\nOuter loop iteration "<< iter << "\n----------------------" << endl;

        // Is problem solved? If so, break.
        double int_gap = mplp.last_obj - mplp.m_best_val;
        if(int_gap < int_gap_thr){
            if (MPLP_DEBUG_MODE) printf("Done! Integrality gap less than %lg\n", int_gap_thr);
            break;
        }

        // Heuristic: when the integrality gap is sufficiently small, allow the algorithm
        // more time to run till convergence

        if(int_gap < 1){
            niter_later = max(niter_later, MPLPIndexType(600));  // TODO opt: don't hard code
            obj_del_thr = min(obj_del_thr, 1e-5);
            if (MPLP_DEBUG_MODE) cout << "Int gap small, so setting niter_later to " << niter_later << " and obj_del_thr to " << obj_del_thr << endl;
        }

        // Keep track of global decoding time and run this frequently, but at most 20% of total runtime
        if(doGlobalDecoding && (((double)clock() - mplp.last_global_decoding_end_time)/CLOCKS_PER_SEC >= mplp.last_global_decoding_total_time*4)) {
            // Alternate between global decoding methods
            if(prevGlobalDecodingWas1) {
                mplp.RunGlobalDecoding(false);
                prevGlobalDecodingWas1 = false;
            }
            else {
                mplp.RunGlobalDecoding2(false);
                prevGlobalDecodingWas1 = true;
            }
        }

        // Tighten LP
        if (MPLP_DEBUG_MODE) cout << "Now attempting to tighten LP relaxation..." << endl;

        clock_t tightening_start_time = clock();
        double bound=0; double bound2 = 0;
        MPLPIndexType nClustersAdded = 0;

        nClustersAdded += TightenTriplet(mplp, nclus_to_add_min, nclus_to_add_max, triplet_set, bound);
        nClustersAdded += TightenCycle(mplp, nclus_to_add_min, triplet_set, bound2, 1);

        if(max(bound, bound2) < MPLP_CLUSTER_THR) {

            if(MPLP_DEBUG_MODE)
                cout << "TightenCycle did not find anything useful! Re-running with FindPartition." << endl;

            nClustersAdded += TightenCycle(mplp, nclus_to_add_min, triplet_set, bound2, 2);
        }

        // Check to see if guaranteed bound criterion was non-trivial.
        // TODO: these bounds are not for the cycles actually added (since many of the top ones are skipped, already being in the relaxation). Modify it to be so.
        bool noprogress = false;
        if(max(bound, bound2) < MPLP_CLUSTER_THR)
            noprogress = true;

        clock_t tightening_end_time = clock();
        double tightening_total_time = (double)(tightening_end_time - tightening_start_time)/CLOCKS_PER_SEC;
        if (MPLP_DEBUG_MODE) {
            cout << " -- Added " << nClustersAdded << " clusters to relaxation. Took " << tightening_total_time << " seconds" << endl;
        }
        if(LOG_MODE) {
            fprintf(log_file, "I added %lu clusters. Took %lg seconds\n", nClustersAdded, tightening_total_time);
        }

        // For CSP instances, 2/3 through run time, start decimation -- OR, when no progress being made
        if((mplp.CSP_instance || noprogress) && ((double)(clock() - start) / CLOCKS_PER_SEC) > time_limit*2/3)
            force_decimation = true;

        /*
      We have done as much as we can with the existing edge intersection sets. Now
      add in all new edge intersection sets for large clusters.
         */
        if(nClustersAdded == 0 && addEdgeIntersections) {
            mplp.AddAllEdgeIntersections();
            addEdgeIntersections = false; // only makes sense to run this code once
        }

        // Not able to tighten relaxation further, so try to see if decoding is the problem
        // Do not run this too often!
        else if((!addEdgeIntersections && nClustersAdded == 0) || force_decimation) {

            // Do one last push to try to find the global assignment!
            if(doGlobalDecoding && (!useDecimation || !decimation_has_started))
                mplp.RunGlobalDecoding3();

            // Do one step of decimation
            if (useDecimation) {
                decimation_has_started = true;

                bool fixed_node = mplp.RunDecimation();
                if(!fixed_node) {
                    if(MPLP_DEBUG_MODE)
                        cout << "Decimation fixed all of the nodes it could... quiting." << endl;
                    break;
                }
            }
        }

        if (MPLP_DEBUG_MODE) cout << "Running MPLP again for " << niter_later << " more iterations" << endl;
        mplp.RunMPLP(niter_later, obj_del_thr, int_gap_thr);

        if(UAIsettings) {
            // For UAI competition: time limit can be up to 1 hour, so kill process if still running.
            time_elapsed = (double)(clock() - start)/ CLOCKS_PER_SEC;
            if (time_elapsed > 4000 && time_elapsed > time_limit + 60) {
                break;    // terminates if alreay running past time limit (this should be very conservative)
            }
        }

        if(LOG_MODE) fflush(log_file);
    }

    if(LOG_MODE) fflush(log_file);
    if(LOG_MODE) fclose(log_file);

    return 0;  //the solver returns 0 as normal exit status
}
