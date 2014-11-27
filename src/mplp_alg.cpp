/*
 *  mplp_alg.cpp
 *  mplp
 *
 *  Created by Amir Globerson and David Sontag on 6/25/08.
 *  Updated by David Sontag, Do Kook Choe, and Yitao Li in 2012.
 *  Copyright 2008 MIT, 2012 NYU. All rights reserved.
 *
 */

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <list>
#include <queue>
#include <stack>

#include <MPLP/mplp_alg.h>

using namespace std;

#define DEBUG_MODE 1

// Gap used within decoding algorithm. TODO: Better algorithm for choosing this (perhaps iteratively).
#define GAP_THR .001

/////////////////////////////////////////////////////////////////////////////////////
// Code to implement the Region object.
/////////////////////////////////////////////////////////////////////////////////////

mplpLib::Region::Region(vector<int> & region_inds, vector<vector<int> > & all_intersects, vector<int> & intersect_inds, vector<int> & var_sizes, int region_intersect): m_region_inds(region_inds), m_intersect_inds(intersect_inds), m_region_intersect(region_intersect)
{
    // Find the indices of each intersection within the region. Also intialize the message into that intersection
    for (int si=0; si<m_intersect_inds.size(); ++si){
        vector<int> tmp_inds_of_intersects;

        vector<int> curr_intersect = all_intersects[m_intersect_inds[si]];

        vector<int> intersect_var_sizes;
        // Go over all variables in the intersection set
        for (int i=0; i< curr_intersect.size(); ++i){
            int var_in_intersect = curr_intersect[i];
            intersect_var_sizes.push_back(var_sizes[var_in_intersect]);

            vector<int>::iterator iter = find(m_region_inds.begin(), m_region_inds.end(), var_in_intersect);


            // Verify that the current intersection variable indeed appears in this region, and get the index where it appears
            if (iter == m_region_inds.end()){
                cerr << "Intersection set contains variable " << var_in_intersect << " which is not in region" << endl;
                return;
            }
            else{
                tmp_inds_of_intersects.push_back(iter-m_region_inds.begin());
            }
        }
        m_inds_of_intersects.push_back(tmp_inds_of_intersects);

        // This will initialize the message and  set it to zero
        MulDimArr curr_msg(intersect_var_sizes);
        curr_msg = 0;

        m_msgs_from_region.push_back(curr_msg);
    }
    // Calculate the size of the region state space (although this should already be in the lambda object, so we should
    // probably avoid this multiplicity)
    for (int i=0; i<region_inds.size(); ++i) {
        m_var_sizes.push_back(var_sizes[region_inds[i]]);
    }
}

void mplpLib::Region::AddIntersectionSet(int intersect_loc, vector<vector<int> > & all_intersects, vector<int> & var_sizes){
    // Find the indices of intersection within the region. Also intialize the message into that intersection
    vector<int> tmp_inds_of_intersects;
    vector<int> curr_intersect = all_intersects[intersect_loc];

    /*
	if(DEBUG_MODE) {
	  for (int i=0; i< curr_intersect.size(); ++i){
	    int var_in_intersect = curr_intersect[i];
	    cout << "   " << var_in_intersect;
	  }
	  cout << endl;
	}
     */

    vector<int> intersect_var_sizes;
    // Go over all variables in the intersection set
    for (int i=0; i< curr_intersect.size(); ++i){
        int var_in_intersect = curr_intersect[i];
        intersect_var_sizes.push_back(var_sizes[var_in_intersect]);

        vector<int>::iterator iter = find(m_region_inds.begin(), m_region_inds.end(), var_in_intersect);

        // Verify that the current intersection variable indeed appears in this region, and get the index where it appears
        if (iter == m_region_inds.end()){
            cerr << "Intersection set contains variable " << var_in_intersect << " which is not in region" << endl;
            return;
        }else{
            tmp_inds_of_intersects.push_back(iter-m_region_inds.begin());
        }
    }

    /*	if(DEBUG_MODE) {
	  cout << "intersection locs: ";
	  for(int i=0; i<tmp_inds_of_intersects.size(); i++) {

	    cout << tmp_inds_of_intersects[i] << " ";
	  }
	  cout << endl;
	  }*/

    // This will initialize the message and	set it to zero
    MulDimArr curr_msg(intersect_var_sizes);
    curr_msg = 0;

    // Adds new intersection set
    m_intersect_inds.push_back(intersect_loc);
    m_inds_of_intersects.push_back(tmp_inds_of_intersects);
    m_msgs_from_region.push_back(curr_msg);
}

/*
 * This algorithm is a modification of the original MPLP algorithm,
 * described in Figure A-1 on Sontag's Ph.D. thesis (page 104). Rather
 * than add redundant edge intersection sets, it makes use of the
 * original edge intersection sets. As a result there are never any
 * edge->edge messages. This difference is only relevant for
 * tightening (before tightening it is identical to MPLP).
 */
void mplpLib::Region::UpdateMsgs(vector<MulDimArr> & sum_into_intersects)
{
    /* First do the expansion:
	1. Take out the message into the intersection set from the current cluster
	2. Expand it to the size of the region
	3. Add this for all intersection sets
     */
    // Set this to be the region's intersection set value
    MulDimArr orig(sum_into_intersects[m_region_intersect]);
    for (int si=0; si<m_intersect_inds.size(); ++si){
        // Take out previous message
        vector<int> & curr_inds_of_intersect = m_inds_of_intersects[si];

        /*
		if(DEBUG_MODE) {
		  cout << "Intersection set { ";
		  for(int i=0; i < curr_inds_of_intersect.size(); i++) {
		    cout << curr_inds_of_intersect[i] << " ";
		  }
		  cout << "}" << endl;
		}
         */

        m_msgs_from_region[si].ExpandAndAdd(orig, curr_inds_of_intersect);
    }
    // Will store the total messages going into the intersection, but not from the Region
    vector<MulDimArr> lam_minus_region;
    for (int si=0; si<m_intersect_inds.size(); ++si){
        int curr_intersect = m_intersect_inds[si];

        lam_minus_region.push_back(MulDimArr());
        lam_minus_region.back() = sum_into_intersects[curr_intersect];
        lam_minus_region.back() -= m_msgs_from_region[si];

        vector<int> & curr_inds_of_intersect = m_inds_of_intersects[si];

        // If the intersection has the same size as the region, we assume they are the same, and therefore no need
        // to expand. NOTE: This may cause problems if the intersection has the same indices but rearranged.
        if (Get_nVars()==curr_inds_of_intersect.size()){
            sum_into_intersects[m_region_intersect]+= sum_into_intersects[curr_intersect];
        }else{
            sum_into_intersects[curr_intersect].ExpandAndAdd(sum_into_intersects[m_region_intersect],curr_inds_of_intersect);
        }
    }
    // Update messages
    sum_into_intersects[m_region_intersect].max_into_multiple_subsets_special(m_inds_of_intersects,m_msgs_from_region); // sets m_msgs_from_region
    int sC = m_intersect_inds.size();
    for (int si=0; si<m_intersect_inds.size(); ++si){
        // Take out previous message
        int curr_intersect = m_intersect_inds[si];
        // Update message
        m_msgs_from_region[si]*= 1.0/sC;
        // Put in current message
        sum_into_intersects[curr_intersect] = m_msgs_from_region[si];
        // Finish updating message
        // msg_new = new - old + msg_old
        m_msgs_from_region[si] -= lam_minus_region[si];
        // Update region intersection set
        vector<int> & curr_inds_of_intersect = m_inds_of_intersects[si];
        m_msgs_from_region[si].ExpandAndSubtract(orig, curr_inds_of_intersect);
    }
    memcpy(sum_into_intersects[m_region_intersect].m_dat, orig.m_dat, orig.m_n_prodsize * sizeof(double));
    return;
}

////////////////////////////////////////////////////////////////////////////////
// Code to read in factor graph and initialize MPLP.
////////////////////////////////////////////////////////////////////////////////

mplpLib::MPLPAlg::MPLPAlg(clock_t start, int time_limit, const std::string model_file, const std::string evid_file, FILE *log_file, bool uaiCompetition) : start(start), time_limit(time_limit), begin(false), m_best_val(-huge), last_obj(huge), total_mplp_iterations(0), previous_run_of_global_decoding(-1), obj_del(huge), _log_file(log_file), m_uaiCompetition(uaiCompetition) {

    size_t n;
    _res_fname = model_file.substr( (n = model_file.find_last_of('/')) == std::string::npos ? 0 : n + 1 ).append(".MPE");

    _ofs_res.open(_res_fname.c_str(), std::ios::out | std::ios::trunc);
    _ofs_res<<"MPE";
    _ofs_res<<std::endl;	       //first write the required "MPE" header
    _ofs_res.flush();

    // TODO: put close statement in deconstructor
    //	ofs_res.close();

    std::string tmp = model_file.substr(model_file.length()-6, 6);

    // We have two input formats.
    if (!tmp.compare("UAI.LG")) {
        Init2(model_file, evid_file);
    } else {
        Init(model_file, evid_file);
    }
}

void mplpLib::MPLPAlg::Init(const std::string fn, const std::string evid_fn){
    if(DEBUG_MODE)
        cout << "Calling Init() code...\n";

    std::vector<int> var_sizes;
    std::vector< std::vector<int> > all_factors;
    std::vector< std::vector<double> > all_lambdas;
    std::cout<<"Reading model file"<<std::endl;

    read_model_file(var_sizes, evidence, all_factors, all_lambdas, fn, evid_fn);

    std::cout<<"Initializing..."<<std::endl;
    Init(var_sizes, all_factors, all_lambdas);
}

void mplpLib::MPLPAlg::Init2(const std::string fn, const std::string evid_fn){
    if(DEBUG_MODE)
        cout << "Calling Init2() code...\n";

    std::vector<int> var_sizes;
    std::vector< std::vector<int> > all_factors;
    std::vector< std::vector<double> > all_lambdas;

    // NOTE: also accepts stereo files with special Potts potentials
    read_model_file2(var_sizes, evidence, all_factors, all_lambdas, fn, evid_fn);

    Init(var_sizes, all_factors, all_lambdas);
}


void mplpLib::MPLPAlg::Init(vector<int> & var_sizes, vector<vector<int> > & all_region_inds, vector<vector<double> > & all_lambdas) {
    // Set m_var_sizes
    m_var_sizes = var_sizes;   //invoking copy constructor

    // Set the intersection sets to be all single nodes and also all regions
    // Initialize sum into intersections.

    // First, add all individual nodes as their own intersection set
    for(int si=0; si < m_var_sizes.size(); ++si) {
        m_all_intersects.push_back(vector<int>(1,si));
        vector<int> subset_size(1,m_var_sizes[si]);    // Initialize sum into intersections to zero for these
        m_sum_into_intersects.push_back(MulDimArr(subset_size) = 0);
        m_single_node_lambdas.push_back(MulDimArr(subset_size) = 0);
    }

    // Next initialize all regions. If not a single node, give them their own intersection set

    for (int ri=0; ri < all_region_inds.size(); ++ri) {
        //vector<int> *region_var_sizes = new vector<int>();
        vector<int> region_var_sizes;
        for (int i=0; i < all_region_inds[ri].size(); ++i){
            region_var_sizes.push_back(m_var_sizes[all_region_inds[ri][i]]);
        }

        /*		if(DEBUG_MODE) {
		  if(all_region_inds[ri].size() == 0) {
		    cout << "Length zero region!" << endl;
		  }
		}*/

        if(all_region_inds[ri].size() == 1 && all_lambdas[ri].size() != 0){
            MulDimArr curr_lambda(region_var_sizes);
            for (int i = 0; i < curr_lambda.m_n_prodsize; ++i){
                curr_lambda[i] = all_lambdas[ri][i];
            }
            // Insert the single node potential into sum_into_intersects
            m_sum_into_intersects[all_region_inds[ri][0]] += curr_lambda;
            m_single_node_lambdas[all_region_inds[ri][0]] += curr_lambda;
        }else{
            vector<int> curr_intersect(all_region_inds[ri]);
            m_all_intersects.push_back(curr_intersect);
            int curr_intersect_loc = m_all_intersects.size() - 1;
            if (all_lambdas[ri].size()!=0) {
                MulDimArr curr_lambda = MulDimArr(region_var_sizes);

                // Assume all_lambdas is given as a "flat" vector and put it into curr_lambda
                for (int i=0; i < curr_lambda.m_n_prodsize; ++i){
                    curr_lambda[i] = all_lambdas[ri][i];
                }
                m_sum_into_intersects.push_back(curr_lambda);

                Region curr_region(all_region_inds[ri], m_all_intersects, all_region_inds[ri], m_var_sizes, curr_intersect_loc);

                m_all_regions.push_back(curr_region);
                m_region_lambdas.push_back(MulDimArr(curr_lambda));
            }else{  // Empty constructor
                Region curr_region(all_region_inds[ri], m_all_intersects, all_region_inds[ri], m_var_sizes, curr_intersect_loc);
                m_all_regions.push_back(curr_region);
                m_region_lambdas.push_back(MulDimArr());
            }
            // If this is an edge, insert into the map
            if(curr_intersect.size() == 2){
                // First sort
                vector<int> tmp_inds(curr_intersect);
                sort(tmp_inds.begin(), tmp_inds.end());
                // Then insert
                m_intersect_map.insert(pair<pair<int,int>,int>(pair<int,int>(tmp_inds[0], tmp_inds[1]), curr_intersect_loc));
            }
        }
    }

    // Initialize output vector
    for (int i=0; i<m_var_sizes.size(); ++i){
        m_decoded_res.push_back(0);
        m_best_decoded_res.push_back(0);
    }

    // Incorporate evidence
    for (std::map<int, int>::iterator it = evidence.begin(); it != evidence.end(); ++it){
        m_decoded_res[it->first] = it->second;
        m_best_decoded_res[it->first] = it->second;    //record evidence values
    }

    if(m_uaiCompetition) {

        // For UAI competition, for some silly reason the all 0's configuration is often a solution of the CSPs.
        // For papers one should keep this commented out, so as to not bias experimental results using our prior knowledge.
        UpdateResult(); // in case the MAP solution was actually all 0's

        // Recognize CSP instances by a dual value of identically zero and no fields
        // In these cases, randomly perturb the single node intersection sets
        double obj = LocalDecode();
        if(obj == 0) {

            // Check to see if fields
            bool isCSP = true;
            for(int si=0; si < m_var_sizes.size(); ++si) {
                for(int loc=0; loc < m_sum_into_intersects[si].m_n_prodsize; loc++) {
                    if(m_sum_into_intersects[si].m_dat[loc] != 0) {
                        isCSP = false;
                        break;
                    }
                }
                if(!isCSP) break;
            }

            if(isCSP) {
                if(DEBUG_MODE) cout << "Likely CSP instance. Randomly perturbing objective." << endl;
                CSP_instance = true;

                for(int si=0; si < m_var_sizes.size(); ++si) {
                    // Randomly perturb
                    for(int loc=0; loc < m_sum_into_intersects[si].m_n_prodsize; loc++) {
                        // TODO: how to set the scale?
                        m_sum_into_intersects[si].m_dat[loc] += .01 * rand() / double(RAND_MAX);
                    }
                }
            }
            else CSP_instance = false;
        }
    }

    last_global_decoding_end_time = 0;
    last_global_decoding_total_time = 0;
}

////////////////////////////////////////////////////////////////////////////////
// Main logic of MPLP (besides UpdateMsgs() which is in the Region class above)
////////////////////////////////////////////////////////////////////////////////

void mplpLib::MPLPAlg::RunMPLP(int niter, double obj_del_thr, double int_gap_thr){
    int ri;

    // Perform the GMPLP updates (Sontag's modified version), not quite as in the GJ NIPS07 paper
    for (int it=0; it<niter; ++it){

        for (ri=0; ri<m_all_regions.size(); ++ri){
            m_all_regions[ri].UpdateMsgs(m_sum_into_intersects);
        }

        total_mplp_iterations++;

        double obj, int_gap;

        obj = LocalDecode();
        obj_del = last_obj-obj;
        last_obj = obj;

        // Run global decoding at least once, a third of the way through
        if(previous_run_of_global_decoding == -1 &&
                ((double)(clock() - start) / CLOCKS_PER_SEC) > time_limit/3) {
            if(DEBUG_MODE)
                cout << "Third of the way! Going to run global decoding once." << endl;
            RunGlobalDecoding(false);
        }

        int_gap = obj - m_best_val;

        if (DEBUG_MODE){
            cout << "Iter=" << (it+1) << " Objective=" << obj <<  " Decoded=" << m_best_val << " ObjDel=" <<  obj_del << " IntGap=" << int_gap << endl;
        }
        if(_log_file != 0){
            fprintf(_log_file, "%.2f %.4f %.4f\n", ((double)(clock()-start)/CLOCKS_PER_SEC), obj, m_best_val);
        }

        if (obj_del<obj_del_thr && it > 16) // TODO: put these choices as parameters to the program
            break;
        if (int_gap<int_gap_thr)
            break;
    }
    return;
}

void mplpLib::MPLPAlg::AddAllEdgeIntersections()
{

    if(DEBUG_MODE) cout << "Adding all edge intersection sets..." << endl;

    // Iterate over all of the regions
    for (int ri=0; ri<m_all_regions.size(); ++ri){

        // We only care about the regions with >2 variables
        if(m_all_regions[ri].m_region_inds.size() <= 2) continue;

        // For each pair of the variables, add the corresponding intersection set
        for(int vi=0; vi < m_all_regions[ri].m_region_inds.size()-1; vi++) {
            int i = m_all_regions[ri].m_region_inds[vi];

            for(int vj= vi + 1; vj < m_all_regions[ri].m_region_inds.size(); vj++) {
                int j = m_all_regions[ri].m_region_inds[vj];

                // Find or add intersection set to the objective
                vector<int> ij_edge;
                ij_edge.push_back(i); ij_edge.push_back(j);
                int ij_intersect_loc = FindIntersectionSet(ij_edge);
                if(ij_intersect_loc == -1) {
                    ij_intersect_loc = AddIntersectionSet(ij_edge);
                    // If edge intersection set didn't exist, then there is no edge region either
                }

                // Add intersection set to the region
                // TODO: Test this more thoroughly.
                m_all_regions[ri].AddIntersectionSet(ij_intersect_loc, m_all_intersects, m_var_sizes);
            }
        }
    }

}

/*
 * Assumes that no intersection set already exists for this region (creates a new one).
 */
int mplpLib::MPLPAlg::AddRegion(vector<int> & inds_of_vars, vector<int> & intersect_inds)
{
    // No potential to go along with the region
    int region_intersection_set = AddIntersectionSet(inds_of_vars);
    Region new_region(inds_of_vars, m_all_intersects, intersect_inds, m_var_sizes, region_intersection_set);
    // This will also initialize the messages to zero, which is what we want
    m_all_regions.push_back(new_region);
    m_region_lambdas.push_back(MulDimArr());

    //  return m_all_regions.size()-1;
    return region_intersection_set;
}

int mplpLib::MPLPAlg::AddIntersectionSet(vector<int> & inds_of_vars)
{
    m_all_intersects.push_back(inds_of_vars);
    // Calculate the sizes of the variables in this set
    vector<int> sizes;
    for (int i=0; i< inds_of_vars.size(); ++i)
        sizes.push_back(m_var_sizes[inds_of_vars[i]]);

    // If this is an edge, insert into the map
    if(inds_of_vars.size() == 2)
    {
        // First sort
        vector<int> tmp_inds(inds_of_vars);
        sort(tmp_inds.begin(), tmp_inds.end());

        //		if(DEBUG_MODE)
        //		  cout << "adding edge intersection set " << tmp_inds[0] << " " << tmp_inds[1] << endl;

        // Then insert
        m_intersect_map.insert(pair<pair<int,int>,int>(pair<int,int>(tmp_inds[0], tmp_inds[1]), m_all_intersects.size()-1));
    }

    MulDimArr new_arr(sizes);
    new_arr = 0;
    m_sum_into_intersects.push_back(new_arr);
    return m_all_intersects.size()-1;
}

/*
 * Returns -1 if the intersection set not found.
 */
int mplpLib::MPLPAlg::FindIntersectionSet(vector<int> & inds_of_vars)
{
    // Sort the indices to make lookup and comparison easy
    vector<int> tmp_inds_of_vars(inds_of_vars);
    sort(tmp_inds_of_vars.begin(), tmp_inds_of_vars.end());

    // Is this an edge? If so, do map lookup
    if(tmp_inds_of_vars.size() == 2)
    {
        map<pair<int,int>, int>::iterator iter = m_intersect_map.find(pair<int,int>(tmp_inds_of_vars[0], tmp_inds_of_vars[1]));
        if( iter != m_intersect_map.end() ) // If the edge is found
            return iter->second;
        else {
            return -1;
        }
    }

    for(int i=0; i < m_all_intersects.size(); ++i)
    {
        // copy, then sort
        vector<int> tmp_inds(m_all_intersects[i]);
        sort(tmp_inds.begin(), tmp_inds.end());

        if(tmp_inds == tmp_inds_of_vars)
            return i;
    }

    return -1;
}

double mplpLib::MPLPAlg::IntVal(vector<int> & assignment) const{
    double int_val = 0;
    for (int ri=0; ri<m_all_regions.size(); ++ri){
        if (m_region_lambdas[ri].m_n_prodsize){
            vector<int> tmpvec;
            for (int vi = 0; vi < m_all_regions[ri].m_region_inds.size(); ++vi){
                tmpvec.push_back(assignment[m_all_regions[ri].m_region_inds[vi]]);
            }
            int_val+= m_region_lambdas[ri].GetVal(tmpvec);
        }
    }
    //This iterates over all singletons
    for (int ni=0; ni<m_var_sizes.size(); ++ni){
        int_val+=m_single_node_lambdas[ni][assignment[ni]];
    }
    return int_val;
}

double mplpLib::MPLPAlg::LocalDecode(void){
    double obj=0;
    int max_at;
    for (int si=0; si<m_sum_into_intersects.size(); ++si){
        obj+= m_sum_into_intersects[si].Max(max_at);
        // If this is a singleton, keep its value (so that we also have an integral assignment).
        // NOTE: Here we assume that all singletons are intersection sets. Otherwise, some variables will not be decoded here
        // evidence.find() is used because we do not want to set the state of a variable whose state
        // is fixed because it is evidence.
        if (m_all_intersects[si].size()==1 && evidence.find(m_all_intersects[si][0]) == evidence.end()){
            m_decoded_res[m_all_intersects[si][0]] = max_at;
        }
    }
    UpdateResult();
    return obj;
}

/*
 * Checks to see if the current integer assignment is better than any found before, and
 * if so writes it to the output file.
 */
double mplpLib::MPLPAlg::UpdateResult(void){
    double int_val;
    if ( (int_val = IntVal(m_decoded_res)) > m_best_val){

        if(DEBUG_MODE)
            cout << "int val: " << int_val << endl;

        if (time_limit + (double)(start - clock()) / CLOCKS_PER_SEC > MIN_APP_TIME){ // Prevent a partial write
            Write(_res_fname.c_str());
        }
        m_best_decoded_res.assign(m_decoded_res.begin(), m_decoded_res.end());
        m_best_val = int_val;
    }
    return int_val;
}

/*
 * Write to output file (containing best MAP assignment found so far).
 * TODO: Modify so that this writes a checkpoint, i.e. the full list of intersection sets, regions (not potentials),
 *       and messages, so we can use in debugging and re-running.
 */
void mplpLib::MPLPAlg::Write(const char *res_fname, const char *msgs_fname, const char *suminto_fname, const char *objhist_fname, const char *inthist_fname, const char *timehist_fname)
{
    stringstream s;

    if (begin){
        s << "-BEGIN-\n";
    }else{
        begin = true;
    }
    //"This will be the number of lines (not include this line) in the solution part."
    s << "1" << "\n" << m_var_sizes.size() << " ";

    for (int vi=0; vi< m_var_sizes.size()-1; ++vi){
        s << m_decoded_res[vi] << " ";
    }

    s << m_decoded_res[m_var_sizes.size()-1];

    //write solution in append mode and avoid writing a partial solution in case of timeout (SIGKILL)
    _ofs_res << s.str();
    _ofs_res << std::endl;
    _ofs_res.flush();
}


//////////////////////////////////////////////////////////////////////////////////////////
//
// The below code is not used in UAI '12 paper, but rather was used in the UAI
// competition. It attempts to find a good integer assignment from the current dual
// solution via several different decoding methods.
//
//////////////////////////////////////////////////////////////////////////////////////////


void mplpLib::MPLPAlg::RunGlobalDecoding(bool exhaustive){

    if(DEBUG_MODE) {
        cout << "Running global decoding..." << endl;
    }

    std::set<int> not_decoded;
    double global_decoding_start_time = (double)clock();

    int m, i, j, k;
    std::vector< std::vector< std::vector<double> > > tmp_msgs;
    std::vector< std::vector<double> > tmp_sums;
    std::map<int, int> tmp_evid = evidence;//, max_at;
    //std::map<int, double> gap_vals;
    int *max_at = new int[m_var_sizes.size()];
    double *gap_vals = new double [m_var_sizes.size()];

    int num_mplp_iters_global_decoding = 0;

    for (i = 0; i < m_all_regions.size(); ++i){
        tmp_msgs.push_back(std::vector< std::vector<double> >());
        for (j = 0; j < m_all_regions[i].m_msgs_from_region.size(); ++j){
            tmp_msgs[i].push_back(std::vector<double>());
            for (k = 0; k < m_all_regions[i].m_msgs_from_region[j].m_n_prodsize; ++k){
                tmp_msgs[i][j].push_back(m_all_regions[i].m_msgs_from_region[j].m_dat[k]);
            }
        }
    }

    for (i = 0; i < m_sum_into_intersects.size(); ++i){
        tmp_sums.push_back(std::vector<double>());
        for (j = 0; j < m_sum_into_intersects[i].m_n_prodsize; ++j){
            tmp_sums[i].push_back(m_sum_into_intersects[i].m_dat[j]);
        }
    }

    // Initialize all variables as not yet decoded
    for (i = 0; i < m_var_sizes.size(); ++i){
        if (evidence.find(i) == evidence.end()){
            not_decoded.insert(i);
        }
    }

    while (!not_decoded.empty()){
        double biggest_gap = -huge;
        for (std::set<int>::iterator s_it = not_decoded.begin(); s_it != not_decoded.end(); ++s_it){
            if( ( gap_vals[*s_it] = gap(*s_it, m) ) >= biggest_gap)
                biggest_gap = gap_vals[*s_it];
            max_at[*s_it] = m;
        }

        if(exhaustive)
            biggest_gap /= 1;
        else
            biggest_gap /= 10;

        if(DEBUG_MODE){
            cout << "biggest gap now set to " << biggest_gap << ", remaining nodes = " << not_decoded.size() << ", new int sol = " << m_best_val << endl;
        }

        // Iterative over nodes not yet decoded (TODO: use a reasonable gap criterion)
        for (std::set<int>::iterator s_it = not_decoded.begin(); s_it != not_decoded.end(); s_it++){
            if (gap_vals[*s_it] >= biggest_gap){   //gap stores argmax of reparametrized local potential in max_at

                evidence[*s_it] = max_at[*s_it];   //note: this is not permanent

                m_decoded_res[*s_it] = max_at[*s_it];
                not_decoded.erase(*s_it);

                for (i = 0; i < max_at[*s_it]; ++i){
                    m_sum_into_intersects[*s_it][i] = -huge;
                }
                while (++i < m_var_sizes[*s_it]){
                    m_sum_into_intersects[*s_it][i] = -huge;
                }

                if(exhaustive)
                    break;
            }
        }

        for (int it=0; it<10; ++it){
            for (int ri=0; ri<m_all_regions.size(); ++ri){
                m_all_regions[ri].UpdateMsgs(m_sum_into_intersects);
            }
            num_mplp_iters_global_decoding++;

            // Re-run local decoding, as some of the unfixed variables' assignments may have changed
            LocalDecode();
            UpdateResult();
        }

        // Give up after 100 rounds (1000 total MPLP iterations)
        //		if(!exhaustive && num_mplp_iters_global_decoding >= 1000)
        //		  break;
    }

    for (i = 0; i < m_all_regions.size(); ++i){
        for (j = 0; j < m_all_regions[i].m_msgs_from_region.size(); ++j){
            for (k = 0; k < m_all_regions[i].m_msgs_from_region[j].m_n_prodsize; ++k){
                m_all_regions[i].m_msgs_from_region[j].m_dat[k] = tmp_msgs[i][j][k];
            }
        }
    }
    for (i = 0; i < m_sum_into_intersects.size(); ++i){
        for (j = 0; j < m_sum_into_intersects[i].m_n_prodsize; ++j){
            m_sum_into_intersects[i].m_dat[j] = tmp_sums[i][j];
        }
    }
    evidence = tmp_evid;

    previous_run_of_global_decoding = total_mplp_iterations;
    last_global_decoding_end_time = (double)clock();
    last_global_decoding_total_time = (last_global_decoding_end_time - global_decoding_start_time)/CLOCKS_PER_SEC;

    delete [] max_at;
    delete [] gap_vals;
}

//note: sorting gap values in ascending order and fixing the node with lowest gap value first (to resolve possible frustration)
// TODO: put some randomness into this! Of all tied states, randomly sample. Or, randomly permute tied nodes, then choose 1. Or both.
void mplpLib::MPLPAlg::RunGlobalDecoding2(bool exhaustive){
    if(DEBUG_MODE) {
        cout << "Running global decoding2..." << endl;
    }

    std::set<int> not_decoded;
    double global_decoding_start_time = (double)clock();

    int m, i, j;
    std::vector< std::vector< std::vector<double> > > tmp_msgs;
    std::vector< std::vector<double> > tmp_sums;
    std::map<int, int> tmp_evid = evidence;//, max_at;
    int *max_at = new int[m_var_sizes.size()];
    double *gap_vals = new double [m_var_sizes.size()];

    int num_mplp_iters_global_decoding = 0;

    //saving current mplp state
    for (i = 0; i < m_all_regions.size(); ++i){
        tmp_msgs.push_back(std::vector< std::vector<double> >());
        for (j = 0; j < m_all_regions[i].m_msgs_from_region.size(); ++j){
            tmp_msgs[i].push_back(std::vector<double>());
            for (int k = 0; k < m_all_regions[i].m_msgs_from_region[j].m_n_prodsize; ++k){
                tmp_msgs[i][j].push_back(m_all_regions[i].m_msgs_from_region[j].m_dat[k]);
            }
        }
    }

    for (i = 0; i < m_sum_into_intersects.size(); ++i){
        tmp_sums.push_back(std::vector<double>());
        for (j = 0; j < m_sum_into_intersects[i].m_n_prodsize; ++j){
            tmp_sums[i].push_back(m_sum_into_intersects[i].m_dat[j]);
        }
    }

    // find all variables that are not yet decoded
    for (int i = 0; i < m_var_sizes.size(); ++i){
        if (evidence.find(i) == evidence.end()){
            not_decoded.insert(i);
        }
    }

    while (!not_decoded.empty()){
        int i;
        double smallest_gap = huge;
        int index_smallest = -1;
        for (std::set<int>::iterator s_it = not_decoded.begin(); s_it != not_decoded.end(); ++s_it){
            if( ( gap_vals[*s_it] = gap(*s_it, m) ) < smallest_gap) {
                smallest_gap = gap_vals[*s_it];
                index_smallest = *s_it;
            }else if (gap_vals[*s_it] == smallest_gap){
                if (rand() % 2){   //some randomness for tied nodes (but this may still not be uniform)
                    if(DEBUG_MODE)
                        cout << "Adding some randomness to globaldecoding2" << endl;
                    index_smallest = *s_it;
                }
            }
            max_at[*s_it] = m;
        }
        // Stopping criterion
        if(!exhaustive && smallest_gap > GAP_THR)
            break;

        if(DEBUG_MODE)
            cout << "Fixing node " << index_smallest << " to value " << max_at[index_smallest] << endl;

        // Fix one at a time
        evidence[index_smallest] = max_at[index_smallest];   //note: this is not permanent
        m_decoded_res[index_smallest] = max_at[index_smallest];
        not_decoded.erase(index_smallest);

        for (i = 0; i < max_at[index_smallest]; ++i){
            m_sum_into_intersects[index_smallest][i] = -huge;
        }
        while (++i < m_var_sizes[index_smallest]){
            m_sum_into_intersects[index_smallest][i] = -huge;
        }

        for (int it=0; it<10; ++it){
            for (int ri=0; ri<m_all_regions.size(); ++ri){
                m_all_regions[ri].UpdateMsgs(m_sum_into_intersects);
            }
            num_mplp_iters_global_decoding++;

            // Re-run local decoding, as some of the unfixed variables' assignments may have changed
            LocalDecode();
            UpdateResult();
        }

        //		if(DEBUG_MODE){
        //		  cout << "new int sol = " << m_best_val << endl;
        //		}

        // Give up after 100 rounds (1000 total MPLP iterations)
        if(!exhaustive && num_mplp_iters_global_decoding >= 1000) {

            if(DEBUG_MODE)
                cout << "Exhausted number of rounds for global decoding. Quitting." << endl;
            break;
        }
    }
    for (int i = 0; i < m_all_regions.size(); ++i){
        for (int j = 0; j < m_all_regions[i].m_msgs_from_region.size(); ++j){
            for (int k = 0; k < m_all_regions[i].m_msgs_from_region[j].m_n_prodsize; ++k){
                m_all_regions[i].m_msgs_from_region[j].m_dat[k] = tmp_msgs[i][j][k];
            }
        }
    }
    for (i = 0; i < m_sum_into_intersects.size(); ++i){
        for (j = 0; j < m_sum_into_intersects[i].m_n_prodsize; ++j){
            m_sum_into_intersects[i].m_dat[j] = tmp_sums[i][j];
        }
    }
    evidence = tmp_evid;

    previous_run_of_global_decoding = total_mplp_iterations;
    last_global_decoding_end_time = (double)clock();
    last_global_decoding_total_time = (last_global_decoding_end_time - global_decoding_start_time)/CLOCKS_PER_SEC;

    delete [] max_at;
    delete [] gap_vals;
}


// Do large numbers of random objective permutations, run 10 iterations of MPLP, restore
void mplpLib::MPLPAlg::RunGlobalDecoding3(void){

    if(DEBUG_MODE) {
        cout << "Running global decoding3..." << endl;
    }

    //saving current mplp state
    std::vector< std::vector< std::vector<double> > > tmp_msgs;
    std::vector< std::vector<double> > tmp_sums;
    for (int i = 0; i < m_all_regions.size(); ++i){
        tmp_msgs.push_back(std::vector< std::vector<double> >());
        for (int j = 0; j < m_all_regions[i].m_msgs_from_region.size(); ++j){
            tmp_msgs[i].push_back(std::vector<double>());
            for (int k = 0; k < m_all_regions[i].m_msgs_from_region[j].m_n_prodsize; ++k){
                tmp_msgs[i][j].push_back(m_all_regions[i].m_msgs_from_region[j].m_dat[k]);
            }
        }
    }
    for (int i = 0; i < m_sum_into_intersects.size(); ++i){
        tmp_sums.push_back(std::vector<double>());
        for (int j = 0; j < m_sum_into_intersects[i].m_n_prodsize; ++j){
            tmp_sums[i].push_back(m_sum_into_intersects[i].m_dat[j]);
        }
    }

    // Do large numbers of random objective permutations, run 4 iterations of MPLP, restore
    for(int trial=0; trial <= 10; trial++) {

        for(int si=0; si < m_var_sizes.size(); ++si) {
            // Randomly perturb single node potentials
            for(int loc=0; loc < m_sum_into_intersects[si].m_n_prodsize; loc++) {
                // TODO: how to set the scale?  perhaps look at objective value.
                // For GRIDS, 10*.01 works well
                //   -- does this just come from the randomness in decoding?
                m_sum_into_intersects[si].m_dat[loc] += 10*.01 * rand() / double(RAND_MAX);
            }
        }

        bool exhaustive = false;
        if(trial == 10)
            exhaustive = true;

        RunGlobalDecoding(exhaustive);
        RunGlobalDecoding2(exhaustive);

        // Restoring MPLP state
        for (int i = 0; i < m_all_regions.size(); ++i){
            for (int j = 0; j < m_all_regions[i].m_msgs_from_region.size(); ++j){
                for (int k = 0; k < m_all_regions[i].m_msgs_from_region[j].m_n_prodsize; ++k){
                    m_all_regions[i].m_msgs_from_region[j].m_dat[k] = tmp_msgs[i][j][k];
                }
            }
        }
        for (int i = 0; i < m_sum_into_intersects.size(); ++i){
            for (int j = 0; j < m_sum_into_intersects[i].m_n_prodsize; ++j){
                m_sum_into_intersects[i].m_dat[j] = tmp_sums[i][j];
            }
        }
    }
}


bool mplpLib::MPLPAlg::RunDecimation(void){
    bool fixed_node = false;

    int i, m;
    double biggest_gap = -huge;
    int *max_at = new int[m_var_sizes.size()];
    double *gap_vals = new double[m_var_sizes.size()];
    std::vector<int> not_decoded;
    if(DEBUG_MODE) {
        cout << "Running decimation..." << endl;
    }
    for (int i = 0; i < m_var_sizes.size(); ++i){
        if (evidence.find(i) == evidence.end()){
            not_decoded.push_back(i);
        }
    }
    for (std::vector<int>::iterator s_it = not_decoded.begin(); s_it != not_decoded.end(); ++s_it){
        if( ( gap_vals[*s_it] = gap(*s_it, m) ) >= biggest_gap){
            biggest_gap = gap_vals[*s_it];
        }
        if(DEBUG_MODE) {
            //		  cout << "gap " << gap_vals[*s_it] << endl;
        }
        max_at[*s_it] = m;
    }

    biggest_gap /= 1;
    for (std::vector<int>::iterator s_it = not_decoded.begin(); s_it != not_decoded.end(); ++s_it){
        if (gap_vals[*s_it] >= biggest_gap){   //gap stores argmax of reparametrized local potential in max_at

            //		  if(DEBUG_MODE) {
            //		    cout << "Fixing node " << *s_it << endl;
            //		  }
            fixed_node = true;

            evidence[*s_it] = max_at[*s_it];   //note: this is permanently fixed for the current instance of MPLP
            m_decoded_res[*s_it] = max_at[*s_it];
            for (i = 0; i < max_at[*s_it]; ++i){
                m_sum_into_intersects[*s_it][i] = -huge;
            }
            while (++i < m_var_sizes[*s_it]){
                m_sum_into_intersects[*s_it][i] = -huge;
            }
        }
    }
    delete [] max_at;
    delete [] gap_vals;

    return fixed_node;
}

double mplpLib::MPLPAlg::gap(int n, int & max_at) const{
    if (m_var_sizes[n] < 2){
        max_at = 0;
        return huge;
    }
    double m, m1 = -huge, m2 = -huge;
    for (int s = 0; s < m_var_sizes[n]; ++s){
        if ((m = m_sum_into_intersects[n][s]) > m1){
            m2 = m1;
            m1 = m;
            max_at = s;
        }else if (m > m2){
            m2 = m;
        }
    }
    return m1 - m2;
}

