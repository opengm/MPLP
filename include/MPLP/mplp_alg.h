/*
 *  mplp_alg.h
 *  mplp
 *
 *  Created by Amir Globerson and David Sontag on 6/25/08.
 *  Copyright 2008 MIT. All rights reserved.
 *
 */
#ifndef MPLP_MPLP_ALG_H
#define MPLP_MPLP_ALG_H

#include <vector>
#include <iostream>
#include <fstream>
#include <set>
#include <map>

#include <MPLP/muldim_arr.h>

#include <MPLP/read_model_file.h>

#define MIN_APP_TIME .0001  //amount of time reserved for appending an answer into a file (to prevent any partially written answers)

class Region
{
public:
    // The variables in the region
    std::vector<int> m_region_inds;

    // Specifies the indices of the intersection sets
    std::vector<int> m_intersect_inds;

    // Specifies the indices corresponding to the position of each
    // intersection set in the region
    std::vector<std::vector<int> > m_inds_of_intersects;

    // Every region has a corresponding intersection set. What is its index?
    int m_region_intersect;

    // Contains the messages from each region to its intersection sets
    std::vector<MulDimArr> m_msgs_from_region;
    std::vector<int> m_var_sizes;

    Region(std::vector<int> & region_inds, std::vector<std::vector<int> > & all_intersects, std::vector<int> & intersect_inds, std::vector<int> & var_sizes, int region_intersect);

    // Adds intersection set to the region
    void AddIntersectionSet(int intersect_loc, std::vector<std::vector<int> > & all_intersects, std::vector<int> & var_sizes);

    void UpdateMsgs(std::vector<MulDimArr> & sum_into_intersects);
    int Get_nVars() {return m_var_sizes.size();};
};

class MPLPAlg
{
public:

    /* Invariant: every region has a corresponding intersection set, and all
     * potentials are inside of the intersection sets. The objective function
     * can be computed simply by summing over the intersection sets.
     * "Region" simply corresponds to an algorithmic concept, specifying both
     * dual variables shared between a larger intersection set and smaller ones,
     * and an update strategy. See Sontag's Ph.D. thesis page 104.
     */
    bool begin;    //whether the current solution is the first 1

    double m_best_val, last_obj, obj_del;   //the best primal objective so far
    int total_mplp_iterations;
    int previous_run_of_global_decoding;
    double last_global_decoding_end_time;
    double last_global_decoding_total_time;

    bool m_uaiCompetition;
    std::vector<std::vector<int> > m_all_intersects;
    std::vector<Region> m_all_regions;
    //	vector<vector<int> > m_all_region_inds;
    std::vector<MulDimArr> m_sum_into_intersects;
    //	vector<MulDimArr> m_all_lambdas;
    //	vector<vector<int> > m_all_region_intersects;
    std::map<int, int> evidence;
    std::vector<int> m_var_sizes;
    std::vector<int> m_decoded_res;
    std::vector<int> m_best_decoded_res;
    std::vector<double> m_objhist;
    std::vector<double> m_inthist;
    std::vector<double> m_timehist;
    std::vector<MulDimArr> m_single_node_lambdas;
    std::vector<MulDimArr> m_region_lambdas;

    // Is this a CSP instance?
    bool CSP_instance;

    // This map allows us to quickly look up the index of edge intersection sets
    std::map<std::pair<int, int>, int> m_intersect_map;

    //create an MPLP instance from the model file and evidence file (if any)
    MPLPAlg(clock_t, int, const std::string, const std::string, FILE *, bool uaiCompetition);

    MPLPAlg(void){};     //for decoding purpose only

    void Init(const std::string, const std::string = "");

    void Init2(const std::string, const std::string = "");

    void Init(std::vector<int> & var_sizes, std::vector<std::vector<int> > & all_region_inds, std::vector<std::vector<double> > & all_lambdas);

    void RunMPLP(int, double, double);

    double IntVal(std::vector<int> & assignment) const;

    double gap(int, int &) const;

    // argument specifies whether to go one by one (true) or to do in blocks (false), which is faster
    void RunGlobalDecoding(bool);
    void RunGlobalDecoding2(bool);
    void RunGlobalDecoding3(void);

    // returns false if no more variables are left to decimate
    bool RunDecimation(void);

    // Add a new region and return its index. intersect_inds refers to the index of the intersection sets that this
    // region intersects with (that is, the index into m_all_intersects)
    int AddRegion(std::vector<int> & inds_of_vars, std::vector<int> & intersect_inds);
    // As in the Matlab code, for now we will assume that an intersetion set is added before adding the regions that
    // intersect with it
    int AddIntersectionSet(std::vector<int> & inds_of_vars);

    // For regions of size >2, remove single node intersection sets and add all edge intersection sets
    void AddAllEdgeIntersections();

    // Find the index number into m_all_intersects of a given set of variables' intersection set.
    // Returns -1 if not found.
    int FindIntersectionSet(std::vector<int> & inds_of_vars);

    void Write(const char *res_fname, const char *msgs_fname = "msgs.txt", const char *suminto_fname = "suminto.txt", const char *objhist_fname = "objhist.txt", const char *inthist_fname = "inthist.txt", const char *timehist_fname = "timehist.txt");
private:
    std::string _res_fname;
    std::ofstream _ofs_res;
    FILE *_log_file; // TODO: make sure not to use when not initialized
    std::ifstream rnd_seed;
    clock_t start;
    int time_limit;
    double LocalDecode(void);      //single node decoding
    double UpdateResult(void);   //returns primal objective of this mplp instance
};

#endif
