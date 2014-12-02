/*
 *
 *  Created by Amir Globerson and David Sontag on 8/10/08.
 *  Updated by David Sontag, Do Kook Choe, and Yitao Li in 2012.
 *  Copyright 2008 MIT, 2012 NYU. All rights reserved.
 *
 */
#ifndef MPLP_CYCLE_H
#define MPLP_CYCLE_H

#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <set>
#include <list>
#include <map>
#include <queue>

#include <MPLP/mplp_config.h>
#include <MPLP/mplp_alg.h>

namespace mplpLib {

#define MPLP_DEBUG_MODE 0

#define MPLP_Inf 9999999999.9
#define MPLP_CYCLE_THRESH .00001  // TODO: Figure out how to set this
#define MPLP_CLUSTER_THR .0000001

typedef std::map<std::pair<MPLPIndexType, MPLPIndexType>, MPLPIndexType> mapType;
typedef std::vector<std::vector<std::pair<MPLPIndexType, double> > > adj_type;

// Structure for storing a candidate triplet cluster for tightening
struct TripletCluster {
    double bound;
    MPLPIndexType i,j,k;
    MPLPIndexType ij_intersect_loc, jk_intersect_loc, ki_intersect_loc;

    bool operator <(const TripletCluster & rhs) const {
        return bound < rhs.bound;
    }
};


/////////////////////////////////////////////////////////////////////////////////
// Code for union-find data structure (used by FindPartition)

struct Node { //node for union-find
    MPLPIndexType bit; //10, 01, or 11
    MPLPIndexType rank;
    Node* parent;
    MPLPIndexType i_size; //number of elements in i
    MPLPIndexType j_size; //number of elements in j

    Node(MPLPIndexType b) : bit(b), rank(0), parent(this) {
        if (b == 2) { //represent elements of i
            i_size = 1;
            j_size = 0;
        }
        else { //reprsent elements of j
            i_size = 0;
            j_size = 1;
        }
    }
};

Node* find(Node* n) {
    if (n != n->parent) {
        n->parent = find(n->parent);
    }
    return n->parent;
} 


Node* merge(Node* x, Node* y) {
    Node* root_x = find(x);
    Node* root_y = find(y);
    if (root_x == root_y) { //x and y have the same head
        return root_x;
    }

    if (root_x->rank > root_y->rank) {
        root_y->parent = root_x; //head x
        root_x->i_size += root_y->i_size; //add number of elements of i to the head node
        root_x->j_size += root_y->j_size; //add number of elements of j to the head node
        return root_x;
    }
    else {
        root_x->parent = root_y; //head y
        if (root_x->rank == root_y->rank) {
            ++(root_y->rank);
        }
        root_y->i_size += root_x->i_size; //add number of elements of i to the head node
        root_y->j_size += root_x->j_size; //add number of elements of j to the head node
        return root_y;
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Code to evaluate how good a cycle cluster is

double maximizeIndependently(std::vector<MulDimArr*> & beliefs)
{
    double sum=0.0;
    MPLPIndexType max_at; // not actually needed
    for(MPLPIndexType i=0; i < beliefs.size(); i++)
    {
        sum += beliefs[i]->Max(max_at);
    }

    return sum;
}

double getValCycle(std::vector<MulDimArr*> & beliefs, std::vector<bool> & b_transpose, std::vector<MPLPIndexType> & assignments)
{
    double sum=0.0;
    //std::vector<MPLPIndexType> inds; inds.push_back(-1); inds.push_back(-1); // temp
    std::vector<MPLPIndexType> inds(2);

    // All except the last edge
    assert(beliefs.size() > 0);
    for(MPLPIndexType i=0; i < beliefs.size()-1; i++)
    {
        inds[b_transpose[i]?1:0] = assignments[i];
        inds[b_transpose[i]?0:1] = assignments[i+1];
        sum += beliefs[i]->GetVal(inds);
    }

    // Now do last edge
    inds[b_transpose[beliefs.size()-1]?1:0] = assignments[beliefs.size()-1];
    inds[b_transpose[beliefs.size()-1]?0:1] = assignments[0];
    sum += beliefs[beliefs.size()-1]->GetVal(inds);

    return sum;
}

double maximizeCycle(std::vector<MulDimArr*> & beliefs, std::vector<bool> & b_transpose)
{
    double max_val = -MPLP_Inf;

    // Fix value of the first variable
    MPLPIndexType first_var_size = beliefs[0]->m_base_sizes[b_transpose[0]?1:0];
    MPLPIndexType second_var_size = beliefs[0]->m_base_sizes[b_transpose[0]?0:1];
    for(MPLPIndexType vo=0; vo < first_var_size; vo++)
    {
        //std::vector<MPLPIndexType> inds; inds.push_back(-1); inds.push_back(-1); // temp
        std::vector<MPLPIndexType> inds(2);
        inds[b_transpose[0]?1:0] = vo;

        // Do first edge (construct initial field)
        std::vector<double> field;
        for(MPLPIndexType v2=0; v2 < second_var_size; v2++)
        {
            inds[b_transpose[0]?0:1] = v2;
            field.push_back(beliefs[0]->GetVal(inds));
        }

        // Go over rest of edges, except last (which has to be treated specially)
        assert(beliefs.size() > 0);
        for(MPLPIndexType i=1; i < beliefs.size()-1; i++)
        {
            std::vector<double> new_field;
            for(MPLPIndexType v2=0; v2 < beliefs[i]->m_base_sizes[b_transpose[i]?0:1]; v2++)
            {
                inds.clear(); inds.push_back(-1); inds.push_back(-1); // temp
                inds[b_transpose[i]?0:1] = v2;

                // Take max
                double tmp_max_val = -MPLP_Inf;
                for(MPLPIndexType v1=0; v1 < field.size(); v1++)
                {
                    inds[b_transpose[i]?1:0] = v1;
                    tmp_max_val = std::max(tmp_max_val, field[v1]+beliefs[i]->GetVal(inds));
                }
                new_field.push_back(tmp_max_val);
            }
            field.clear(); // necessary?
            field = new_field;
        }

        // Do last edge (fix endpoint value to vo)
        //inds.clear(); inds.push_back(-1); inds.push_back(-1); // temp
        inds[b_transpose[b_transpose.size()-1]?0:1] = vo;

        // Take max
        double tmp_max_val = -MPLP_Inf;
        for(MPLPIndexType v1=0; v1 < field.size(); v1++)
        {
            inds[b_transpose[b_transpose.size()-1]?1:0] = v1;
            tmp_max_val = std::max(tmp_max_val, field[v1]+beliefs[beliefs.size()-1]->GetVal(inds));
        }

        max_val = std::max(max_val, tmp_max_val);
    }
    return max_val;
}

/////////////////////////////////////////////////////////////////////////////////
// Implementation of UAI 2008 algorithm (just for triplets; square functionality removed)

MPLPIndexType TightenTriplet(MPLPAlg & mplp, MPLPIndexType nclus_to_add_min, MPLPIndexType nclus_to_add_max, std::map<std::vector<MPLPIndexType>, bool >& triplet_set, double & promised_bound) {

    MPLPIndexType nClustersAdded = 0;
    MPLPIndexType nNewClusters = 0;

    if(MPLP_DEBUG_MODE)
        std::cout << "Doing pre-processing for adding triplet clusters." << std::endl;

    // Initialize adjacency list (filled in later) TODO: only do this when needed
    std::vector<MPLPIndexType>* adjacency_list = new std::vector<MPLPIndexType>[mplp.m_var_sizes.size()];

    // Construct adjacency list for the graph
    // Iterate over all of the edges (we do this by looking at the edge intersection sets)
    for(mapType::const_iterator it = mplp.m_intersect_map.begin(); it != mplp.m_intersect_map.end(); ++it)
    {
        // Get the two nodes i & j
        MPLPIndexType i=it->first.first; MPLPIndexType j=it->first.second;
        adjacency_list[i].push_back(j);
        adjacency_list[j].push_back(i);
    }

    // Sort the adjacency list, for fast intersections later
    for(MPLPIndexType i=0; i < sizeof(adjacency_list)/sizeof(std::vector<MPLPIndexType>); i++)
    {
        sort(adjacency_list[i].begin(), adjacency_list[i].end());
    }

    // Count the number of triangles
    std::vector<MPLPIndexType>::iterator intersects_iter_end;
    std::vector<MPLPIndexType> commonNodes(mplp.m_var_sizes.size());
    for(mapType::const_iterator it = mplp.m_intersect_map.begin(); it != mplp.m_intersect_map.end(); ++it)
    {

        // Get the two nodes i & j
        MPLPIndexType i=it->first.first; MPLPIndexType j=it->first.second;

        // Now find all neighbors of both i and j to see where the triangles are
        intersects_iter_end = set_intersection(adjacency_list[i].begin(), adjacency_list[i].end(), adjacency_list[j].begin(), adjacency_list[j].end(), commonNodes.begin());

        for(std::vector<MPLPIndexType>::const_iterator n=commonNodes.begin(); n != intersects_iter_end; ++n)
        {
            // Since a triplet shows up three times as an edge plus
            // a node, we only consider it for the case when n<i and n<j
            if(*n < i && *n < j)
                nNewClusters++;
        }
    }

    if(nNewClusters == 0) {
        if(MPLP_DEBUG_MODE)
            std::cout << "nNewClusters = 0. Returning." << std::endl;

        delete []adjacency_list;
        return 0;
    }

    if(MPLP_DEBUG_MODE)
        std::cout << "Looking for triangle clusters to add (" << nNewClusters << " triplets) " << std::endl;

    // TODO: put this elsewhere so that the space isn't re-allocated continuously?
    // Enumerate over all of the edges
    TripletCluster* newCluster = new TripletCluster[nNewClusters];

    MPLPIndexType index=0;

    // Iterate over all of the edge intersection sets
    std::vector<MPLPIndexType> tripAssignment; tripAssignment.push_back(-1); tripAssignment.push_back(-1); tripAssignment.push_back(-1);
    for(mapType::const_iterator it = mplp.m_intersect_map.begin(); it != mplp.m_intersect_map.end(); ++it)
    {

        // Get the two nodes i & j, and the edge intersection index
        MPLPIndexType i=it->first.first; MPLPIndexType j=it->first.second;
        MPLPIndexType ij_intersect_loc = it->second;

        // Now find all neighbors of both i and j to see where the triangles are
        // TEMP TEMP -- fails at i=0, j=1, on i==3.
        intersects_iter_end = set_intersection(adjacency_list[i].begin(), adjacency_list[i].end(), adjacency_list[j].begin(), adjacency_list[j].end(), commonNodes.begin());

        for(std::vector<MPLPIndexType>::const_iterator n=commonNodes.begin(); n != intersects_iter_end; ++n)
        {
            MPLPIndexType k = *n;

            // Since a triplet shows up three times as an edge plus
            // a node, we only consider it for the case when k<i and k<j
            if(!(k < i && k < j))
                continue;

            newCluster[index].i = i;
            newCluster[index].j = j;
            newCluster[index].k = k;

            // Find the intersection sets for this triangle
            newCluster[index].ij_intersect_loc = ij_intersect_loc;
            std::vector<MPLPIndexType> jk_edge; jk_edge.push_back(newCluster[index].j); jk_edge.push_back(newCluster[index].k);
            newCluster[index].jk_intersect_loc = mplp.FindIntersectionSet(jk_edge);
            std::vector<MPLPIndexType> ki_edge; ki_edge.push_back(newCluster[index].k); ki_edge.push_back(newCluster[index].i);
            newCluster[index].ki_intersect_loc = mplp.FindIntersectionSet(ki_edge);

            // Construct the beliefs for each edge, which will be maximized below
            std::vector<MulDimArr*> beliefs;
            std::vector<bool> b_transpose;
            beliefs.push_back(&mplp.m_sum_into_intersects[newCluster[index].ij_intersect_loc]);
            b_transpose.push_back(mplp.m_all_intersects[newCluster[index].ij_intersect_loc][0] != newCluster[index].i); // i first
            beliefs.push_back(&mplp.m_sum_into_intersects[newCluster[index].jk_intersect_loc]);
            b_transpose.push_back(mplp.m_all_intersects[newCluster[index].jk_intersect_loc][0] != newCluster[index].j); // then 'j'
            beliefs.push_back(&mplp.m_sum_into_intersects[newCluster[index].ki_intersect_loc]);
            b_transpose.push_back(mplp.m_all_intersects[newCluster[index].ki_intersect_loc][0] != newCluster[index].k); // then 'k'

            double bound_indep = maximizeIndependently(beliefs);

            // Before doing expensive joint maximization, see if we can quickly find optimal assignment
            tripAssignment[0] = mplp.m_decoded_res[i]; tripAssignment[1] = mplp.m_decoded_res[j]; tripAssignment[2] = mplp.m_decoded_res[k];
            double bound_quick = getValCycle(beliefs, b_transpose, tripAssignment);
            if(bound_indep == bound_quick)
            {
                newCluster[index].bound = 0;
            } else {
                // Do expensive joint maximization
                newCluster[index].bound = bound_indep - maximizeCycle(beliefs, b_transpose);
            }

            index++;
        }
    }

    // TODO opt: have a class for a cluster, so we can have different types and sort by bound,
    //       choosing best one by bound.
    //       Make the sorting and adding independent of the type of graph...

    // Sort the clusters by the bound
    std::sort(newCluster, newCluster+nNewClusters);

    if(MPLP_DEBUG_MODE)
        std::cout << " -- Considered " << nNewClusters << " clusters, smallest bound " << newCluster[std::max(static_cast<int>(nNewClusters)-static_cast<int>(nclus_to_add_max), 0)].bound << ", largest bound " << newCluster[nNewClusters-1].bound << std::endl;

    promised_bound = newCluster[nNewClusters-1].bound;

    // Add the top nclus_to_add clusters to the relaxation
    assert(nNewClusters > 0);
    for(MPLPIndexType clusterId = nNewClusters-1; clusterId >= 0 && nClustersAdded < nclus_to_add_max && (nClustersAdded < nclus_to_add_min || ((newCluster[clusterId].bound >= newCluster[nNewClusters-1].bound/5) && newCluster[clusterId].bound >= MPLP_CLUSTER_THR)) ; clusterId--)
    {
        // Check to see if this triplet is already being used
        std::vector<MPLPIndexType> temp;
        temp.push_back(newCluster[clusterId].i);
        temp.push_back(newCluster[clusterId].j);
        temp.push_back(newCluster[clusterId].k);
        sort(temp.begin(), temp.end());

        std::map<std::vector<MPLPIndexType>, bool >::iterator t_itr = triplet_set.find(temp);
        if (t_itr == triplet_set.end())
            triplet_set.insert(std::pair<std::vector<MPLPIndexType>, bool >(temp, true));
        else {
            //	if(MPLP_DEBUG_MODE)
            //	  cout << "   Triplet was already present. Skipping..." << endl;
            continue;
        }

        // Now add cluster ijk
        std::vector<MPLPIndexType> ijk_inds;
        ijk_inds.push_back(newCluster[clusterId].i); ijk_inds.push_back(newCluster[clusterId].j); ijk_inds.push_back(newCluster[clusterId].k);

        std::vector<MPLPIndexType> ijk_intersect_inds;
        ijk_intersect_inds.push_back(newCluster[clusterId].ij_intersect_loc);
        ijk_intersect_inds.push_back(newCluster[clusterId].jk_intersect_loc);
        ijk_intersect_inds.push_back(newCluster[clusterId].ki_intersect_loc);

        mplp.AddRegion(ijk_inds, ijk_intersect_inds);

        if(MPLP_DEBUG_MODE)
            std::cout << "Cluster added on nodes " << newCluster[clusterId].i << ", " << newCluster[clusterId].j << ", " << newCluster[clusterId].k << std::endl;
        // TODO: log which clusters are chosen...

        nClustersAdded++;
    }

    delete []newCluster;
    delete []adjacency_list;

    return nClustersAdded;
}


/////////////////////////////////////////////////////////////////////////////////
// Everything that follows is for the UAI 2012 cycle finding algorithm which
// can find arbitrary-length cycles to use in tightening the relaxation.

bool edge_sorting(std::list<MPLPIndexType> i, std::list<MPLPIndexType> j) {
    return i.back() > j.back();
}

bool edge_sorting3(std::pair<std::pair<MPLPIndexType, MPLPIndexType>, double> i, std::pair<std::pair<MPLPIndexType, MPLPIndexType>, double> j) {
    return i.second > j.second;
}

// De-allocate memory relating to the projection graph (TODO: finish writing this)
void delete_projection_graph(MPLPIndexType num_vars, std::vector<std::vector<MPLPIndexType> > &projection_map, std::vector<MPLPIndexType> &projection_imap_var, adj_type &projection_adjacency_list, double* &array_of_sij) {

    //for(MPLPIndexType node=0; node < num_vars; node++)
    //  delete []projection_map[node];
    //delete []projection_map;

    //delete projection_imap_var;

    //delete []projection_adjacency_list;

    delete []array_of_sij;
}

// Compute a single edge weight in the projection graph
double find_smn(const std::vector<MPLPIndexType>& partition_i, MPLPIndexType var_i_size, const std::vector<MPLPIndexType>& partition_j, MPLPIndexType var_j_size, MulDimArr* edge_belief) {

    MPLPIndexType* whole_i = new MPLPIndexType[var_i_size];
    MPLPIndexType* whole_j = new MPLPIndexType[var_j_size];
    for (MPLPIndexType i = 0; i < var_i_size; i++) {
        whole_i[i] = 0;
    }
    for (MPLPIndexType i = 0; i < var_j_size; i++) {
        whole_j[i] = 0;
    }

    double smn = -MPLP_Inf;
    for (MPLPIndexType i = 0; i < partition_i.size(); i++) {
        whole_i[partition_i[i]] = 1;
        for (MPLPIndexType j = 0; j < partition_j.size(); j++) {
            whole_j[partition_j[j]] = 1;
        }
    }

    double sec_max = -MPLP_Inf;
    for (MPLPIndexType i = 0; i < var_i_size; i++) {
        for (MPLPIndexType j = 0; j < var_j_size; j++) {
            if (whole_i[i] == whole_j[j]) {
                std::vector<MPLPIndexType> inds; inds.push_back(i); inds.push_back(j);
                double temp_val = edge_belief->GetVal(inds);
                if (smn < temp_val) smn = temp_val;
            }
            else {
                std::vector<MPLPIndexType> inds; inds.push_back(i); inds.push_back(j);
                double temp_val = edge_belief->GetVal(inds);
                if (sec_max < temp_val) sec_max = temp_val;
            }
        }
    }
    smn -= sec_max;

    delete [] whole_i;
    delete [] whole_j;

    return smn;
}

// Compute a single edge weight in the projection graph (more efficiently)
double find_smn_state_i(MPLPIndexType single_i, MPLPIndexType var_i_size, const std::vector<MPLPIndexType>& partition_j, MPLPIndexType var_j_size, MulDimArr* edge_belief, std::vector<std::vector<double> >& max_i_bij_not_xi) {
    double max, sec_max;
    max = sec_max = -MPLP_Inf;
    //double whole_i[var_i_size];
    double whole_j[var_j_size];
    /*for (MPLPIndexType i = 0; i < var_i_size; i++) {
        whole_i[i] = 0;
    }*/
    for (MPLPIndexType j = 0; j < var_j_size; j++) {
        whole_j[j] = 0;
    }
    //whole_i[single_i] = 1;
    for (MPLPIndexType j = 0; j < partition_j.size(); j++) {
        //brute force max_{pi(x_i)=pi(x_j)=1}bij(x_i,x_j)
        whole_j[partition_j[j]] = 1;
        std::vector<MPLPIndexType> inds; inds.push_back(single_i); inds.push_back(partition_j[j]);
        double tmp = edge_belief->GetVal(inds);
        if (max < tmp) max = tmp;
    }

    for (MPLPIndexType j = 0; j < var_j_size; j++) {
        if (whole_j[j] == 0) {
            //bruteforce max_{pi(x_i)=pi(x_j)=0}bij(x_i,x_j)
            std::vector<MPLPIndexType> inds; inds.push_back(single_i); inds.push_back(j);
            double tmp = edge_belief->GetVal(inds);
            if (sec_max < tmp) sec_max = tmp;

            if (max < max_i_bij_not_xi[single_i][j]) max = max_i_bij_not_xi[single_i][j];
        }
        else {
            if (sec_max < max_i_bij_not_xi[single_i][j]) sec_max = max_i_bij_not_xi[single_i][j];
        }
    }

    return max - sec_max;
}

// Compute a single edge weight in the projection graph (more efficiently)
double find_smn_state_j(const std::vector<MPLPIndexType>& partition_i, MPLPIndexType var_i_size, MPLPIndexType single_j, MPLPIndexType var_j_size, MulDimArr* edge_belief, std::vector<std::vector<double> >& max_j_bij_not_xj) {
    double max, sec_max;
    max = sec_max = -MPLP_Inf;
    double whole_i[var_i_size];
    //double whole_j[var_j_size];
    for (MPLPIndexType i = 0; i < var_i_size; i++) {
        whole_i[i] = 0;
    }
    /*for (MPLPIndexType j = 0; j < var_j_size; j++) {
        whole_j[j] = 0;
    }
    whole_j[single_j] = 1;*/
    for (MPLPIndexType i = 0; i < partition_i.size(); i++) {
        //brute force max_{pi(x_i)=pi(x_j)=1}bij(x_i,x_j)
        whole_i[partition_i[i]] = 1;
        std::vector<MPLPIndexType> inds; inds.push_back(partition_i[i]); inds.push_back(single_j);
        double tmp = edge_belief->GetVal(inds);
        if (max < tmp) max = tmp;
    }

    for (MPLPIndexType i = 0; i < var_i_size; i++) {
        if (whole_i[i] == 0) {
            //bruteforce max_{pi(x_i)=pi(x_j)=0}bij(x_i,x_j)
            std::vector<MPLPIndexType> inds; inds.push_back(i); inds.push_back(single_j);
            double tmp = edge_belief->GetVal(inds);
            if (sec_max < tmp) sec_max = tmp;

            if (max < max_j_bij_not_xj[i][single_j]) max = max_j_bij_not_xj[i][single_j];
        }
        else {
            if (sec_max < max_j_bij_not_xj[i][single_j]) sec_max = max_j_bij_not_xj[i][single_j];
        }
    }

    return max - sec_max;
}


// This will find the best partitioning of the states of a single edge
// Implements the FindPartition algorithm described in supplementary material of UAI 2012 paper.
void find_partition(std::vector<std::map<std::vector<MPLPIndexType>, MPLPIndexType> >& partition_set, MPLPAlg& mplp, MulDimArr* edge_belief, MPLPIndexType index_i, MPLPIndexType index_j) {

    MPLPIndexType size_i = mplp.m_var_sizes[index_i];
    MPLPIndexType size_j = mplp.m_var_sizes[index_j];

    std::vector<Node*> x_i; //partition for node i
    std::vector<Node*> x_j; //partition for node j

    for (MPLPIndexType i = 0; i < size_i; i++) {
        Node* t = new Node(2);
        x_i.push_back(t);
    }
    for (MPLPIndexType j = 0; j < size_j; j++) {
        Node* t = new Node(1);
        x_j.push_back(t);
    }

    // sort edges by their weights in decreasing order
    std::vector<std::pair<std::pair<MPLPIndexType, MPLPIndexType>, double> > sorted_edge_set;
    for (MPLPIndexType i = 0; i < size_i; i++) {
        for (MPLPIndexType j = 0; j < size_j; j++) {
            std::vector<MPLPIndexType> inds; inds.push_back(i); inds.push_back(j);
            double temp_val = edge_belief->GetVal(inds);
            sorted_edge_set.push_back(std::pair<std::pair<MPLPIndexType, MPLPIndexType>, double>(std::make_pair(std::make_pair(i, j), temp_val)));
        }
    }
    sort(sorted_edge_set.begin(), sorted_edge_set.end(), edge_sorting3);

    std::vector<std::pair<std::pair<MPLPIndexType, MPLPIndexType>, double> >::iterator it = sorted_edge_set.begin();
    double val_s = it->second; //val_s = bij(x_i^*, x_j^*);

    // union find step;
    for ( ; it != sorted_edge_set.end(); it++) {
        MPLPIndexType ind_i, ind_j;
        ind_i = it->first.first; //state of x_i
        ind_j = it->first.second; //staet of x_j
        Node* i_root = find(x_i[ind_i]); //find head of ind_i
        Node* j_root = find(x_j[ind_j]); //find head of ind_j
        if (i_root == j_root) continue; // if i, j belong to the same partition already, nothing to do

        MPLPIndexType bit_i = i_root->bit; //bit of ind_i node
        MPLPIndexType bit_j = j_root->bit; //bit of ind_j node
        if (bit_i == 2 && bit_j == 1) { //two singletons
            ;
        }
        if (bit_i == 2 && bit_j == 3) {
            if (size_i == 2) { //if number of partition is 2, then break
                val_s -= it->second; //val_s = bij(x_i^*, x_j^*) - max_{pi(x_i) != pi(x_j)} bij(x_i, x_j)
                break;
            }
            assert(size_i > 0);
            size_i--;
        }
        if (bit_i == 3 && bit_j == 1) {
            if (size_j == 2) {
                val_s -= it->second;
                break;
            }
            assert(size_j > 0);
            size_j--;
        }
        if (bit_i == 3 && bit_j == 3) {
            if (size_i == 2 || size_j == 2) {
                val_s -= it->second;
                break;
            }
            assert(size_i > 0); assert(size_j > 0);
            size_i--; size_j--;
        }
        Node* new_root = merge(i_root, j_root); //if # of partiton > 2, merge two two partitions
        new_root->bit = 3;
    }

    if (size_i != 2 && size_j != 2) {
        std::cout << "something is wrong" << std::endl;
    }

    if (size_i == 1 || size_j == 1) {
        std::cout << "should not happen" << std::endl;
        return;
    }
    std::vector<MPLPIndexType> part_i, part_j;

    //find the smallest partition and use it as an index
    // If there's a tie, use the partition with the smallest index. Keep it sorted.
    MPLPIndexType min = mplp.m_var_sizes[index_i] + 1;
    Node* head = NULL;
    for (MPLPIndexType i = 0; i < mplp.m_var_sizes[index_i]; i++) {
        Node* t = find(x_i[i]);
        if ((size_i == 2 || t->bit == 3) && min > t->i_size) {
            min = t->i_size;
            head = t;
        }
    }
    for (MPLPIndexType i = 0; i < mplp.m_var_sizes[index_i]; i++) {
        Node* t = find(x_i[i]);
        if (head == t)
            part_i.push_back(i);
    }

    min = mplp.m_var_sizes[index_j] + 1;
    for (MPLPIndexType j = 0; j < mplp.m_var_sizes[index_j]; j++) {
        Node* t = find(x_j[j]);
        if ((size_j == 2 || t->bit == 3) && min > t->j_size) {
            min = t->j_size;
            head = t;
        }
    }
    for (MPLPIndexType j = 0; j < mplp.m_var_sizes[index_j]; j++) {
        Node* t = find(x_j[j]);
        if (head == t)
            part_j.push_back(j);
    }

    // We have one map per variable.

    if (val_s != 0) {
        // Don't add a partition that already exists
        if (partition_set[index_i].find(part_i) == partition_set[index_i].end()) {

            // Add this partition
            // Keep track of the number of times that this partition is used
            partition_set[index_i][part_i] = 1;
        }
        else
            // Keep track of the number of times that this partition is used
            partition_set[index_i][part_i]++;

        if (partition_set[index_j].find(part_j) == partition_set[index_j].end()) {

            // Add this partition. Keep track of the number of times that this partition is used
            partition_set[index_j][part_j] = 1;
        }
        else
            partition_set[index_j][part_j]++;
    }

    for (MPLPIndexType i = 0; i < mplp.m_var_sizes[index_i]; i++) {
        delete [] x_i[i];
    }
    for (MPLPIndexType j = 0; j < mplp.m_var_sizes[index_j]; j++) {
        delete [] x_j[j];
    }
}


// Create the expanded projection graph by including all singleton partitions and also
// all partitions found by calling FindPartition on all edges.
MPLPIndexType create_expanded_projection_graph(MPLPAlg& mplp, std::vector<MPLPIndexType>& projection_imap_var, std::vector<std::vector<std::pair<MPLPIndexType, double> > >& projection_adjacency_list, std::map<std::pair<MPLPIndexType, MPLPIndexType>, double>&projection_edge_weights, double* &array_of_sij, MPLPIndexType& array_of_sij_size, std::vector<std::vector<MPLPIndexType> >& partition_imap) {
    // projection_imap_var maps from projection node to variable
    // partition_imap maps from projection node to vector of states

    MPLPIndexType num_of_vars = mplp.m_var_sizes.size();
    std::vector<std::map<std::vector<MPLPIndexType>,MPLPIndexType> > partition_set;
    std::set<double> set_of_sij;
    MPLPIndexType num_projection_nodes = 0;

    for (MPLPIndexType i = 0; i < num_of_vars; i++) {
        std::map<std::vector<MPLPIndexType>,MPLPIndexType> p;
        partition_set.push_back(p);

        // push all singleton partitions
        for (MPLPIndexType i_state = 0; i_state < mplp.m_var_sizes[i]; i_state++) {
            std::vector<MPLPIndexType> single;
            single.push_back(i_state);
            partition_set[i][single] = 1;
        }
    }

    bool partition = true;
    clock_t find_partition_start_time;
    if (partition) {
        find_partition_start_time = clock();
        //THIS IS THE NEW PARTITIONING ALGORITHM

        for (mapType::const_iterator it = mplp.m_intersect_map.begin(); it != mplp.m_intersect_map.end(); ++it) {
            MPLPIndexType i=it->first.first; MPLPIndexType j=it->first.second;
            MPLPIndexType ij_intersect_loc = it->second;

            MulDimArr* edge_belief = &mplp.m_sum_into_intersects[ij_intersect_loc];
            if(mplp.m_all_intersects[ij_intersect_loc][0] != i)  { // swap i and j
                MPLPIndexType tmp_i = i;
                i = j;
                j = tmp_i;
            }

            // Check to see if i and j have at least two states each -- otherwise, cannot be part of any frustrated edge
            if(mplp.m_var_sizes[i] <= 1 || mplp.m_var_sizes[j] <= 1)
                continue;

            find_partition(partition_set, mplp, edge_belief, i, j);
        }
        clock_t find_partition_end_time = clock();
        double find_partition_total_time = (double)(find_partition_end_time - find_partition_start_time)/CLOCKS_PER_SEC;
        if (MPLP_DEBUG_MODE) {
            printf(" -- find_partition. Took %lg seconds\n", find_partition_total_time);
        }
    }

    // Create one projection node for each remaining partition. Overload value of partition_set
    // to now refer to the node number in the projection graph.
    for(MPLPIndexType i=0; i < partition_set.size(); i++) {
        for(std::map<std::vector<MPLPIndexType>, MPLPIndexType>::iterator it = partition_set[i].begin(); it != partition_set[i].end(); it++) {
            it->second = num_projection_nodes++;
            projection_imap_var.push_back(i);
            partition_imap.push_back(it->first);

            std::vector<std::pair<MPLPIndexType, double> > temp;
            projection_adjacency_list.push_back(temp);
        }
    }

    clock_t find_smn_start_time = clock();

    // Create projection graph edges for each edge of original graph and each pair of partitions
    for (mapType::const_iterator it = mplp.m_intersect_map.begin(); it != mplp.m_intersect_map.end(); ++it) {
        MPLPIndexType i = it->first.first; MPLPIndexType j = it->first.second;
        MPLPIndexType ij_intersect_loc = it->second;

        // Do some pre-processing to make the case of single state partitions very fast

        MulDimArr* edge_belief = &mplp.m_sum_into_intersects[ij_intersect_loc];
        if (mplp.m_all_intersects[ij_intersect_loc][0] != i) {
            MPLPIndexType tmp_i = i;
            i = j;
            j = tmp_i;
        }

        if (mplp.m_var_sizes[i] <= 1 || mplp.m_var_sizes[j] <= 1) continue;

        std::vector<std::vector<double> > max_i_bij_not_xi;
        std::vector<std::vector<double> > max_j_bij_not_xj;

        for(MPLPIndexType state1=0; state1 < mplp.m_var_sizes[i]; state1++) {

            // Find max over state2
            double largest_val = -MPLP_Inf;
            MPLPIndexType largest_ind = 0;
            for(MPLPIndexType state2=0; state2 < mplp.m_var_sizes[j]; state2++) {

                std::vector<MPLPIndexType> inds; inds.push_back(state1); inds.push_back(state2);
                double tmp_val = edge_belief->GetVal(inds);

                if(tmp_val > largest_val) {
                    largest_val = tmp_val;
                    largest_ind = state2;
                }
            }

            // Find second largest val over state2
            double sec_largest_val = -MPLP_Inf;
            for(MPLPIndexType state2=0; state2 < mplp.m_var_sizes[j]; state2++) {

                std::vector<MPLPIndexType> inds; inds.push_back(state1); inds.push_back(state2);
                double tmp_val = edge_belief->GetVal(inds);

                if(tmp_val > sec_largest_val && state2 != largest_ind) {
                    sec_largest_val = tmp_val;
                }
            }

            // Assign values
            for(MPLPIndexType state2=0; state2 < mplp.m_var_sizes[j]; state2++) {
                std::vector<double> state_j; max_j_bij_not_xj.push_back(state_j);
                max_j_bij_not_xj[state1].push_back(largest_val);
            }
            max_j_bij_not_xj[state1][largest_ind] = sec_largest_val;
        }


        for(MPLPIndexType state1=0; state1 < mplp.m_var_sizes[j]; state1++) {

            // Find max over state2
            double largest_val = -MPLP_Inf;
            MPLPIndexType largest_ind = 0;
            for(MPLPIndexType state2=0; state2 < mplp.m_var_sizes[i]; state2++) {

                std::vector<MPLPIndexType> inds; inds.push_back(state2); inds.push_back(state1);
                double tmp_val = edge_belief->GetVal(inds);

                if(tmp_val > largest_val) {
                    largest_val = tmp_val;
                    largest_ind = state2;
                }
            }

            // Find second largest val over state2
            double sec_largest_val = -MPLP_Inf;
            for(MPLPIndexType state2=0; state2 < mplp.m_var_sizes[i]; state2++) {

                std::vector<MPLPIndexType> inds; inds.push_back(state2); inds.push_back(state1);
                double tmp_val = edge_belief->GetVal(inds);

                if(tmp_val > sec_largest_val && state2 != largest_ind) {
                    sec_largest_val = tmp_val;
                }
            }

            // Assign values
            for(MPLPIndexType state2=0; state2 < mplp.m_var_sizes[i]; state2++) {
                std::vector<double> state_i; max_i_bij_not_xi.push_back(state_i);
                max_i_bij_not_xi[state2].push_back(largest_val);
            }
            max_i_bij_not_xi[largest_ind][state1] = sec_largest_val;
        }


        // decompose: max_{x_j!=x_j}max_{x_i != x_i'}. Then, use the above computations.
        double max_ij_bij_not_xi_xj[mplp.m_var_sizes[i]][mplp.m_var_sizes[j]];
        for(MPLPIndexType state1=0; state1 < mplp.m_var_sizes[i]; state1++) {

            // Find max over state2
            double largest_val = -MPLP_Inf;
            MPLPIndexType largest_ind = 0;
            for(MPLPIndexType state2=0; state2 < mplp.m_var_sizes[j]; state2++) {
                double tmp_val = max_i_bij_not_xi[state1][state2];

                if(tmp_val > largest_val) {
                    largest_val = tmp_val;
                    largest_ind = state2;
                }
            }

            // Find second largest val over state2
            double sec_largest_val = -MPLP_Inf;
            for(MPLPIndexType state2=0; state2 < mplp.m_var_sizes[j]; state2++) {
                double tmp_val = max_i_bij_not_xi[state1][state2];

                if(tmp_val > sec_largest_val && state2 != largest_ind) {
                    sec_largest_val = tmp_val;
                }
            }

            // Assign values
            for(MPLPIndexType state2=0; state2 < mplp.m_var_sizes[j]; state2++)
                max_ij_bij_not_xi_xj[state1][state2] = largest_val;
            max_ij_bij_not_xi_xj[state1][largest_ind] = sec_largest_val;
        }

        // Now, for each partition of node i and each partition of node j, compute edge weights
        // If the edge weight is non-zero, insert edge into adjacency list

        for(std::map<std::vector<MPLPIndexType>, MPLPIndexType>::iterator it_i = partition_set[i].begin(); it_i != partition_set[i].end(); it_i++) {
            MPLPIndexType n = it_i->second;

            for(std::map<std::vector<MPLPIndexType>, MPLPIndexType>::iterator it_j = partition_set[j].begin(); it_j != partition_set[j].end(); it_j++) {
                MPLPIndexType m = it_j->second;

                double smn = 0;
                if (it_i->first.size() == 1 && it_j->first.size() == 1) {
                    MPLPIndexType xi = it_i->first[0]; MPLPIndexType xj = it_j->first[0];
                    std::vector<MPLPIndexType> inds; inds.push_back(xi); inds.push_back(xj);
                    double tmp_val = edge_belief->GetVal(inds);
                    smn = std::max(tmp_val, max_ij_bij_not_xi_xj[xi][xj]) - std::max(max_i_bij_not_xi[xi][xj], max_j_bij_not_xj[xi][xj]);
                }
                else if (it_i->first.size() == 1) {
                    MPLPIndexType xi = it_i->first[0];
                    smn = find_smn_state_i(xi, mplp.m_var_sizes[i], it_j->first, mplp.m_var_sizes[j], edge_belief, max_i_bij_not_xi);
                }
                else if (it_j->first.size() == 1) {
                    MPLPIndexType xj = it_j->first[0];
                    smn = find_smn_state_j(it_i->first, mplp.m_var_sizes[i], xj, mplp.m_var_sizes[j], edge_belief, max_j_bij_not_xj);
                }
                else {
                    // This computes smn
                    smn = find_smn(it_i->first, mplp.m_var_sizes[i], it_j->first, mplp.m_var_sizes[j], edge_belief);
                }

                if (smn != 0) {
                    set_of_sij.insert(fabs(smn));

                    projection_adjacency_list[m].push_back(std::make_pair(n, smn));
                    projection_adjacency_list[n].push_back(std::make_pair(m, smn));

                    projection_edge_weights.insert(std::pair<std::pair<MPLPIndexType, MPLPIndexType>, double>(std::pair<MPLPIndexType, MPLPIndexType>(n, m), smn));
                    projection_edge_weights.insert(std::pair<std::pair<MPLPIndexType, MPLPIndexType>, double>(std::pair<MPLPIndexType, MPLPIndexType>(m, n), smn));
                }

            }
        }
    }
    clock_t find_smn_end_time = clock();
    double find_smn_total_time = (double)(find_smn_end_time - find_smn_start_time)/CLOCKS_PER_SEC;
    if (MPLP_DEBUG_MODE) {
        printf(" -- find_smn. Took %lg seconds\n", find_smn_total_time);
    }

    array_of_sij = new double[set_of_sij.size()];
    array_of_sij_size = 0;

    for (std::set<double>::iterator set_iter = set_of_sij.begin(); set_iter != set_of_sij.end(); set_iter++) {  //NOTE: PASSED AS FUNCTION ARG, DO NOT DELETE HERE!!
        array_of_sij[array_of_sij_size++] = *set_iter;
    }

    return num_projection_nodes;
}


// Creates the k-projection graph (just a single partition per variable)
void create_k_projection_graph(MPLPAlg& mplp, std::vector<std::vector<MPLPIndexType> > &projection_map, MPLPIndexType& num_projection_nodes, std::vector<MPLPIndexType> &projection_imap_var, std::vector<std::vector<MPLPIndexType> > &partition_imap, std::map<std::pair<MPLPIndexType, MPLPIndexType>, double>& projection_edge_weights, adj_type &projection_adjacency_list, double* &array_of_sij, MPLPIndexType& array_of_sij_size) {

    // TODO: make sure for binary variables that there is only one node per variable (rather than 2).
    // TODO: projection_edge_weights can likely be removed from this function and elsewhere.
    // TODO: most of these ints can be changed to be unsigned and/or fewer bits. Look into memory allocation.

    // Initialize the projection graph
    MPLPIndexType projection_node_iter = 0;
    for(MPLPIndexType node=0; node < mplp.m_var_sizes.size(); node++) {
        std::vector<MPLPIndexType> i;

        for(MPLPIndexType state=0; state < mplp.m_var_sizes[node]; state++) {
            std::vector<std::pair<MPLPIndexType, double> > temp;
            projection_adjacency_list.push_back(temp);
            i.push_back(projection_node_iter++);
        }
        projection_map.push_back(i);
    }
    // Construct inverse map
    num_projection_nodes = projection_node_iter;

    projection_node_iter = 0;
    for(MPLPIndexType node=0; node < mplp.m_var_sizes.size(); node++) {
        for(MPLPIndexType state=0; state < mplp.m_var_sizes[node]; state++) {
            projection_imap_var.push_back(node);

            std::vector<MPLPIndexType> tmp_state_vector;
            tmp_state_vector.push_back(state);
            partition_imap.push_back(tmp_state_vector);

            projection_node_iter++;
        }
    }

    // Iterate over all of the edges (we do this by looking at the edge intersection sets)
    std::set<double> set_of_sij;
    for(mapType::const_iterator it = mplp.m_intersect_map.begin(); it != mplp.m_intersect_map.end(); ++it) {
        // Get the two nodes i & j and the edge intersection set. Put in right order.
        MPLPIndexType i=it->first.first; MPLPIndexType j=it->first.second;
        MPLPIndexType ij_intersect_loc = it->second;

        MulDimArr* edge_belief = &mplp.m_sum_into_intersects[ij_intersect_loc];
        if(mplp.m_all_intersects[ij_intersect_loc][0] != i)  { // swap i and j
            MPLPIndexType tmp_i = i;
            i = j;
            j = tmp_i;
        }

        // Check to see if i and j have at least two states each -- otherwise, cannot be part of any frustrated edge
        if(mplp.m_var_sizes[i] <= 1 || mplp.m_var_sizes[j] <= 1)
            continue;

        // Do some pre-processing for speed.
        double max_j_bij_not_xj[mplp.m_var_sizes[i]][mplp.m_var_sizes[j]];
        double max_i_bij_not_xi[mplp.m_var_sizes[i]][mplp.m_var_sizes[j]];

        for(MPLPIndexType state1=0; state1 < mplp.m_var_sizes[i]; state1++) {

            // Find max over state2
            double largest_val = -MPLP_Inf;
            MPLPIndexType largest_ind = 0;
            for(MPLPIndexType state2=0; state2 < mplp.m_var_sizes[j]; state2++) {

                std::vector<MPLPIndexType> inds; inds.push_back(state1); inds.push_back(state2);
                double tmp_val = edge_belief->GetVal(inds);

                if(tmp_val > largest_val) {
                    largest_val = tmp_val;
                    largest_ind = state2;
                }
            }

            // Find second largest val over state2
            double sec_largest_val = -MPLP_Inf;
            for(MPLPIndexType state2=0; state2 < mplp.m_var_sizes[j]; state2++) {

                std::vector<MPLPIndexType> inds; inds.push_back(state1); inds.push_back(state2);
                double tmp_val = edge_belief->GetVal(inds);

                if(tmp_val > sec_largest_val && state2 != largest_ind) {
                    sec_largest_val = tmp_val;
                }
            }

            // Assign values
            for(MPLPIndexType state2=0; state2 < mplp.m_var_sizes[j]; state2++)
                max_j_bij_not_xj[state1][state2] = largest_val;
            max_j_bij_not_xj[state1][largest_ind] = sec_largest_val;
        }


        for(MPLPIndexType state1=0; state1 < mplp.m_var_sizes[j]; state1++) {

            // Find max over state2
            double largest_val = -MPLP_Inf;
            MPLPIndexType largest_ind = 0;
            for(MPLPIndexType state2=0; state2 < mplp.m_var_sizes[i]; state2++) {

                std::vector<MPLPIndexType> inds; inds.push_back(state2); inds.push_back(state1);
                double tmp_val = edge_belief->GetVal(inds);

                if(tmp_val > largest_val) {
                    largest_val = tmp_val;
                    largest_ind = state2;
                }
            }

            // Find second largest val over state2
            double sec_largest_val = -MPLP_Inf;
            for(MPLPIndexType state2=0; state2 < mplp.m_var_sizes[i]; state2++) {

                std::vector<MPLPIndexType> inds; inds.push_back(state2); inds.push_back(state1);
                double tmp_val = edge_belief->GetVal(inds);

                if(tmp_val > sec_largest_val && state2 != largest_ind) {
                    sec_largest_val = tmp_val;
                }
            }

            // Assign values
            for(MPLPIndexType state2=0; state2 < mplp.m_var_sizes[i]; state2++)
                max_i_bij_not_xi[state2][state1] = largest_val;
            max_i_bij_not_xi[largest_ind][state1] = sec_largest_val;
        }

        // decompose: max_{x_j!=x_j}max_{x_i != x_i'}. Then, use the above computations.
        double max_ij_bij_not_xi_xj[mplp.m_var_sizes[i]][mplp.m_var_sizes[j]];
        for(MPLPIndexType state1=0; state1 < mplp.m_var_sizes[i]; state1++) {

            // Find max over state2
            double largest_val = -MPLP_Inf;
            MPLPIndexType largest_ind = 0;
            for(MPLPIndexType state2=0; state2 < mplp.m_var_sizes[j]; state2++) {
                double tmp_val = max_i_bij_not_xi[state1][state2];

                if(tmp_val > largest_val) {
                    largest_val = tmp_val;
                    largest_ind = state2;
                }
            }

            // Find second largest val over state2
            double sec_largest_val = -MPLP_Inf;
            for(MPLPIndexType state2=0; state2 < mplp.m_var_sizes[j]; state2++) {
                double tmp_val = max_i_bij_not_xi[state1][state2];

                if(tmp_val > sec_largest_val && state2 != largest_ind) {
                    sec_largest_val = tmp_val;
                }
            }

            // Assign values
            for(MPLPIndexType state2=0; state2 < mplp.m_var_sizes[j]; state2++)
                max_ij_bij_not_xi_xj[state1][state2] = largest_val;
            max_ij_bij_not_xi_xj[state1][largest_ind] = sec_largest_val;
        }

        // For each of their states
        for(MPLPIndexType xi=0; xi < mplp.m_var_sizes[i]; xi++) {
            MPLPIndexType m = projection_map[i][xi];

            for(MPLPIndexType xj=0; xj < mplp.m_var_sizes[j]; xj++) {
                MPLPIndexType n = projection_map[j][xj];

                std::vector<MPLPIndexType> inds; inds.push_back(xi); inds.push_back(xj);
                double tmp_val = edge_belief->GetVal(inds);

                // Compute s_mn for this edge
                double val_s = std::max(tmp_val, max_ij_bij_not_xi_xj[xi][xj]) - std::max(max_i_bij_not_xi[xi][xj], max_j_bij_not_xj[xi][xj]);

                // TODO: use threshold here, to make next stage faster
                if(val_s != 0) {
                    projection_adjacency_list[m].push_back(std::make_pair(n, val_s));
                    projection_adjacency_list[n].push_back(std::make_pair(m, val_s));
                    set_of_sij.insert(fabs(val_s));
                }

                // Insert into edge weight map
                projection_edge_weights.insert(std::pair<std::pair<MPLPIndexType,MPLPIndexType>,double>(std::pair<MPLPIndexType,MPLPIndexType>(n, m), val_s));
                projection_edge_weights.insert(std::pair<std::pair<MPLPIndexType,MPLPIndexType>,double>(std::pair<MPLPIndexType,MPLPIndexType>(m, n), val_s));
            }
        }

    }

    // Sort list_of_sij and remove duplicates
    array_of_sij = new double[set_of_sij.size()];
    array_of_sij_size = 0;

    for(std::set<double>::iterator set_iter = set_of_sij.begin(); set_iter != set_of_sij.end(); ++set_iter) {
        array_of_sij[array_of_sij_size++] = *set_iter;
    }
}


// Does binary search over the edge weights to find largest edge weight
// such that there is an odd-signed cycle.
double find_optimal_R(adj_type &projection_adjacency_list, double* &array_of_sij, MPLPIndexType array_of_sij_size) {

    // Do binary search over sij
    MPLPIndexType bin_search_lower_bound = 0;
    MPLPIndexType bin_search_upper_bound = array_of_sij_size;
    double sij_min = -1;
    MPLPIndexType num_projection_nodes = projection_adjacency_list.size();

    while(bin_search_lower_bound <= bin_search_upper_bound) {

        // Compute mid-point
        MPLPIndexType R_pos = floor((bin_search_lower_bound + bin_search_upper_bound)/2);
        double R = array_of_sij[R_pos];

        // Does there exist an odd signed cycle using just edges with sij >= R? If so, go up. If not, go down.
        bool found_odd_signed_cycle = false;

        // Initialize
        // TODO: do not allocate memory for this every time! Just re-initialize.

        int node_sign[num_projection_nodes];
        for(MPLPIndexType m=0; m<num_projection_nodes; m++) {
            node_sign[m] = 0; // denotes "not yet seen"
        }

        // Graph may be disconnected, so check from all nodes
        for(MPLPIndexType i = 0; i < num_projection_nodes && !found_odd_signed_cycle; i++) {
            if(node_sign[i] == 0) {
                node_sign[i] = 1; // root node
                std::queue<MPLPIndexType> q;
                q.push(i);

                while (!q.empty() && !found_odd_signed_cycle) {
                    MPLPIndexType current = q.front();
                    q.pop();

                    for (MPLPIndexType j = 0; j < projection_adjacency_list[current].size(); j++) {
                        double smn = projection_adjacency_list[current][j].second;
                        // Ignore edges with weight less than R
                        if (fabs(smn) < R)
                            continue;

                        MPLPIndexType next = projection_adjacency_list[current][j].first;
                        int sign_of_smn = (smn > 0) - (smn < 0);

                        if (node_sign[next] == 0) {
                            node_sign[next] = node_sign[current] * sign_of_smn;
                            q.push(next);
                        }
                        else if(node_sign[next] == -node_sign[current] * sign_of_smn) {
                            // Found an odd-signed cycle! Can quit.
                            found_odd_signed_cycle = true;
                            break;
                        }
                    }
                }
            }
        }

        if(found_odd_signed_cycle) {
            sij_min = R;
            bin_search_lower_bound = R_pos+1;
        }
        else
            bin_search_upper_bound = R_pos-1;
    }

    return sij_min;
}


// Returns an array which is a random permutation of the numbers 0 through n-1.
MPLPIndexType* random_permutation(MPLPIndexType n) {
    MPLPIndexType *p = new MPLPIndexType[n];
    for (MPLPIndexType i = 0; i < n; ++i) {
        MPLPIndexType j = rand() % (i + 1);
        p[i] = p[j];
        p[j] = i;
    }
    return p;
}


// Given an undirected graph, finds an odd-signed cycle.
// This works by breath-first search. Better might be to find minimal depth tree.
// NOTE: this function depends on the random seed becase it creates a random spanning tree.
void FindCycles(std::vector<std::list<MPLPIndexType> > &cycle_set, double optimal_R, MPLPIndexType ncycles_to_add, adj_type &projection_adjacency_list) {

    double R = optimal_R;
    MPLPIndexType num_projection_nodes = projection_adjacency_list.size();

    // Initialize  (NOTE: uses heap allocation)
    int *node_sign = new int[num_projection_nodes];
    MPLPIndexType *node_parent = new MPLPIndexType[num_projection_nodes], *node_depth = new MPLPIndexType[num_projection_nodes];

    for(MPLPIndexType i = 0; i < num_projection_nodes; i++) {
        node_sign[i] = 0; // denotes "not yet seen"
        node_depth[i] = 0;
        node_parent[i] = 0;
    }

    // construct the rooted spanning tree(s) -- randomly choose a root!
    MPLPIndexType* random_projection_node = random_permutation(num_projection_nodes);
    //  if(MPLP_DEBUG_MODE) {
    //    cout << "First random node is: " << random_projection_node[0] << endl;
    //  }

    for(MPLPIndexType ri = 0; ri < num_projection_nodes; ri++) {
        MPLPIndexType i = random_projection_node[ri];
        if(node_sign[i] == 0) {
            node_sign[i] = 1; // root node
            node_parent[i] = i;
            node_depth[i] = 0;

            std::queue<MPLPIndexType> q;
            q.push(i);

            while (!q.empty()) {
                MPLPIndexType current = q.front();
                q.pop();

                // Randomize the adjacency list
                MPLPIndexType* random_nbrs = random_permutation(projection_adjacency_list[current].size());
                for (MPLPIndexType rj = 0; rj < projection_adjacency_list[current].size(); rj++) {
                    MPLPIndexType j = random_nbrs[rj];
                    MPLPIndexType next = projection_adjacency_list[current][j].first;
                    if (node_sign[next] == 0) {
                        double smn = projection_adjacency_list[current][j].second;
                        if (fabs(smn) < R) {
                            ;
                        }
                        else {
                            int sign_of_smn = (smn > 0) - (smn < 0);
                            node_sign[next] = node_sign[current] * sign_of_smn;
                            node_parent[next] = current;
                            node_depth[next] = node_depth[current] + 1;
                            q.push(next);
                        }
                    }
                }
                delete [] random_nbrs;
            }
        }
    }
    delete [] random_projection_node;

    // construct edge set that contains edges that are not parts of the tree
    // TODO: This can be made faster! Can be performed in CONSTANT time.
    std::vector<std::list<MPLPIndexType> > edge_set;
    std::map<std::pair<MPLPIndexType, MPLPIndexType>, MPLPIndexType> edge_map;

    for (MPLPIndexType i = 0; i < num_projection_nodes; i++) {
        for (MPLPIndexType j = 0; j < projection_adjacency_list[i].size(); j++) {
            if (node_parent[i] == j || node_parent[j] == i) {
                continue;
            }
            double smn = projection_adjacency_list[i][j].second;

            if (fabs(smn) < R) {
                continue;
            }
            MPLPIndexType jj = projection_adjacency_list[i][j].first;
            int sign_of_smn = (smn > 0) - (smn < 0);

            if (node_sign[i] == -node_sign[jj] * sign_of_smn) { //cycle found
                MPLPIndexType depth_of_i, depth_of_j, temp_di, temp_dj;
                depth_of_i = temp_di = node_depth[i];
                depth_of_j = temp_dj = node_depth[jj];

                MPLPIndexType anc_i, anc_j;
                anc_i = i;
                anc_j = jj;

                while (temp_dj > temp_di) {
                    anc_j = node_parent[anc_j];
                    assert(temp_dj > 0);
                    temp_dj--;
                }
                while (temp_di > temp_dj) {
                    anc_i = node_parent[anc_i];
                    assert(temp_di > 0);
                    temp_di--;
                }
                while (temp_di >= 0) {
                    if (anc_i == anc_j) { //least common ancestor found
                        std::list<MPLPIndexType> temp;
                        temp.push_back(i);
                        temp.push_back(jj);
                        assert(depth_of_i >= temp_di); assert(depth_of_j >= temp_dj);
                        temp.push_back((depth_of_i - temp_di) + (depth_of_j - temp_dj));
                        std::map<std::pair<MPLPIndexType, MPLPIndexType>, MPLPIndexType>::iterator m_it = edge_map.find(std::make_pair(i, jj));
                        if (m_it == edge_map.end()) {
                            edge_map[std::pair<MPLPIndexType, MPLPIndexType>(i, jj)] = (depth_of_i - temp_di) + (depth_of_j - temp_dj);
                            edge_set.push_back(temp);
                        }
                        break;
                    }
                    anc_i = node_parent[anc_i];
                    anc_j = node_parent[anc_j];
                    assert(temp_di > 0);
                    temp_di--;
                    assert(temp_dj > 0);
                    temp_dj--;
                }
            }
        }
    }

    // sort the edge by the distance between two nodes in the tree
    sort(edge_set.begin(), edge_set.end(), edge_sorting);

    // Note: here we do not check for duplicate variables
    for (MPLPIndexType i = 0; i < edge_set.size() && cycle_set.size() < ncycles_to_add; i++) {

        // Get next edge (and depth) from sorted list
        std::list<MPLPIndexType> temp = edge_set.back();
        edge_set.pop_back();

        // Find least common ancestor
        MPLPIndexType left, right, left_d, right_d, ancestor, left_anc, right_anc;
        left = temp.front();
        temp.pop_front();
        right = temp.front();
        left_d = node_depth[left];
        right_d = node_depth[right];
        left_anc = left;
        right_anc = right;

        // find the first common ancestor
        while (left_d > right_d) {
            left_anc = node_parent[left_anc];
            assert(left_d > 0);
            left_d--;
        }
        while (right_d > left_d) {
            right_anc = node_parent[right_anc];
            assert(right_d > 0);
            right_d--;
        }
        while (left_anc != right_anc) {
            left_anc = node_parent[left_anc];
            right_anc = node_parent[right_anc];
        }
        ancestor = left_anc;

        std::list<MPLPIndexType> cycle;

        // backtrace the found cycle
        while (node_parent[left] != ancestor) {
            cycle.push_back(left);
            left = node_parent[left];
        }
        cycle.push_back(left);
        cycle.push_back(node_parent[left]);
        while (node_parent[right] != ancestor) {
            cycle.push_front(right);
            right = node_parent[right];
        }
        cycle.push_front(right);
        cycle_set.push_back(cycle);
    }
    delete [] node_sign;
    delete [] node_parent;
    delete [] node_depth;
}


// Check to see if there are any duplicate variables, and, if so, shortcut
// NOTE: this is not currently being used.
void shortcut(std::list<MPLPIndexType> &cycle, std::vector<MPLPIndexType> &projection_imap_var, std::map<std::pair<MPLPIndexType, MPLPIndexType>, double>& projection_edge_weights, MPLPIndexType& num_projection_nodes) {

    bool exist_duplicates = true;
    MPLPIndexType cycle_array[cycle.size()]; MPLPIndexType tmp_ind = 0;
    int cycle_sign[cycle.size()];

    for (std::list<MPLPIndexType>::iterator it=cycle.begin(); it!=cycle.end(); ++it) {
        cycle_array[tmp_ind++] = *it;
    }

    MPLPIndexType cycle_start = 0;
    assert(cycle.size() > 0);
    MPLPIndexType cycle_end = cycle.size()-1;

    // first, shorten the cycle
    while(exist_duplicates) {
        exist_duplicates = false;

        // This map allows us to quickly find duplicates and their locations within the cycle
        std::map<MPLPIndexType, MPLPIndexType> duplicates_map;

        // Must keep track of sign.
        cycle_sign[cycle_start] = 1;

        duplicates_map.insert(std::pair<MPLPIndexType,MPLPIndexType>(projection_imap_var[cycle_array[cycle_start]], cycle_start));
        for(MPLPIndexType i=cycle_start+1; i <= cycle_end; i++) {

            // Get edge weight
            std::map< std::pair<MPLPIndexType,MPLPIndexType>, double>::iterator weights_iter = projection_edge_weights.find(std::make_pair(cycle_array[i-1], cycle_array[i]));

            double smn = weights_iter->second;
            int sign_of_smn = (smn > 0) - (smn < 0);
            cycle_sign[i] = cycle_sign[i-1]*sign_of_smn;

            // Does node i already exist in the map?
            std::map<MPLPIndexType, MPLPIndexType>::iterator dups_iter = duplicates_map.find(projection_imap_var[cycle_array[i]]);
            if( dups_iter != duplicates_map.end() ) { // Duplicate found
                MPLPIndexType first_occurence_index = dups_iter->second;

                // Look to see if the first half of the cycle is violated.
                weights_iter = projection_edge_weights.find(std::make_pair(cycle_array[first_occurence_index], cycle_array[i-1]));
                double edge_weight = weights_iter->second;
                int sign_of_edge = (edge_weight > 0) - (edge_weight < 0);

                // NOTE: Bound decrease could be LESS than promised!
                if(cycle_sign[first_occurence_index] == -cycle_sign[i-1] * sign_of_edge) {

                    // Found a violated cycle! Can quit
                    cycle_start = first_occurence_index;
                    assert(i > 0);
                    cycle_end = i-1;
                }
                else {

                    // Otherwise, cut out the part of the cycle from first_occurence_index to cycle_end, and continue
                    assert(first_occurence_index > 0);
                    for(MPLPIndexType pos=first_occurence_index-1; pos>= cycle_start; pos--)
                        cycle_array[pos+i-first_occurence_index] = cycle_array[pos];
                    cycle_start += i-first_occurence_index;

                    exist_duplicates = true;
                }

                break;
            }
            else {
                // Insert into the map
                duplicates_map.insert(std::pair<MPLPIndexType,MPLPIndexType>(projection_imap_var[cycle_array[i]], i));
            }
        }
    }

    // Modify the cycle
    cycle.clear();
    for(MPLPIndexType i=cycle_start; i <= cycle_end; i++)
        cycle.push_back(cycle_array[i]);
}


MPLPIndexType add_cycle(MPLPAlg& mplp, std::list<MPLPIndexType> &cycle, std::vector<MPLPIndexType> &projection_imap_var, std::map<std::vector<MPLPIndexType>, bool >& triplet_set, MPLPIndexType& num_projection_nodes) {

    MPLPIndexType nClustersAdded = 0;

    // Number of clusters we're adding is length_cycle - 2
    assert(cycle.size() > 1);
    MPLPIndexType nNewClusters = cycle.size() - 2;
    TripletCluster newCluster[nNewClusters];
    MPLPIndexType cluster_index = 0;

    // Convert cycle to array
    MPLPIndexType cycle_array[cycle.size()]; MPLPIndexType tmp_ind = 0;
    for (std::list<MPLPIndexType>::iterator it=cycle.begin(); it!=cycle.end(); ++it) { // this is unnecesary
        //    if (MPLP_DEBUG_MODE) {
        //      if (*it >= num_projection_nodes) {
        //        cout << "Cycle uses non-trivial partitionining." << endl;
        //      }
        //    }
        cycle_array[tmp_ind++] = *it;
    }

    // Found violated cycle, now triangulate and add to the relaxation!
    //for(MPLPIndexType i=0; i+1 < cycle.size()-2-i; i++) {
    for(MPLPIndexType i=0; (2*i)+3 < cycle.size(); i++) {
        // Add projection_imap_var applied to [i, i+1, cycle.size()-2-i]
        newCluster[cluster_index].i = projection_imap_var[cycle_array[i]];
        newCluster[cluster_index].j = projection_imap_var[cycle_array[i+1]];
        newCluster[cluster_index].k = projection_imap_var[cycle_array[cycle.size()-2-i]];

        cluster_index++;
    }

    //for(MPLPIndexType i=cycle.size()-1; i-1 > cycle.size()-1-i; i--) {
    for(MPLPIndexType i=cycle.size()-1; 2*i > cycle.size(); i--) {
        // Add projection_imap_var applied to [i, i-1, cycle.size()-1-i]
        newCluster[cluster_index].i = projection_imap_var[cycle_array[i]];
        newCluster[cluster_index].j = projection_imap_var[cycle_array[i-1]];
        newCluster[cluster_index].k = projection_imap_var[cycle_array[cycle.size()-1-i]];

        cluster_index++;
    }

    // Add the top nclus_to_add clusters to the relaxation
    for(MPLPIndexType clusterId = 0; clusterId < nNewClusters; clusterId++) {
        // Check that these clusters and intersection sets haven't already been added
        std::vector<MPLPIndexType> temp;
        temp.push_back(newCluster[clusterId].i);
        temp.push_back(newCluster[clusterId].j);
        temp.push_back(newCluster[clusterId].k);
        sort(temp.begin(), temp.end());

        // Check to see if this cluster involves two of the same variables
        // (could happen because we didn't shortcut)
        if(temp[0] == temp[1] || temp[0] == temp[2] || temp[1] == temp[2]) {
            //      if(MPLP_DEBUG_MODE)
            //	cout << "skipping this triplet because it is an edge." << endl;
            continue;
        }

        std::map<std::vector<MPLPIndexType>, bool >::iterator t_itr = triplet_set.find(temp);
        if (t_itr == triplet_set.end())
            triplet_set.insert(std::pair<std::vector<MPLPIndexType>, bool >(temp, true));
        else {
            //	if(MPLP_DEBUG_MODE)
            //	  cout << "   Triplet was already present. Skipping..." << endl;
            continue;
        }

        // Find the intersection sets for this triangle
        std::vector<MPLPIndexType> ij_edge; ij_edge.push_back(newCluster[clusterId].i); ij_edge.push_back(newCluster[clusterId].j);
        const int temp_ij_intersect_loc = mplp.FindIntersectionSet(ij_edge);
        // This edge intersection set may not already exist, in which case we should add it
        if(temp_ij_intersect_loc == -1) {
            newCluster[clusterId].ij_intersect_loc = mplp.AddIntersectionSet(ij_edge);
        } else {
            newCluster[clusterId].ij_intersect_loc = static_cast<MPLPIndexType>(temp_ij_intersect_loc);
        }

        std::vector<MPLPIndexType> jk_edge; jk_edge.push_back(newCluster[clusterId].j); jk_edge.push_back(newCluster[clusterId].k);
        const int temp_jk_intersect_loc = mplp.FindIntersectionSet(jk_edge);
        // This edge intersection set may not already exist, in which case we should add it
        if(temp_jk_intersect_loc == -1) {
            newCluster[clusterId].jk_intersect_loc = mplp.AddIntersectionSet(jk_edge);
        } else {
            newCluster[clusterId].jk_intersect_loc = static_cast<MPLPIndexType>(temp_jk_intersect_loc);
        }

        std::vector<MPLPIndexType> ki_edge; ki_edge.push_back(newCluster[clusterId].k); ki_edge.push_back(newCluster[clusterId].i);
        const int temp_ki_intersect_loc = mplp.FindIntersectionSet(ki_edge);
        // This edge intersection set may not already exist, in which case we should add it
        if(temp_ki_intersect_loc == -1) {
            newCluster[clusterId].ki_intersect_loc = mplp.AddIntersectionSet(ki_edge);
        } else {
            newCluster[clusterId].ki_intersect_loc = temp_ki_intersect_loc;
        }

        // Now add cluster ijk
        std::vector<MPLPIndexType> ijk_inds;
        ijk_inds.push_back(newCluster[clusterId].i); ijk_inds.push_back(newCluster[clusterId].j); ijk_inds.push_back(newCluster[clusterId].k);

        std::vector<MPLPIndexType> ijk_intersect_inds;
        ijk_intersect_inds.push_back(newCluster[clusterId].ij_intersect_loc);
        ijk_intersect_inds.push_back(newCluster[clusterId].jk_intersect_loc);
        ijk_intersect_inds.push_back(newCluster[clusterId].ki_intersect_loc);

        mplp.AddRegion(ijk_inds, ijk_intersect_inds);

        // TODO: log which clusters are chosen...

        nClustersAdded++;
    }

    return nClustersAdded;
}


/**
 * Main function implementing UAI 2012 cycle tightening algorithm.
 * 
 * method=1: use create_k_projection_graph
 * method=2: use create_expanded_projection_graph
 */
MPLPIndexType TightenCycle(MPLPAlg & mplp, MPLPIndexType nclus_to_add,  std::map<std::vector<MPLPIndexType>, bool >& triplet_set, double & promised_bound, MPLPIndexType method) {

    MPLPIndexType nClustersAdded = 0;
    //MPLPIndexType nNewClusters;

    if (MPLP_DEBUG_MODE) std::cout << "Finding the most violated cycle...." << std::endl;

    // This map allows us to quickly look up the edge weights
    std::map<std::pair<MPLPIndexType, MPLPIndexType>, double> projection_edge_weights;
    MPLPIndexType num_projection_nodes;
    std::vector<std::vector<MPLPIndexType> > projection_map;
    std::vector<MPLPIndexType> projection_imap_var;
    std::vector<std::vector<std::pair<MPLPIndexType, double> > > projection_adjacency_list;
    std::vector<std::vector<MPLPIndexType> > partition_imap;

    double* array_of_sij; MPLPIndexType array_of_sij_size;

    // Define the projection graph and all edge weights
    if(method == 2)
        num_projection_nodes = create_expanded_projection_graph(mplp, projection_imap_var, projection_adjacency_list, projection_edge_weights, array_of_sij, array_of_sij_size, partition_imap);
    else if(method == 1)
        create_k_projection_graph(mplp, projection_map, num_projection_nodes, projection_imap_var, partition_imap, projection_edge_weights, projection_adjacency_list, array_of_sij, array_of_sij_size);
    else {
        std::cout << "ERROR: method not defined." << std::endl;
        return 0;
    }

    std::vector<std::list<MPLPIndexType> > cycle_set;
    double optimal_R = find_optimal_R(projection_adjacency_list, array_of_sij, array_of_sij_size);
    if (MPLP_DEBUG_MODE) std::cout << "R_optimal = " << optimal_R << std::endl;

    promised_bound = optimal_R;

    // Look for cycles. Some will be discarded.
    clock_t start_time = clock();

    if (optimal_R > 0) {
        // TODO: this is almost certainly doing more computation than necessary. Might want to change
        // nclus_to_add*10 to nclus_to_add, and comment out all but the top 3.
        FindCycles(cycle_set, optimal_R, nclus_to_add*10, projection_adjacency_list);
        FindCycles(cycle_set, optimal_R/2, nclus_to_add*10, projection_adjacency_list);
        FindCycles(cycle_set, optimal_R/4, nclus_to_add*10, projection_adjacency_list);
        FindCycles(cycle_set, optimal_R/8, nclus_to_add*10, projection_adjacency_list);
        FindCycles(cycle_set, optimal_R/16, nclus_to_add*10, projection_adjacency_list);
        FindCycles(cycle_set, optimal_R/32, nclus_to_add*10, projection_adjacency_list);
        FindCycles(cycle_set, optimal_R/64, nclus_to_add*10, projection_adjacency_list);
        FindCycles(cycle_set, optimal_R/128, nclus_to_add*10, projection_adjacency_list);
    }

    clock_t end_time = clock();
    double total_time = (double)(end_time - start_time)/CLOCKS_PER_SEC;
    if (MPLP_DEBUG_MODE) {
        printf(" -- FindCycles. Took %lg seconds\n", total_time);
    } 


    // Add all cycles that we've found to the relaxation!
    start_time = clock();
    for (MPLPIndexType z = 0; z < cycle_set.size() && nClustersAdded < nclus_to_add; z++) {

        // Check to see if there are any duplicate nodes, and, if so, shortcut
        //    shortcut(cycle_set[z], projection_imap_var, projection_edge_weights, num_projection_nodes);

        // Output the cycle
        if (MPLP_DEBUG_MODE){
            for (std::list<MPLPIndexType>::iterator it=cycle_set[z].begin(); it!=cycle_set[z].end(); ++it) {
                std::cout << projection_imap_var[*it] << "(";
                std::vector<MPLPIndexType> temp = partition_imap[*it];
                assert(temp.size() > 0);
                for(MPLPIndexType s=0; s < temp.size()-1; s++)
                    std::cout << temp[s] << ",";
                std::cout << temp[temp.size()-1] << "), ";
            }
            std::cout << std::endl;
        }

        // Add cycle to the relaxation
        nClustersAdded += add_cycle(mplp, cycle_set[z], projection_imap_var, triplet_set, num_projection_nodes);
    }

    end_time = clock();
    total_time = (double)(end_time - start_time)/CLOCKS_PER_SEC;
    if (MPLP_DEBUG_MODE) {
        printf(" -- shortcut + add_cycles. Took %lg seconds\n", total_time);
    }

    delete_projection_graph(mplp.m_var_sizes.size(), projection_map, projection_imap_var, projection_adjacency_list, array_of_sij);
    return nClustersAdded;
}

} // namespace mplpLib

#endif
