/*
 *  muldim_arr.cpp
 *  mplp
 *
 *  Created by Amir Globerson on 6/24/08.
 *  Copyright 2008 MIT. All rights reserved.
 *
 */

#include <math.h>
#include <float.h>
#include <cassert>

#include <MPLP/muldim_arr.h>

using namespace std;

void print_int_vec(vector<mplpLib::MPLPIndexType> v){
    for (mplpLib::MPLPIndexType i=0; i<v.size(); ++i){
        cout << v[i] << " ";
    }
}


mplpLib::MulDimArr::MulDimArr(vector<MPLPIndexType> & base_sizes)
{
    m_base_sizes = base_sizes;
    m_n_prodsize = 1;
    for (MPLPIndexType i=0; i<base_sizes.size(); i++)
        m_n_prodsize*= base_sizes[i];
    m_dat= new double[m_n_prodsize];

    // Initialize array to zero
    for(MPLPIndexType i=0; i< m_n_prodsize; i++)
        m_dat[i] = 0.0;

    m_ep = &m_dat[m_n_prodsize];	
}

mplpLib::MulDimArr::MulDimArr(const MulDimArr & v)
{
    m_base_sizes = v.m_base_sizes;
    m_n_prodsize = v.m_n_prodsize;
    m_dat= new double[m_n_prodsize];
    m_ep = &m_dat[m_n_prodsize];	
    memcpy(m_dat,v.m_dat,m_n_prodsize*sizeof(double));
}

void mplpLib::MulDimArr::print(void) const
{
    for (MPLPIndexType i=0; i<m_n_prodsize; i++)
    {
        cout << m_dat[i] << " ";
    }
    cout << endl;
}

// Print elements along with their indices
void mplpLib::MulDimArr::print_with_inds(void) const
{
    vector<MPLPIndexType> inds_for_big;

    for (MPLPIndexType i=0; i< m_base_sizes.size(); i++)
        inds_for_big.push_back(0);

    for (MPLPIndexType i=0; i<m_n_prodsize; i++)
    {
        cout << "[";
        print_int_vec(inds_for_big);
        cout << "]: ";
        cout << m_dat[i] << endl;
        BaseInc(inds_for_big);
    }
}

mplpLib::MulDimArr & mplpLib::MulDimArr::operator=(const MulDimArr & v)
{
    m_base_sizes = v.m_base_sizes;
    m_n_prodsize = v.m_n_prodsize;
    if (m_dat!=NULL)
        delete [] m_dat;
    m_dat= new double[m_n_prodsize];
    m_ep = &m_dat[m_n_prodsize];
    memcpy(m_dat,v.m_dat,m_n_prodsize*sizeof(double));
    return (*this);
}

// Return the index in the flat multi-dimensional array corresponding to
// the given multi-index
inline mplpLib::MPLPIndexType mplpLib::MulDimArr::GetFlatInd(vector<MPLPIndexType> & base_inds) const{
    MPLPIndexType y=0;
    MPLPIndexType fact = 1;
    MPLPIndexType nx = m_base_sizes.size();

    assert(nx > 0);
    for (MPLPIndexType i=nx;i>0; i--) {
        y += base_inds[i-1]*fact;
        fact*= m_base_sizes[i-1];
    }
    return y;
}

//this is for decoding purposes
void mplpLib::MulDimArr::GetInds(MPLPIndexType ind, vector<MPLPIndexType> & inds) const{
    MPLPIndexType c, N = m_base_sizes.size();
    inds.reserve(N);
    while (N--) {
        inds[N] = c = ind % m_base_sizes[N];
        //              inds.push_back(c = ind % m_base_sizes[N]);
        ind = (ind - c) / m_base_sizes[N];
    }
}


inline mplpLib::MPLPIndexType mplpLib::MulDimArr::GetFlatIndFromBig(vector<MPLPIndexType> big_base_inds, vector<MPLPIndexType> inds_in_big) const{
    MPLPIndexType y=0;
    MPLPIndexType fact = 1;
    MPLPIndexType nx = m_base_sizes.size();

    assert(nx > 0);
    for (MPLPIndexType i=nx;i>0; i--) {
        y += big_base_inds[inds_in_big[i-1]]*fact;
        fact*= m_base_sizes[i-1];
    }
    return y;
}



double mplpLib::MulDimArr::GetVal(vector<MPLPIndexType> & indices) const
{
    return m_dat[GetFlatInd(indices)];
}


inline void mplpLib::MulDimArr::BaseInc(vector<MPLPIndexType> & inds) const{
    MPLPIndexType nx = m_base_sizes.size();
    assert(nx > 0);
    MPLPIndexType pos = nx-1;
    inds[pos]++;
    while (inds[pos]==m_base_sizes[pos]) {
        inds[pos] = 0;
        if (pos==0)
            break;
        pos--;
        inds[pos]++;
    }
}


/* 
 *	Expand the current vector into a new (larger) multi dimensional vector.
 *   inds_in_big gives the indices of the small array variables in the big array
 *
 */
void mplpLib::MulDimArr::ExpandAndAdd(MulDimArr & big_to_add_to, vector<MPLPIndexType> & inds_of_small_in_big)
{
    // TODO: Verify that var_sizes_big agrees with the var_sizes of the smaller (this) array
    //	MulDimArr big_arr(var_sizes_big);
    // This will index the big array
    vector<MPLPIndexType> inds_for_big;

    // Initialize vector of indices for big array to zeros.
    for (MPLPIndexType i=0; i<big_to_add_to.m_base_sizes.size(); i++)
        inds_for_big.push_back(0);

    // We now go over the big indices one by one, and for each, take the value from the indices
    // corresponding to the small (this) array
    double *big_dat = &big_to_add_to.m_dat[0];
    for (MPLPIndexType vi=0; vi<big_to_add_to.m_n_prodsize; vi++)
    {
        // Get the flat index in the small array
        MPLPIndexType ind = GetFlatIndFromBigSpecial(inds_for_big,inds_of_small_in_big);
        // Put value in big array (note that the "safe" way to do this would be to translate inds_for_big into a flat index, but
        // we don't do this to save computation
        //		big_to_add_to.m_dat[vi] += m_dat[ind];
        (*big_dat)+=m_dat[ind];
        big_dat++;
        // Move to the next indices for big (note again that we assume this corresponds to a flat index of vi+1)
        big_to_add_to.BaseIncSpecial(inds_for_big);
    }
}	


/* 
 *	Expand the current vector into a new (larger) multi dimensional vector.
 *   inds_in_big gives the indices of the small array variables in the big array
 *
 */
void mplpLib::MulDimArr::ExpandAndSubtract(MulDimArr & big_to_sub_from, vector<MPLPIndexType> & inds_of_small_in_big)
{
    // TODO: Verify that var_sizes_big agrees with the var_sizes of the smaller (this) array
    //	MulDimArr big_arr(var_sizes_big);
    // This will index the big array
    vector<MPLPIndexType> inds_for_big;

    // Initialize vector of indices for big array to zeros.
    for (MPLPIndexType i=0; i<big_to_sub_from.m_base_sizes.size(); i++)
        inds_for_big.push_back(0);

    // We now go over the big indices one by one, and for each, take the value from the indices
    // corresponding to the small (this) array
    double *big_dat = &big_to_sub_from.m_dat[0];
    for (MPLPIndexType vi=0; vi<big_to_sub_from.m_n_prodsize; vi++)
    {
        // Get the flat index in the small array
        MPLPIndexType ind = GetFlatIndFromBigSpecial(inds_for_big,inds_of_small_in_big);
        // Put value in big array (note that the "safe" way to do this would be to translate inds_for_big into a flat index, but
        // we don't do this to save computation
        //		big_to_sub_from.m_dat[vi] += m_dat[ind];
        (*big_dat)-=m_dat[ind];
        big_dat++;
        // Move to the next indices for big (note again that we assume this corresponds to a flat index of vi+1)
        big_to_sub_from.BaseIncSpecial(inds_for_big);
    }
}	


/* 
 *	Expand the current vector into a new (larger) multi dimensional vector.
 *   inds_in_big gives the indices of the small array variables in the big array
 *
 */
mplpLib::MulDimArr mplpLib::MulDimArr::Expand(vector<MPLPIndexType> & var_sizes_big, vector<MPLPIndexType> & inds_of_small_in_big)
{
    // TODO: Verify that var_sizes_big agrees with the var_sizes of the smaller (this) array
    MulDimArr big_arr(var_sizes_big);
    // This will index the big array
    vector<MPLPIndexType> inds_for_big;

    // Initialize vector of indices for big array to zeros.
    for (MPLPIndexType i=0; i<var_sizes_big.size(); i++)
        inds_for_big.push_back(0);

    // We now go over the big indices one by one, and for each, take the value from the indices
    // corresponding to the small (this) array
    for (MPLPIndexType vi=0; vi<big_arr.m_n_prodsize; vi++)
    {
        // Get the flat index in the small array
        MPLPIndexType ind = GetFlatIndFromBigSpecial(inds_for_big,inds_of_small_in_big);
        // Put value in big array (note that the "safe" way to do this would be to translate inds_for_big into a flat index, but
        // we don't do this to save computation
        big_arr.m_dat[vi] = m_dat[ind];
        // Move to the next indices for big (note again that we assume this corresponds to a flat index of vi+1)
        big_arr.BaseIncSpecial(inds_for_big);
    }

    return big_arr;
}	


mplpLib::MulDimArr & mplpLib::MulDimArr::operator*=(double val)
        {
    double *p;
    for (p=m_dat; p<m_ep; p++)
        (*p)*=val;
    return (*this);
        }

mplpLib::MulDimArr & mplpLib::MulDimArr::operator=(double val)
{
    double *p;
    for (p=m_dat; p<m_ep; p++)
        (*p)=val;
    return (*this);
}

mplpLib::MulDimArr & mplpLib::MulDimArr::operator+=(MulDimArr & v)
        {
    double *p,*pv;
    for (p=m_dat, pv=v.m_dat; p<m_ep; p++, pv++)
        (*p)+=(*pv);
    return (*this);
        }

mplpLib::MulDimArr & mplpLib::MulDimArr::operator-=(MulDimArr & v)
        {
    double *p,*pv;
    for (p=m_dat, pv=v.m_dat; p<m_ep; p++, pv++)
        (*p)-=(*pv);
    return (*this);
        }

double mplpLib::MulDimArr::Max(MPLPIndexType &max_at) const
{
    double m = m_dat[0];
    max_at = 0;
    for (MPLPIndexType i=1;i<m_n_prodsize; i++) {
        if (m_dat[i]>=m)
        {
            max_at = i;
            m=m_dat[i];
        }
    }
    return m;
}


double mplpLib::MulDimArr::Entropy() const
{
    double sum = 0, e = 0, y;
    for (MPLPIndexType i = 0; i < m_n_prodsize; ++i){
        sum += y = exp(m_dat[i]);
        e -= y * m_dat[i];
    }
    return e / sum + log(sum);
}


double mplpLib::MulDimArr::Entropy_over_free_variables(const std::vector<MPLPIndexType> & fixed_variables, const std::vector<MPLPIndexType> & assignment) const{
    double s = 0, e = 0;
    _Entropy_over_free_variables(m_base_sizes.size(), fixed_variables.size(), 0, 1, fixed_variables, assignment, s, e);
    return e / s + log(s);
}

void mplpLib::MulDimArr::_Entropy_over_free_variables(MPLPIndexType l, MPLPIndexType p, MPLPIndexType base, MPLPIndexType fact, const std::vector<MPLPIndexType> & fixed_variables, const std::vector<MPLPIndexType> & assignment, double & sum, double & ent) const{
    double y;
    while (l--){
        if (p == 0 || l != fixed_variables[p - 1]){  //not a fixed variable
            for (MPLPIndexType i = 0; i < m_base_sizes[l]; ++i){
                _Entropy_over_free_variables(l, p, base + i * fact, fact * m_base_sizes[l], fixed_variables, assignment, sum, ent);
            }
            return;
        }
        --p;  //variable is fixed
        base += assignment[l] * fact;
        fact *= m_base_sizes[l];
    }
    sum += y = exp(m_dat[base]);
    ent -= y * m_dat[base];
}

double mplpLib::MulDimArr::max_over_free_variables(const vector<MPLPIndexType> & fixed_variables, vector<MPLPIndexType> & max_assignment) const{
    return _max_over_free_variables(m_base_sizes.size(), fixed_variables.size(), 0, 1, fixed_variables, max_assignment);
}

double mplpLib::MulDimArr::_max_over_free_variables(MPLPIndexType l, MPLPIndexType p, MPLPIndexType base, MPLPIndexType fact, const vector<MPLPIndexType> & fixed_variables, vector<MPLPIndexType> & max_assignment) const{
    double max = -DBL_MAX, _max;     //sort(fixed_variables.begin(), fixed_variables.end());
    while (l--){
        if (p == 0 || l != fixed_variables[p - 1]){  //not a fixed variable
            for (MPLPIndexType i = 0; i < m_base_sizes[l]; ++i){
                vector<MPLPIndexType> _max_assignment(max_assignment);
                _max_assignment[l] = i;
                if ( (_max = _max_over_free_variables(l, p, base + i * fact, fact * m_base_sizes[l], fixed_variables, _max_assignment)) > max ){
                    max = _max;
                    max_assignment = _max_assignment;
                }
            }
            return max;
        }
        --p;  //variable is fixed
        base += max_assignment[l] * fact;
        fact *= m_base_sizes[l];
    }
    return m_dat[base];   //all variables are fixed
}

double mplpLib::MulDimArr::gap_over_free_variables(const vector<MPLPIndexType> & fixed_variables, vector<MPLPIndexType> & max_assignment, double & m1, double & m2) const{
    MPLPIndexType l = m_base_sizes.size(), p = fixed_variables.size(), base = 0, fact = 1;
    double m;
    m1 = -MPLP_huge, m2 = -MPLP_huge;
    while (l--){
        if (p == 0 || l != fixed_variables[p - 1]){  //not a fixed variable
            for (MPLPIndexType i = 0; i < m_base_sizes[l]; ++i){
                vector<MPLPIndexType> _max_assignment(max_assignment);
                _max_assignment[l] = i;
                if ( (m = _max_over_free_variables(l, p, base + i * fact, fact * m_base_sizes[l], fixed_variables, _max_assignment)) > m1 ){
                    m2 = m1;
                    m1 = m;
                    max_assignment = _max_assignment;
                }else if (m > m2){
                    m2 = m;
                }
            }
            return m1 - m2;
        }
        --p;  //variable is fixed
        base += max_assignment[l] * fact;
        fact *= m_base_sizes[l];
    }
    return HUGE;   //all variables are fixed
}

// For a given subset of variables, for any assignment to the subset maximize over
// the value of all variables outside the subset
// TODO: Do this for multiple subsets simultaneously. Should be much more efficient!!
void mplpLib::MulDimArr::max_into_multiple_subsets_special(vector<vector<MPLPIndexType> > & all_subset_inds, vector<MulDimArr> & all_max_res) const
{
    // Prepare a MulDimArr into which we'll maximize
    // First set its size
    MPLPIndexType nSubsets = all_subset_inds.size();

    bool *b_need_max = new bool[nSubsets];

    for (MPLPIndexType si=0; si<nSubsets; si++)
    {
        // If the subset equals the big array	then maximizing would give us the subset (assuming there is no reordering, which I think we can)
        if (all_subset_inds[si].size()==m_base_sizes.size())
        {
            all_max_res[si] = (*this);
            b_need_max[si] = 0;
        }
        else
        {
            all_max_res[si] = -1e9;
            b_need_max[si] = 1;
        }
    }
    // Go over all values of the big array (this). For each check if its
    // value on the subset is larger than the current max
    // NOTE: We make the same assumption here as in the Expand function, regaring
    // the correspondence between the flat index (vi) and the the index vector inds_for_big
    vector<MPLPIndexType> inds_for_big;
    for (MPLPIndexType i=0; i< m_base_sizes.size(); i++)
        inds_for_big.push_back(0);

    for (MPLPIndexType vi=0; vi<m_n_prodsize; vi++)
    {
        for (MPLPIndexType si=0; si<nSubsets; si++)
        {
            if (!b_need_max[si])
                continue;
            MPLPIndexType flat_subind = all_max_res[si].GetFlatIndFromBigSpecial(inds_for_big,all_subset_inds[si]);
            all_max_res[si][flat_subind] = max(all_max_res[si][flat_subind],m_dat[vi]);
        }
        BaseIncSpecial(inds_for_big);
    }

    delete [] b_need_max;
}

void mplpLib::MulDimArr::Write(ofstream & ofs)
{
    for (MPLPIndexType i=0; i<m_n_prodsize; i++)
        ofs << m_dat[i] << " ";
    ofs << endl;
}

void mplpLib::MulDimArr::Read(ifstream & ifs)
{
    for (MPLPIndexType i=0; i<m_n_prodsize; i++)
        ifs >> m_dat[i];
}
