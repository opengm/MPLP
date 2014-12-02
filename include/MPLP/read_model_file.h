#ifndef MPLP_READ_MODEL_FILE_H
#define MPLP_READ_MODEL_FILE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <limits.h>
#include <stdlib.h>
#include <math.h>

#include <MPLP/mplp_config.h>

namespace mplpLib {

#define MPLP_eps 1e-40
#define MPLP_huge 1e40

MPLPIndexType read_model_file(std::vector<MPLPIndexType> &, std::vector< std::vector<MPLPIndexType> > &, std::vector< std::vector<double> > &, const std::string, const std::string = "");  //at most 1 evidence
MPLPIndexType read_model_file(std::vector<MPLPIndexType> &, std::map<MPLPIndexType, MPLPIndexType> &, std::vector< std::vector<MPLPIndexType> > &, std::vector< std::vector<double> > &, const std::string, const std::string);
MPLPIndexType read_model_file2(std::vector<MPLPIndexType> &, std::map<MPLPIndexType, MPLPIndexType> &, std::vector< std::vector<MPLPIndexType> > &, std::vector< std::vector<double> > &, const std::string, const std::string);

MPLPIndexType read_stereo_file(std::vector<MPLPIndexType> &, std::vector< std::vector<MPLPIndexType> > &, std::vector< std::vector<double> > &, const std::string);

void get_lambdas(MPLPIndexType, MPLPIndexType, MPLPIndexType, MPLPIndexType, const std::vector<MPLPIndexType> &, const std::map<MPLPIndexType, MPLPIndexType> &, const std::vector<MPLPIndexType> &, const std::vector<MPLPIndexType> &, const std::vector<double> &, std::vector<double> &);

} // namespace mplpLib

#endif
