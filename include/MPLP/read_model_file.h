#ifndef MPLP_READ_MODEL_FILE_H
#define MPLP_READ_MODEL_FILE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <limits.h>
#include <stdlib.h>
#include <math.h>

namespace mplpLib {

#define eps 1e-40
#define huge 1e40

int read_model_file(std::vector<int> &, std::vector< std::vector<int> > &, std::vector< std::vector<double> > &, const std::string, const std::string = "");  //at most 1 evidence
int read_model_file(std::vector<int> &, std::map<int, int> &, std::vector< std::vector<int> > &, std::vector< std::vector<double> > &, const std::string, const std::string);
int read_model_file2(std::vector<int> &, std::map<int, int> &, std::vector< std::vector<int> > &, std::vector< std::vector<double> > &, const std::string, const std::string);

int read_stereo_file(std::vector<int> &, std::vector< std::vector<int> > &, std::vector< std::vector<double> > &, const std::string);

void get_lambdas(int, int, int, int, const std::vector<int> &, const std::map<int, int> &, const std::vector<int> &, const std::vector<int> &, const std::vector<double> &, std::vector<double> &);

} // namespace mplpLib

#endif
