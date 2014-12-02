#include <MPLP/read_model_file.h>

#define MPLP_DEBUG_MODE 0

mplpLib::MPLPIndexType mplpLib::read_model_file(std::vector<MPLPIndexType> & var_sizes, std::vector< std::vector<MPLPIndexType> > & all_factors, std::vector< std::vector<double> > & all_lambdas, const std::string fn, const std::string evid_fn){
    MPLPIndexType i, j, v, nvars, nfactors, nevid, evid_size, curr_var, curr_val;
    double val;
    //	vector<int> factor_size, dim;
    //vector **f;
    //double **l;
    std::string s, t;
    std::ifstream fstr;
    fstr.open(fn.c_str(), std::fstream::in);
    if (fstr.fail()){
        std::cerr<<"Error opening model file"<<std::endl;
        return 0;
    }
    fstr>>s;
    if (s != "MARKOV"){
        std::cerr<<"Error: file is not in MARKOV format."<<std::endl;
    }
    fstr>>nvars;
    if (MPLP_DEBUG_MODE){ std::cout<<"nvars = "<<nvars<<std::endl; }
    for (i = 0; i < nvars; ++i){
        fstr>>v;
        var_sizes.push_back(v);
        if (MPLP_DEBUG_MODE){ std::cout<<"var_sizes["<<i<<"] = "<<var_sizes[i]<<std::endl; }
    }
    fstr>>nfactors;
    if (MPLP_DEBUG_MODE){ std::cout<<"nfactors = "<<nfactors<<std::endl; }
    all_factors.reserve(nfactors);
    all_lambdas.reserve(nfactors);
    for (i = 0; i < nfactors; ++i){
        fstr>>nvars;
        all_factors.push_back(std::vector<MPLPIndexType>());
        all_lambdas.push_back(std::vector<double>());
        if (MPLP_DEBUG_MODE){ std::cout<<"factor_size["<<i<<"] = "<<nvars<<std::endl; }
        for (j = 0; j < nvars; ++j){
            fstr>>v;
            all_factors[i].push_back(v);
            if (MPLP_DEBUG_MODE){ std::cout<<"all_factors["<<i<<"]["<<j<<"] = "<<v<<std::endl; }
        }
    }//"all_lambdas": the lambdas (messages) are just log(function table entry) 's
    for (i = 0; i < nfactors; ++i){
        fstr>>nvars;  //the PIC file does not contain empty lines
        for (j = 0; j < nvars; ++j){
            fstr>>val;
            // We work in log space, so take the log of the factors' potentials
            all_lambdas[i].push_back(val > 0 ? log(val) : -MPLP_huge);
            if (MPLP_DEBUG_MODE){ std::cout<<"e^all_lambdas["<<i<<"]["<<j<<"] = "<<val<<std::endl; }
            if (MPLP_DEBUG_MODE){ std::cout<<"all_lambdas["<<i<<"]["<<j<<"] = "<<all_lambdas[i][j]<<std::endl; }
        }
    }
    fstr.close();
    if (evid_fn != ""){
        fstr.open(evid_fn.c_str(), std::fstream::in);
        if (fstr.fail()){
            std::cerr<<"Error opening evidence file"<<std::endl;
            return 0;
        }
        fstr>>nevid;   //for PIC we have at most 1 evidence
        if (nevid){    //add a local potential for the evidence
            fstr>>evid_size;
            while (evid_size--){
                fstr>>curr_var;
                fstr>>curr_val;
                if (MPLP_DEBUG_MODE){ std::cout<<"(1) evidence variable: observing var["<<curr_var<<"] == "<<curr_val<<";"<<std::endl; }
                // We will account for evidence by adding a new single node potential.
                all_factors.push_back(std::vector<MPLPIndexType>(1, curr_var));
                std::vector<double> l_i = std::vector<double>(var_sizes[curr_var], -MPLP_huge);
                l_i[curr_val] = 0.0;
                all_lambdas.push_back(l_i);
            }
        }
    }
    return nevid;
}


// TODO: I believe there is a bug here: when a factor has a single variable and it is observed as evidence, an empty factor is created
mplpLib::MPLPIndexType mplpLib::read_model_file(std::vector<MPLPIndexType> & var_sizes, std::map<MPLPIndexType, MPLPIndexType> & evidence, std::vector< std::vector<MPLPIndexType> > & all_factors, std::vector< std::vector<double> > & all_lambdas, const std::string fn, const std::string evid_fn){
    MPLPIndexType i, j, v, nvars, nfactors, nevid, evid_size, curr_var, curr_val;
    double val;
    std::vector< std::vector<MPLPIndexType> > f;
    std::vector< std::vector<double> > l;
    std::string s, t;
    std::ifstream fstr;
    if (evid_fn != ""){
        fstr.open(evid_fn.c_str(), std::fstream::in);
        if (fstr.fail()){
            std::cerr<<"Error opening evidence file"<<std::endl;
        }else{
            fstr>>nevid;   //for PIC we have at most 1 evidence
            if (nevid){
                fstr>>evid_size;
                while (evid_size--){
                    fstr>>curr_var;
                    fstr>>curr_val;
                    if (MPLP_DEBUG_MODE){ std::cout<<"(2) evidence variable: observing var["<<curr_var<<"] == "<<curr_val<<";"<<std::endl; }
                    evidence[curr_var] = curr_val;
                }
            }
        }
        fstr.close();
    }
    fstr.open(fn.c_str(), std::fstream::in);
    if (fstr.fail()){
        std::cerr<<"Error opening model file"<<std::endl;
        return nevid;
    }
    fstr>>s;
    if (s != "MARKOV"){
        std::cerr<<"Error: file is not in MARKOV format."<<std::endl;
    }
    fstr>>nvars;
    if (MPLP_DEBUG_MODE){ std::cout<<"nvars = "<<nvars<<std::endl; }
    for (i = 0; i < nvars; ++i){
        fstr>>v;
        var_sizes.push_back(v);
        if (MPLP_DEBUG_MODE){ std::cout<<"var_sizes["<<i<<"] = "<<var_sizes[i]<<std::endl; }
    }
    fstr>>nfactors;
    if (MPLP_DEBUG_MODE){ std::cout<<"nfactors = "<<nfactors<<std::endl; }
    for (i = 0; i < nfactors; ++i){
        fstr>>nvars;
        f.push_back(std::vector<MPLPIndexType>());
        all_factors.push_back(std::vector<MPLPIndexType>());
        l.push_back(std::vector<double>());
        if (MPLP_DEBUG_MODE){ std::cout<<"factor_size["<<i<<"] = "<<nvars<<std::endl; }
        for (j = 0; j < nvars; ++j){
            fstr>>curr_var;  //those are the cliques
            f[i].push_back(curr_var);
            if (MPLP_DEBUG_MODE){ std::cout<<"all_factors["<<i<<"]["<<j<<"] = "<<f[i][j]<<std::endl; }
            if (evidence.find(curr_var) == evidence.end()){
                all_factors[i].push_back(curr_var);
            }
        }
    }//"all_lambdas": the lambdas (messages) are just log(function table entry) 's
    for (i = 0; i < nfactors; ++i){
        MPLPIndexType fact = 1, base = 0;
        std::vector<MPLPIndexType> curr_evid;   //list of evidence node contained in a factor (ordered accordingly)
        std::vector<double> lambdas;
        fstr>>nvars;  //the PIC file does not contain empty lines
        for (j = 0; j < nvars; ++j){
            fstr>>val;
            // We work in log space, so take the log of the factors' potentials
            l[i].push_back(val > 0 ? log(val) : -MPLP_huge);
            if (MPLP_DEBUG_MODE){ std::cout<<"e^all_lambdas["<<i<<"]["<<j<<"] = "<<val<<std::endl; }
            if (MPLP_DEBUG_MODE){ std::cout<<"all_lambdas["<<i<<"]["<<j<<"] = "<<l[i][j]<<std::endl; }
        }
        for (j = 0; j < f[i].size(); ++j){
            if (evidence.find(f[i][j]) != evidence.end()){
                curr_evid.push_back(f[i][j]);
            }
        }
        while (--j){
            if (evidence.find(f[i][j]) != evidence.end()){
                base += evidence[f[i][j]] * fact;
            }
            //std::cout<<"nvars == "<<nvars<<", i == "<<i<<", j == "<<j<<", f[i][j] == "<<f[i][j]<<std::endl;
            fact *= var_sizes[f[i][j]];
        }
        if (evidence.find(f[i][0]) != evidence.end()){
            base += evidence[f[i][0]] * fact;
        }
        get_lambdas(base, fact, 0, 0, var_sizes, evidence, curr_evid, f[i], l[i], lambdas);
        all_lambdas.push_back(lambdas);

        if(MPLP_DEBUG_MODE) {
            if(lambdas.size() != l[i].size()) {
                std::cout << "THERE WAS EVIDENCE AND WE SHRUNK THE FACTOR! New size: " << lambdas.size() << ", Old size: " << l[i].size() << std::endl;
            }
        }

    }
    fstr.close();
    return nevid;
}

mplpLib::MPLPIndexType mplpLib::read_model_file2(std::vector<MPLPIndexType> & var_sizes, std::map<MPLPIndexType, MPLPIndexType> & evidence, std::vector< std::vector<MPLPIndexType> > & all_factors, std::vector< std::vector<double> > & all_lambdas, const std::string fn, const std::string evid_fn){
    MPLPIndexType i, j, v, nvars, nfactors, nevid, evid_size, curr_var, curr_val;
    double val;
    std::vector< std::vector<MPLPIndexType> > f;
    std::vector< std::vector<double> > l;
    std::string s, t;
    std::ifstream fstr;
    if (evid_fn != ""){
        fstr.open(evid_fn.c_str(), std::fstream::in);
        if (fstr.fail()){
            std::cerr<<"Error opening evidence file"<<std::endl;
        }else{
            fstr>>nevid;   //for PIC we have at most 1 evidence
            if (nevid){
                fstr>>evid_size;
                while (evid_size--){
                    fstr>>curr_var;
                    fstr>>curr_val;
                    if (MPLP_DEBUG_MODE){ std::cout<<"evidence variable: observing var["<<curr_var<<"] == "<<curr_val<<";"<<std::endl; }
                    evidence[curr_var] = curr_val;
                }
            }
        }
        fstr.close();
    }
    fstr.open(fn.c_str(), std::fstream::in);
    if (fstr.fail()){
        std::cerr<<"Error opening model file"<<std::endl;
        return nevid;
    }
    fstr>>s;
    if (s != "MARKOV"){
        std::cerr<<"Error: file is not in MARKOV format."<<std::endl;
    }
    fstr>>nvars;
    if (MPLP_DEBUG_MODE){ std::cout<<"nvars = "<<nvars<<std::endl; }
    for (i = 0; i < nvars; ++i){
        fstr>>v;
        var_sizes.push_back(v);
        if (MPLP_DEBUG_MODE){ std::cout<<"var_sizes["<<i<<"] = "<<var_sizes[i]<<std::endl; }
    }
    fstr>>nfactors;
    if (MPLP_DEBUG_MODE){ std::cout<<"nfactors = "<<nfactors<<std::endl; }
    for (i = 0; i < nfactors; ++i){
        fstr>>nvars;
        f.push_back(std::vector<MPLPIndexType>());
        all_factors.push_back(std::vector<MPLPIndexType>());
        l.push_back(std::vector<double>());
        if (MPLP_DEBUG_MODE){ std::cout<<"factor_size["<<i<<"] = "<<nvars<<std::endl; }
        for (j = 0; j < nvars; ++j){
            fstr>>curr_var;  //those are the cliques
            f[i].push_back(curr_var);
            if (MPLP_DEBUG_MODE){ std::cout<<"all_factors["<<i<<"]["<<j<<"] = "<<f[i][j]<<std::endl; }
            if (evidence.find(curr_var) == evidence.end()){
                all_factors[i].push_back(curr_var);
            }
        }
    }//"all_lambdas": the lambdas (messages) are just log(function table entry) 's
    for (i = 0; i < nfactors; ++i){
        MPLPIndexType fact = 1, base = 0;
        std::vector<MPLPIndexType> curr_evid;   //list of evidence node contained in a factor (ordered accordingly)
        std::vector<double> lambdas;
        int temp_nvars;
        fstr>>temp_nvars;  //the PIC file does not contain empty lines

        // Potts model potential denoted by -1. Assume that all variables have same number of states!
        if(temp_nvars == -1) {
            fstr>>val;

            for (MPLPIndexType state1 = 0; state1 < var_sizes[0]; ++state1) {
                for (MPLPIndexType state2 = 0; state2 < var_sizes[0]; ++state2) {
                    l[i].push_back(state1 == state2 ? 0 : val);
                }
            }
        }
        else {
            nvars = static_cast<MPLPIndexType>(temp_nvars);
            for (j = 0; j < nvars; ++j){
                fstr>>val;
                // We work in log space, so take the log of the factors' potentials
                //l[i].push_back(val > 0 ? log(val) : -MPLP_huge);
                l[i].push_back(val);
                if (MPLP_DEBUG_MODE){ std::cout<<"e^all_lambdas["<<i<<"]["<<j<<"] = "<<val<<std::endl; }
                if (MPLP_DEBUG_MODE){ std::cout<<"all_lambdas["<<i<<"]["<<j<<"] = "<<l[i][j]<<std::endl; }
            }
        }

        for (j = 0; j < f[i].size(); ++j){
            if (evidence.find(f[i][j]) != evidence.end()){
                curr_evid.push_back(f[i][j]);
            }
        }
        while (--j){
            if (evidence.find(f[i][j]) != evidence.end()){
                base += evidence[f[i][j]] * fact;
            }
            //std::cout<<"nvars == "<<nvars<<", i == "<<i<<", j == "<<j<<", f[i][j] == "<<f[i][j]<<std::endl;
            fact *= var_sizes[f[i][j]];
        }
        if (evidence.find(f[i][0]) != evidence.end()){
            base += evidence[f[i][0]] * fact;
        }
        get_lambdas(base, fact, 0, 0, var_sizes, evidence, curr_evid, f[i], l[i], lambdas);
        all_lambdas.push_back(lambdas);
    }
    fstr.close();
    return nevid;
}

void mplpLib::get_lambdas(MPLPIndexType base, MPLPIndexType fact, MPLPIndexType ind, MPLPIndexType p, const std::vector<MPLPIndexType> & var_sizes, const std::map<MPLPIndexType, MPLPIndexType> & evidence, const std::vector<MPLPIndexType> & curr_evid, const std::vector<MPLPIndexType> & f, const std::vector<double> & l, std::vector<double> & lambdas){
    while (ind < f.size()){
        MPLPIndexType _fact = (ind == f.size() - 1 ? 0 : fact / var_sizes[f[ind + 1]]);
        if (p >= curr_evid.size() || f[ind] != curr_evid[p]){  //not in the evidence set
            for (MPLPIndexType i = 0; i < var_sizes[f[ind]]; ++i){
                get_lambdas(base + i * fact, _fact, ind + 1, p, var_sizes, evidence, curr_evid, f, l, lambdas);
            }
            return;
        }else{
            fact = _fact;
            ++ind;
            ++p;
        }
    }
    lambdas.push_back(l[base]);
}

/*
int main(int argc, char** argv){
	std::vector<int> var_sizes;
	std::vector< std::vector<int> > all_factors;
	std::vector< std::vector<double> > all_lambdas;
	if (argc == 2){
		read_model_file(var_sizes, all_factors, all_lambdas, argv[1]);   //should replace with meaningful parameters
		//read_model_file(NULL, NULL, NULL, "", 0, "");  //optional args
	}else{
		read_model_file(var_sizes, all_factors, all_lambdas, "grid3x3.uai", "grid3x3.uai.evid");
	}
	return 0;
}
 */
