#pragma once
#include <iostream>
#include <cmath>
#include "read_input_file.h"
#include "random.h"
#include "inparams.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <string>


class mpt_decorr{

   private :
    input_paramters &iparam ; 
    read_input_file* rif ;
    random_gen* rand;
    int Nevents;
    int part ;
    int yflag ;
    float ptmin ;
    float ptmax ; 
    float etamin ; 
    float etamax ; 
    int Neta ; 
    float deta ; 
    float* eta_bin_centers ; 
    float* IptI ; // single event mean-pt : [pt]
    inline float get_IptI(int ievents, int ieta){ return IptI[ievents*Neta+ieta];}
    inline void set_IptI(int ievents, int ieta, float xx){ IptI[ievents*Neta+ieta]=xx;}
    inline int get_rap_bin_index(double eta_val){ return (eta_val - etamin) / deta ;}
    std::vector<int> get_an_event_ensemble();
    void calculate_mean_and_variance_of_pt_of_one_ensemble( int ieta,
      std::vector<int> event_ID_ens, double&  Mpt, double&  var_Mpt);
    void calculate_covariance_of_mean_pt_of_one_ensemble(int ieta1, int ieta2, 
      std::vector<int> event_ID_ens, double& cov, double& RMpt);
    void calculate_r_mean_pt_of_one_ensemble(int ieta1, int ieta2, 
      std::vector<int> event_ID_ens, double& rMpt);

   public :
     mpt_decorr(input_paramters &iparam_, read_input_file*, 
      int _part, int _yflag, float _ptmin, float _ptmax );
     ~mpt_decorr();
     void write_mean_and_variance_of_pt_with_eta();
     void write_covariance_of_meanpt();

};
