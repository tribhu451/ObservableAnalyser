#pragma once
#include<iostream>
#include<cmath>
#include "read_input_file.h"
#include "random.h"
#include "inparams.h"
#include <vector>


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
    float* eta_bin_center ; 
    float* IptI ; // single event mean-pt : [pt]
    inline float get_IptI(int ievents, int ieta){ return IptI[ievents*Neta+ieta];}
    inline void set_IptI(int ievents, int ieta, float xx){ IptI[ievents*Neta+ieta]=xx;}

   public :
     mpt_decorr(input_paramters &iparam_, read_input_file*, 
      int _part, int _yflag, float _ptmin, float _ptmax );
     ~mpt_decorr();

};
