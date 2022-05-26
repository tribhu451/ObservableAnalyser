#pragma once
#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TProfile.h"
#include "TMath.h"
#include "events.h"
#include "read_input_file.h"
#include "inparams.h"

class observables{

  public :

    observables(input_paramters &iparam_, read_input_file* );
    void calculate_dnchdeta_eta(double pT_min, double pT_max);
    void calculate_dndy_y( double pT_min, double pT_max);
    void calculate_invariant_yield_vs_pt( int yflag, double Rap_min, double Rap_max );
    void calculate_v1_vs_y_or_eta(int yflag, double psi1,  double pT_min, double pT_max );
    void calculate_v2_vs_y_or_eta(int yflag, double psi1,  double pT_min, double pT_max );
    void calculate_v2_pt( int yflag, double Rap_min, double Rap_max );


  private :

    input_paramters &iparam ; 
    read_input_file* rif ;
    double fit_a_straight_line_and_get_slope(int n, double *x, double *y) ;     

};
