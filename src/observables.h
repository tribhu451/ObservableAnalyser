#pragma once
#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TProfile.h"
#include "TMath.h"
#include "events.h"


#define _N_HISTOGRAMS_DNCHDETA_ETA 7
#define _N_HISTOGRAMS_V2_PT 5
#define _N_HISTOGRAMS_V1_Y 7
#define _N_HISTOGRAMS_V1_ETA 1
#define _N_HISTOGRAMS_SPECTRA 6

class observables{

  public :
  observables();
  ~observables();
  void fill_histogram_of_dnchdeta_eta(events* , double, double);
  void calculate_dnchdeta_eta(int );
  void fill_histogram_of_v2_pt(events* , double , double );
  void calculate_v2_pt();
  void calculate_v1_eta();
  void fill_histogram_of_v1_eta(events* Event, double pT_min, double pT_max );
  void calculate_v1_y();
  void fill_histogram_of_v1_y(events* Event, double pT_min, double pT_max );
  void fill_histogram_of_invariant_yield_pt(events* Event, double eta_min, double eta_max ) ;
  void calculate_invariant_yield_pt(int Nevents, double eta_min, double eta_max);

  private :
   
   float eta_min ;
   float eta_max ;
   int   eta_bins ;
   TH1D*       H1D_DNCHDETA_ETA[_N_HISTOGRAMS_DNCHDETA_ETA] ; 
   TProfile*   PROFILE_V2_PT[_N_HISTOGRAMS_V2_PT] ; 
   TProfile*   PROFILE_V1_Y[_N_HISTOGRAMS_V1_Y] ; 
   TProfile*   PROFILE_V1_ETA[_N_HISTOGRAMS_V1_ETA] ; 
   TH1D*       H1D_INVYLD_PT[_N_HISTOGRAMS_SPECTRA] ; 

   float pT_min   ; 
   float pT_max   ; 
   int   pT_bins  ; 

   float y_min ;
   float y_max ;
   int   y_bins ;

};
