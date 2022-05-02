#include<iostream>
#include<fstream>
#include <vector>
#include "events.h"
#include "read_input_file.h"
#include "observables.h"
#include "inparams.h"
#include "read_pdg.h"
#include "inparams.h"


int main(){


  input_paramters iparams ;
  read_parameters* r = new read_parameters();
  r->read_parameters_from_file(iparams, "input_parameters"); 

  read_pdg* RPDG = new read_pdg(iparams); 
  RPDG->read_and_store_particle_properties_with_decay_channels("PDG/pdg-urqmd_v3.3+_weak.dat");

  reso_decays* RD = new reso_decays(RPDG);
  read_input_file* RIF = new read_input_file(iparams, RD);

  int nEvents = iparams.nEvents ; 
  
  if(iparams.input_read_mode == 0 ){
    RIF->read_input_file_iSS_OSCAR(nEvents); }
  else if(iparams.input_read_mode == 1 ){
    RIF->read_particle_list_dat_from_urqmd(nEvents); }
  else if(iparams.input_read_mode == 2 ){
    RIF->read_particle_list_dat_from_urqmd_binary(nEvents); }
  else {
    std::cout << "reading mode not specified. Exiting ..." << std::endl ;  
    exit(1); }


  observables* OBJ = new observables(iparams,RIF);
  OBJ->calculate_dnchdeta_eta(0.01,3);
  OBJ->calculate_dndy_y(0.01,3);
  OBJ->calculate_invariant_yield_vs_pt(0, -0.5, 0.5);
  OBJ->calculate_invariant_yield_vs_pt(1, -0.5, 0.5);
  OBJ->calculate_v1_vs_y_or_eta(0, 0, 0.1, 3 );
  OBJ->calculate_v1_vs_y_or_eta(1, 0, 0.1, 3 );
  OBJ->calculate_v2_pt( 0, -0.5, 0.5 );
  OBJ->calculate_v2_pt( 1, -0.5, 0.5 );
  /*
  for(int ii=0; ii<nEvents; ii++){
    OBJ->fill_histogram_of_dnchdeta_eta(RIF->get_event(ii), 0.01, 3.0 );
    OBJ->fill_histogram_of_v2_pt(RIF->get_event(ii), -0.5, 0.5);
    OBJ->fill_histogram_of_v1_y(RIF->get_event(ii), 0.01, 3.0 );
    OBJ->fill_histogram_of_v1_eta(RIF->get_event(ii), 0.01, 3.0 );
    OBJ->fill_histogram_of_invariant_yield_pt(RIF->get_event(ii), -0.5 ,0.5) ;
  }

   OBJ->calculate_dnchdeta_eta(nEvents);
   OBJ->calculate_v2_pt();
   OBJ->calculate_v1_eta();
   OBJ->calculate_v1_y();
   OBJ->calculate_invariant_yield_pt(nEvents, -0.5, 0.5) ;
   */

  return 0;
}





