#include<iostream>
#include<fstream>
#include <vector>
#include "events.h"
#include "read_input_file.h"
#include "observables.h"
#include "inparams.h"
#include "read_pdg.h"
#include "inparams.h"
#include "reconstruction.h"


int main(int argc, char **argv){

  if(argc != 3){
    std::cout << "2 arguments required ..." << std::endl ;
    std::cout << "path of files and number of files to be analysed ... " << std::endl ; 
    exit(1); 
  }

  input_paramters iparams ;
  read_parameters* r = new read_parameters();
  r->read_parameters_from_file(iparams, "input_parameters"); 

  read_pdg* RPDG = new read_pdg(iparams); 
  if(iparams.pdg_type==1){
    RPDG->read_and_store_particle_properties_with_decay_channels("PDG/pdg-urqmd_v3.3+_weak.dat");
  }
  else{
    RPDG->read_and_store_particle_properties_with_decay_channels("PDG/pdg05_music_wd_off.dat");
  }

  reso_decays* RD = new reso_decays(RPDG);
  read_input_file* RIF = new read_input_file(iparams, RD, argv[1], atof(argv[2]) );

  int nEvents = iparams.nEvents ; 
  
  // reading the input files //
  if(iparams.input_read_mode == 0 ){
    RIF->read_input_file_iSS_OSCAR(nEvents); }
  else if(iparams.input_read_mode == 1 ){
    RIF->read_particle_list_dat_from_urqmd(nEvents); }
  else if(iparams.input_read_mode == 2 ){
    RIF->read_particle_list_dat_from_urqmd_binary(nEvents); }
  else if(iparams.input_read_mode == 3 ){
    RIF->read_particle_list_dat_from_iSS_binary(nEvents); }
  else {
    std::cout << "reading mode not specified. Exiting ..." << std::endl ;  
    exit(1); 
  }

  // Reconstructing the resonances //
  reconstruction* REC = new reconstruction(iparams,RIF);
  if(iparams.reconstruct_phi_flag > 0){
    std::cout << "reconstructing phi meson" 
              << " by invariant mass method ..." 
              << std::endl ; 
    REC->reconstruct_phi();
  }
  if(iparams.reconstruct_kstar0_flag > 0){
    std::cout << "reconstructing kstar-0 meson" 
              << " by invariant mass method ..." 
              << std::endl ; 
    REC->reconstruct_kstar0();
  }
  if(iparams.reconstruct_kstar0_bar_flag > 0){
    std::cout << "reconstructing kstar-0-bar meson" 
              << " by invariant mass method ..." 
              << std::endl ; 
    REC->reconstruct_kstar0_bar();
  }


  // calculating the observables //
  observables* OBJ = new observables(iparams,RIF);
  //OBJ->calculate_dnchdeta_eta(0.01,3);
  //OBJ->calculate_dndy_y(0.01,3);
  //OBJ->calculate_invariant_yield_vs_pt(0, -0.5, 0.5);
  //OBJ->calculate_invariant_yield_vs_pt(1, -0.5, 0.5);
  //OBJ->calculate_v1_vs_y_or_eta(0, 0, 0.2, 2 );
  //OBJ->calculate_v1_vs_y_or_eta(1, 0, 0.2, 2 );
  //OBJ->calculate_v1_vs_y_or_eta(0, 0, 0.4, 2 );
  //OBJ->calculate_v1_vs_y_or_eta(1, 0, 0.4, 2 );
  //OBJ->calculate_v2_pt( 0, -1.0, 1.0 );
  //OBJ->calculate_v2_pt( 1, -0.5, 0.5 );

  //OBJ->calculate_amn(  211,  0,  -0.5,  0.5, 0.1, 2);
  //OBJ->calculate_amn_from_smeared_grid(  211,  0,  -0.5,  0.5, 0.1, 2);
  //OBJ->calculate_amn_of_charged_hadrons( 0,  -0.5,  0.5, 0.1, 2);
  //OBJ->calculate_amn_of_charged_hadrons( 0,  -0.5,  0.5, 0.01, 2);
  OBJ->calculate_amn_vs_rapidity( 211, 1,  0.2, 2);
  OBJ->calculate_amn_vs_rapidity( 321, 1,  0.2, 2);
  OBJ->calculate_amn_vs_rapidity( -321, 1,  0.2, 2);
  OBJ->calculate_amn_vs_rapidity( 2212, 1,  0.2, 2);
  OBJ->calculate_amn_vs_rapidity( -2212, 1,  0.2, 2);
  return 0;
}





