#pragma once
#include <iostream>
#include <fstream>
#include "events.h"
#include <vector>
#include <sstream>
#include <string>
#include "inparams.h"
#include "reso_decays.h"

using std::istringstream ; 
class read_input_file{
  public :

    read_input_file(input_paramters &iparam_, reso_decays*) ;
    ~read_input_file();
    void read_input_file_iSS_OSCAR(int );
    void read_particle_list_dat_from_urqmd(int );
    void read_particle_list_dat_from_urqmd_binary(int TotalEvents) ; 

    events* get_event(int xx){
      return &event_vector[xx] ; 
    }

  private :
    istringstream* iss;
    input_paramters &iparam ; 
    reso_decays* resonance_decays ; 
    char buff[400];
    std::vector<events> event_vector ;
    int get_PID_from_urqmd_MCID(int mcid, int iso) ; 


};
