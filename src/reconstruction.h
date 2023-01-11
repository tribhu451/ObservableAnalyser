#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include "TMath.h"
#include "events.h"
#include "read_input_file.h"
#include "inparams.h"

class reconstruction{

  public :

    reconstruction(input_paramters &iparam_, read_input_file* );
    ~reconstruction();
    void reconstruct_phi();
    void reconstruct_kstar0();
    void reconstruct_kstar0_bar();

  private :

    input_paramters &iparam ; 
    read_input_file* rif ;

};









