#include "mpt_decorr.h"
#include "random.h"

mpt_decorr::mpt_decorr(input_paramters &iparam_, read_input_file* rif_, 
 int _part, int _yflag, float _ptmin, float _ptmax) : iparam(iparam_), rif(rif_), 
 part(_part), yflag(_yflag), ptmin(_ptmin), ptmax(_ptmax){
  rand = new random_gen();
  Nevents = rif->get_event_buffer_size() ; 

  etamin = -5.25 ; 
  etamax = 5.25 ;
  Neta = 21 ;  
  const float deta = ( etamax - etamin ) / Neta ; 
  eta_bin_center = new float[Neta];
  for (int ii = 0; ii < Neta; ii++) {
     eta_bin_center[ii] = etamin + (ii+0.5) * deta ;
     std::cout << etamin + (ii) * deta  << "   " << eta_bin_center[ii] << "   " <<  etamin + (ii+1) * deta  << std::endl ; 
   }

  // storing [pt] in event x EtaBin matrix 
  IptI = new float[Nevents*Neta]; // (i-event, j-etabin)

  // calculate meanpt of each event at each rapidity window and store it.

}


mpt_decorr::~mpt_decorr(){
  delete[] eta_bin_center ; 
  delete[] IptI ; 
}
