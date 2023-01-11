#include "reconstruction.h"

reconstruction::reconstruction(input_paramters &iparam_, read_input_file* rif_) : iparam(iparam_), rif(rif_) {
}

reconstruction::~reconstruction(){
}


void reconstruction::reconstruct_phi(){
  double phi_meson_mass  = 1.01899999 ;
  double phi_meson_width = 0.0001 ;

  phi_meson_mass = iparam.reconst_phi_mass ;
  phi_meson_width = iparam.reconst_phi_width ;

  std::cout << "reconstructing assuming phi_meson_mass = " 
            << phi_meson_mass << std::endl ; 
  std::cout << "reconstructing assuming phi_meson_width = " 
            << phi_meson_width << std::endl ; 


  const int maxpart = 5000;
  int kaon_plus_index[maxpart];
  int kaon_minus_index[maxpart];
  int k0_index[maxpart];
  int k0_bar_index[maxpart];

  int number_of_kaon_plus   = 0 ; 
  int number_of_kaon_minus  = 0 ; 
  int number_of_k0          = 0 ; 
  int number_of_k0_bar      = 0 ; 

  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){

    number_of_kaon_plus    = 0 ; 
    number_of_kaon_minus   = 0 ; 
    number_of_k0           = 0 ; 
    number_of_k0_bar       = 0 ; 

    for(int jj=0; jj<maxpart; jj++){
      kaon_plus_index[jj]   = -9999 ; 
      kaon_minus_index[jj]  = -9999 ; 
      k0_index[jj]          = -9999 ; 
      k0_bar_index[jj]      = -9999 ; 
    }


    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
       int    PID = Event->get_particle(jj)->get_pid() ; 
       if(PID == 321 ){
         kaon_plus_index[number_of_kaon_plus] = jj ;
         number_of_kaon_plus  += 1 ; 
       }
       if(PID == -321 ){
         kaon_minus_index[number_of_kaon_minus] = jj ;
         number_of_kaon_minus += 1 ; 
       }
       if(PID == 311 ){
         k0_index[number_of_k0] = jj ;
         number_of_k0 += 1 ; 
       }
       if(PID == -311 ){
         k0_bar_index[number_of_k0_bar] = jj ;
         number_of_k0_bar += 1 ; 
       }
    } // particle loop

     // check for all k^+ + k^- channels
     for(int ij=0; ij<number_of_kaon_plus; ij++){
       for(int ik=0; ik<number_of_kaon_minus; ik++){
         double Px1  = Event->get_particle(kaon_plus_index[ij])->get_px()  ; 
         double Py1  = Event->get_particle(kaon_plus_index[ij])->get_py()  ; 
         double Pz1  = Event->get_particle(kaon_plus_index[ij])->get_pz()  ; 
         double E1   = Event->get_particle(kaon_plus_index[ij])->get_e()   ; 
         double t1   = Event->get_particle(kaon_plus_index[ij])->get_t()   ; 
         double x1   = Event->get_particle(kaon_plus_index[ij])->get_x()   ; 
         double y1   = Event->get_particle(kaon_plus_index[ij])->get_y()   ; 
         double z1   = Event->get_particle(kaon_plus_index[ij])->get_z()   ; 

         double Px2  = Event->get_particle(kaon_minus_index[ik])->get_px()  ; 
         double Py2  = Event->get_particle(kaon_minus_index[ik])->get_py()  ; 
         double Pz2  = Event->get_particle(kaon_minus_index[ik])->get_pz()  ; 
         double E2   = Event->get_particle(kaon_minus_index[ik])->get_e()   ; 
         double t2   = Event->get_particle(kaon_minus_index[ik])->get_t()   ; 
         double x2   = Event->get_particle(kaon_minus_index[ik])->get_x()   ; 
         double y2   = Event->get_particle(kaon_minus_index[ik])->get_y()   ; 
         double z2   = Event->get_particle(kaon_minus_index[ik])->get_z()   ; 

         // compute the invariant mass
         double E = E1 + E2;
         double Px = Px1 + Px2;
         double Py = Py1 + Py2;
         double Pz = Pz1 + Pz2;
         double invariant_mass = sqrt(E*E - Px*Px - Py*Py - Pz*Pz);

         if (std::abs(invariant_mass - phi_meson_mass)// phi(1020) resonance found
                                      < phi_meson_width) {
            double spatial_distance = sqrt(
                 (t1 - t2)*(t1 - t2) + (x1 - x2)*(x1 - x2)
                   + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)) ;
             if (spatial_distance < 0.001) {
                 Event->add_particle(333, t1, x1, y1, z1, E, Px, Py, Pz,true, 1) ;
             }
         } // if phi found
   
       } // loop over kaon minus
     } // loop over kaon plus


     // check for all k^0 + k^0-bar channels
     for(int ij=0; ij<number_of_k0; ij++){
       for(int ik=0; ik<number_of_k0_bar; ik++){
         double Px1  = Event->get_particle(k0_index[ij])->get_px()  ; 
         double Py1  = Event->get_particle(k0_index[ij])->get_py()  ; 
         double Pz1  = Event->get_particle(k0_index[ij])->get_pz()  ; 
         double E1   = Event->get_particle(k0_index[ij])->get_e()   ; 
         double t1   = Event->get_particle(k0_index[ij])->get_t()   ; 
         double x1   = Event->get_particle(k0_index[ij])->get_x()   ; 
         double y1   = Event->get_particle(k0_index[ij])->get_y()   ; 
         double z1   = Event->get_particle(k0_index[ij])->get_z()   ; 

         double Px2  = Event->get_particle(k0_bar_index[ik])->get_px()  ; 
         double Py2  = Event->get_particle(k0_bar_index[ik])->get_py()  ; 
         double Pz2  = Event->get_particle(k0_bar_index[ik])->get_pz()  ; 
         double E2   = Event->get_particle(k0_bar_index[ik])->get_e()   ; 
         double t2   = Event->get_particle(k0_bar_index[ik])->get_t()   ; 
         double x2   = Event->get_particle(k0_bar_index[ik])->get_x()   ; 
         double y2   = Event->get_particle(k0_bar_index[ik])->get_y()   ; 
         double z2   = Event->get_particle(k0_bar_index[ik])->get_z()   ; 

         // compute the invariant mass
         double E = E1 + E2;
         double Px = Px1 + Px2;
         double Py = Py1 + Py2;
         double Pz = Pz1 + Pz2;
         double invariant_mass = sqrt(E*E - Px*Px - Py*Py - Pz*Pz);

         if (std::abs(invariant_mass - phi_meson_mass)// phi(1020) resonance found
                                      < phi_meson_width) {
            double spatial_distance = sqrt(
                 (t1 - t2)*(t1 - t2) + (x1 - x2)*(x1 - x2)
                   + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)) ;
             if (spatial_distance < 0.001) {
                 Event->add_particle(333, t1, x1, y1, z1, E, Px, Py, Pz,true,1) ;
             }
         } // if phi found
   
       } // loop over k0
     } // loop over k0-bar

  } // event loop

}



void reconstruction::reconstruct_kstar0(){

  double kstar0_mass  = 0.893 ;
  //double kstar0_width = 0.050 ;
  double kstar0_width = 0.00001 ;

  kstar0_mass  = iparam.reconst_kstar0_mass ; 
  kstar0_width = iparam.reconst_kstar0_width ; 

  std::cout << "reconstructing assuming kstar0_mass = " 
            << kstar0_mass << std::endl ; 
  std::cout << "reconstructing assuming kstar0_width = " 
            << kstar0_width << std::endl ; 


  const int maxpart = 5000;
  int pion_minus_index[maxpart];
  int kaon_plus_index[maxpart];

  int number_of_pion_minus  = 0 ; 
  int number_of_kaon_plus   = 0 ; 

  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){

    number_of_pion_minus   = 0 ; 
    number_of_kaon_plus    = 0 ;

    for(int jj=0; jj<maxpart; jj++){
      pion_minus_index[jj]  = -9999 ; 
      kaon_plus_index[jj]   = -9999 ; 
    }


    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
       int    PID = Event->get_particle(jj)->get_pid() ; 
       if(PID == -211 ){
         pion_minus_index[number_of_pion_minus] = jj ;
         number_of_pion_minus += 1 ; 
       }
       if(PID == 321 ){
         kaon_plus_index[number_of_kaon_plus] = jj ;
         number_of_kaon_plus  += 1 ; 
       }
    } // particle loop

     // check for all k^+ + pi^- channels
     for(int ij=0; ij<number_of_kaon_plus; ij++){
       for(int ik=0; ik<number_of_pion_minus; ik++){
         double Px1  = Event->get_particle(kaon_plus_index[ij])->get_px()  ; 
         double Py1  = Event->get_particle(kaon_plus_index[ij])->get_py()  ; 
         double Pz1  = Event->get_particle(kaon_plus_index[ij])->get_pz()  ; 
         double E1   = Event->get_particle(kaon_plus_index[ij])->get_e()   ; 
         double t1   = Event->get_particle(kaon_plus_index[ij])->get_t()   ; 
         double x1   = Event->get_particle(kaon_plus_index[ij])->get_x()   ; 
         double y1   = Event->get_particle(kaon_plus_index[ij])->get_y()   ; 
         double z1   = Event->get_particle(kaon_plus_index[ij])->get_z()   ; 

         double Px2  = Event->get_particle(pion_minus_index[ik])->get_px()  ; 
         double Py2  = Event->get_particle(pion_minus_index[ik])->get_py()  ; 
         double Pz2  = Event->get_particle(pion_minus_index[ik])->get_pz()  ; 
         double E2   = Event->get_particle(pion_minus_index[ik])->get_e()   ; 
         double t2   = Event->get_particle(pion_minus_index[ik])->get_t()   ; 
         double x2   = Event->get_particle(pion_minus_index[ik])->get_x()   ; 
         double y2   = Event->get_particle(pion_minus_index[ik])->get_y()   ; 
         double z2   = Event->get_particle(pion_minus_index[ik])->get_z()   ; 

         // compute the invariant mass
         double E = E1 + E2;
         double Px = Px1 + Px2;
         double Py = Py1 + Py2;
         double Pz = Pz1 + Pz2;
         double invariant_mass = sqrt(E*E - Px*Px - Py*Py - Pz*Pz);

         if (std::abs(invariant_mass - kstar0_mass)// k*0(892) resonance found
                                      < kstar0_width) {
            double spatial_distance = sqrt(
                 (t1 - t2)*(t1 - t2) + (x1 - x2)*(x1 - x2)
                   + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)) ;
             if (spatial_distance < 0.001) {
                 Event->add_particle(313, t1, x1, y1, z1, E, Px, Py, Pz,true,1) ;
             }
         } // if k*0(892) found
   
       } // loop over kaon plus
     } // loop over pion minus

  } // event loop

}



void reconstruction::reconstruct_kstar0_bar(){

  double kstar0_bar_mass  = 0.893 ;
  double kstar0_bar_width = 0.00001 ;

  kstar0_bar_mass  = iparam.reconst_kstar0_bar_mass ; 
  kstar0_bar_width = iparam.reconst_kstar0_bar_width ; 

  std::cout << "reconstructing assuming kstar0_bar_mass = " 
            << kstar0_bar_mass << std::endl ; 
  std::cout << "reconstructing assuming kstar0_bar_width = " 
            << kstar0_bar_width << std::endl ; 


  const int maxpart = 5000;
  int pion_plus_index[maxpart];
  int kaon_minus_index[maxpart];

  int number_of_pion_plus   = 0 ; 
  int number_of_kaon_minus  = 0 ; 

  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){

    number_of_pion_plus    = 0 ; 
    number_of_kaon_minus   = 0 ; 

    for(int jj=0; jj<maxpart; jj++){
      pion_plus_index[jj]   = -9999 ; 
      kaon_minus_index[jj]  = -9999 ; 
    }


    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
       int    PID = Event->get_particle(jj)->get_pid() ; 

       if(PID == 211 ){
         pion_plus_index[number_of_pion_plus] = jj ;
         number_of_pion_plus  += 1 ; 
       }
       if(PID == -321 ){
         kaon_minus_index[number_of_kaon_minus] = jj ;
         number_of_kaon_minus += 1 ; 
       }
    } // particle loop

     // check for all k^- + pi^+ channels
     for(int ij=0; ij<number_of_kaon_minus; ij++){
       for(int ik=0; ik<number_of_pion_plus; ik++){
         double Px1  = Event->get_particle(kaon_minus_index[ij])->get_px()  ; 
         double Py1  = Event->get_particle(kaon_minus_index[ij])->get_py()  ; 
         double Pz1  = Event->get_particle(kaon_minus_index[ij])->get_pz()  ; 
         double E1   = Event->get_particle(kaon_minus_index[ij])->get_e()   ; 
         double t1   = Event->get_particle(kaon_minus_index[ij])->get_t()   ; 
         double x1   = Event->get_particle(kaon_minus_index[ij])->get_x()   ; 
         double y1   = Event->get_particle(kaon_minus_index[ij])->get_y()   ; 
         double z1   = Event->get_particle(kaon_minus_index[ij])->get_z()   ; 

         double Px2  = Event->get_particle(pion_plus_index[ik])->get_px()  ; 
         double Py2  = Event->get_particle(pion_plus_index[ik])->get_py()  ; 
         double Pz2  = Event->get_particle(pion_plus_index[ik])->get_pz()  ; 
         double E2   = Event->get_particle(pion_plus_index[ik])->get_e()   ; 
         double t2   = Event->get_particle(pion_plus_index[ik])->get_t()   ; 
         double x2   = Event->get_particle(pion_plus_index[ik])->get_x()   ; 
         double y2   = Event->get_particle(pion_plus_index[ik])->get_y()   ; 
         double z2   = Event->get_particle(pion_plus_index[ik])->get_z()   ; 

         // compute the invariant mass
         double E = E1 + E2;
         double Px = Px1 + Px2;
         double Py = Py1 + Py2;
         double Pz = Pz1 + Pz2;
         double invariant_mass = sqrt(E*E - Px*Px - Py*Py - Pz*Pz);

         if (std::abs(invariant_mass - kstar0_bar_mass)// k*0(892)-bar resonance found
                                      < kstar0_bar_width) {
            double spatial_distance = sqrt(
                 (t1 - t2)*(t1 - t2) + (x1 - x2)*(x1 - x2)
                   + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)) ;
             if (spatial_distance < 0.001) {
                 Event->add_particle(-313, t1, x1, y1, z1, E, Px, Py, Pz,true,1) ;
             }
         } // if k*0(892)-bar found
   
       } // loop over kaon minus
     } // loop over pion plus
  } // event loop

}



