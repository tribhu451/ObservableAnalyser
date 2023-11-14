#include "observables.h"

observables::observables(input_paramters &iparam_, read_input_file* rif_) : iparam(iparam_), rif(rif_) {
}


void observables::calculate_dnchdeta_eta( double pT_min, double pT_max ){

  const int _N_HISTOGRAMS_DNCHDETA_ETA = 12 ; 

  double Eta_min   = -8. ; 
  double Eta_max   =  8. ; 
  int    Eta_bins  =  48 ; 

  TH1D*                 H1D_DNCHDETA_ETA[_N_HISTOGRAMS_DNCHDETA_ETA] ; 
  H1D_DNCHDETA_ETA[0] = new TH1D("H0", "charged_particle_eta_differential_yield", Eta_bins, Eta_min, Eta_max);
  H1D_DNCHDETA_ETA[1] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H1");  H1D_DNCHDETA_ETA[1]->SetTitle("pion_plus");
  H1D_DNCHDETA_ETA[2] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H2");  H1D_DNCHDETA_ETA[2]->SetTitle("pion_minus");
  H1D_DNCHDETA_ETA[3] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H3");  H1D_DNCHDETA_ETA[3]->SetTitle("kaon_plus");
  H1D_DNCHDETA_ETA[4] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H4");  H1D_DNCHDETA_ETA[4]->SetTitle("kaon_minus");
  H1D_DNCHDETA_ETA[5] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H5");  H1D_DNCHDETA_ETA[5]->SetTitle("proton");
  H1D_DNCHDETA_ETA[6] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H6");  H1D_DNCHDETA_ETA[6]->SetTitle("anti_proton");
  H1D_DNCHDETA_ETA[7] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H7");  H1D_DNCHDETA_ETA[7]->SetTitle("lambda");
  H1D_DNCHDETA_ETA[8] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H8");  H1D_DNCHDETA_ETA[8]->SetTitle("anti_lambda");
  H1D_DNCHDETA_ETA[9] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H9");  H1D_DNCHDETA_ETA[9]->SetTitle("kstar0");
  H1D_DNCHDETA_ETA[10] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H10");  H1D_DNCHDETA_ETA[10]->SetTitle("kstra0_bar");
  H1D_DNCHDETA_ETA[11] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H11");  H1D_DNCHDETA_ETA[11]->SetTitle("phi");

  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double W   = Event->get_particle(jj)->get_weight()   ; 
      double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double Eta = 0.5 * TMath::Log( ( P + Pz ) / (P - Pz) );
      double Pt = sqrt( Px * Px + Py * Py )  ;
      if(Pt > pT_max || Pt < pT_min)
        continue ; 

     if(PID ==  211  || PID == -211  || PID ==  321 || PID == -321  || 
        PID == 2212  || PID == -2212 || PID == 3122 || PID == -3122 ||
        PID == 3222  || PID == -3222 || PID == 3112 || PID == -3112 ||
        PID == 3312  || PID == -3312 || PID == 3334 || PID == -3334  ){
        H1D_DNCHDETA_ETA[0]->Fill(Eta,1.);
      }

      if(PID == 211){
        H1D_DNCHDETA_ETA[1]->Fill(Eta,W);
      }
      if(PID == -211){
        H1D_DNCHDETA_ETA[2]->Fill(Eta,W);
      }
      if(PID == 321){
        H1D_DNCHDETA_ETA[3]->Fill(Eta,W);
      }
      if(PID == -321){
        H1D_DNCHDETA_ETA[4]->Fill(Eta,W);
      }
      if(PID == 2212){
        H1D_DNCHDETA_ETA[5]->Fill(Eta,W);
      }
      if(PID == -2212){
        H1D_DNCHDETA_ETA[6]->Fill(Eta,W);
      }
      if(PID == 3122){
        H1D_DNCHDETA_ETA[7]->Fill(Eta,W);
      }
      if(PID == -3122){
        H1D_DNCHDETA_ETA[8]->Fill(Eta,W);
      }
      if(PID == 313){
        H1D_DNCHDETA_ETA[9]->Fill(Eta,W);
      }
      if(PID == -313){
        H1D_DNCHDETA_ETA[10]->Fill(Eta,W);
      }
      if(PID == 333){
        H1D_DNCHDETA_ETA[11]->Fill(Eta,W);
      }

     } // particle loop
    } // event loop


  double dX = ( Eta_max - Eta_min ) / Eta_bins ; 
  for(int ii=0; ii<_N_HISTOGRAMS_DNCHDETA_ETA; ii++){
        H1D_DNCHDETA_ETA[ii]->Scale(1.0 / (nEvents * dX));
  }

  std::ofstream mFile;
  std::stringstream output_filename;

  std::string hadron_name[_N_HISTOGRAMS_DNCHDETA_ETA] = {"hpm", "211", "-211", "321", "-321", "2212", "-2212", "3122", "-3122", "313", "-313", "333" };
  int        hadron_index[_N_HISTOGRAMS_DNCHDETA_ETA] = {  0,     1 ,     2 ,    3 ,     4 ,     5 ,      6 ,     7 ,      8  ,   9  ,    10 ,   11  };


 for(int ix =0; ix < _N_HISTOGRAMS_DNCHDETA_ETA; ix++){
   output_filename.str("");
   output_filename << "results/dnchdeta_eta-" << hadron_name[ix];
   output_filename << "_pT_" << pT_min << "_" << pT_max ;
   output_filename << ".dat";

   mFile.open(output_filename.str().c_str(), std::ios::out );
   for(int i=1; i<= H1D_DNCHDETA_ETA[hadron_index[ix]]->GetNbinsX(); i++){
       mFile << H1D_DNCHDETA_ETA[hadron_index[ix]]->GetBinCenter(i) << "\t" << H1D_DNCHDETA_ETA[hadron_index[ix]]->GetBinContent(i) 
              << "\t" << H1D_DNCHDETA_ETA[hadron_index[ix]]->GetBinError(i) << std::endl;
   }
   mFile.close();
 }

   for(int ixx=0; ixx < _N_HISTOGRAMS_DNCHDETA_ETA ; ixx++){
     H1D_DNCHDETA_ETA[ixx]->Clear(); 
   }

delete H1D_DNCHDETA_ETA[_N_HISTOGRAMS_DNCHDETA_ETA] ; 

}



void observables::calculate_dndy_y( double pT_min, double pT_max ){

  const int _N_HISTOGRAMS_DNDY_Y = 12 ; 

  double Y_min   = -8. ; 
  double Y_max   =  8. ; 
  int    Y_bins  =  48 ; 

  TH1D*                 H1D_DNDY_Y[_N_HISTOGRAMS_DNDY_Y] ; 
  H1D_DNDY_Y[0] = new TH1D("YH0", "charged_particle_eta_differential_yield", Y_bins, Y_min, Y_max);
  H1D_DNDY_Y[1] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH1");  H1D_DNDY_Y[1]->SetTitle("pion_plus");
  H1D_DNDY_Y[2] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH2");  H1D_DNDY_Y[2]->SetTitle("pion_minus");
  H1D_DNDY_Y[3] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH3");  H1D_DNDY_Y[3]->SetTitle("kaon_plus");
  H1D_DNDY_Y[4] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH4");  H1D_DNDY_Y[4]->SetTitle("kaon_minus");
  H1D_DNDY_Y[5] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH5");  H1D_DNDY_Y[5]->SetTitle("proton");
  H1D_DNDY_Y[6] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH6");  H1D_DNDY_Y[6]->SetTitle("anti_proton");
  H1D_DNDY_Y[7] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH7");  H1D_DNDY_Y[7]->SetTitle("lambda");
  H1D_DNDY_Y[8] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH8");  H1D_DNDY_Y[8]->SetTitle("anti_lambda");
  H1D_DNDY_Y[9] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH9");  H1D_DNDY_Y[9]->SetTitle("kstar0");
  H1D_DNDY_Y[10] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH10");  H1D_DNDY_Y[10]->SetTitle("kstar0_bar");
  H1D_DNDY_Y[11] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH11");  H1D_DNDY_Y[11]->SetTitle("phi");

  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ; 
      double W   = Event->get_particle(jj)->get_weight() ; 
      double Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
      double Pt = sqrt( Px * Px + Py * Py )  ;
      if(Pt > pT_max || Pt < pT_min)
        continue ; 

     if(PID ==  211  || PID == -211  || PID ==  321 || PID == -321  || 
        PID == 2212  || PID == -2212 || PID == 3122 || PID == -3122 ||
        PID == 3222  || PID == -3222 || PID == 3112 || PID == -3112 ||
        PID == 3312  || PID == -3312 || PID == 3334 || PID == -3334  ){
        H1D_DNDY_Y[0]->Fill(Rap,W);
      }

      if(PID == 211){
        H1D_DNDY_Y[1]->Fill(Rap,W);
      }
      if(PID == -211){
        H1D_DNDY_Y[2]->Fill(Rap,W);
      }
      if(PID == 321){
        H1D_DNDY_Y[3]->Fill(Rap,W);
      }
      if(PID == -321){
        H1D_DNDY_Y[4]->Fill(Rap,W);
      }
      if(PID == 2212){
        H1D_DNDY_Y[5]->Fill(Rap,W);
      }
      if(PID == -2212){
        H1D_DNDY_Y[6]->Fill(Rap,W);
      }
      if(PID == 3122){
        H1D_DNDY_Y[7]->Fill(Rap,W);
      }
      if(PID == -3122){
        H1D_DNDY_Y[8]->Fill(Rap,W);
      }
      if(PID == 313){
        H1D_DNDY_Y[9]->Fill(Rap,W);
      }
      if(PID == -313){
        H1D_DNDY_Y[10]->Fill(Rap,W);
      }
      if(PID == 333){
        H1D_DNDY_Y[11]->Fill(Rap,W);
      }
     } // particle loop
    } // event loop


  double dX = ( Y_max - Y_min ) / Y_bins ; 
  for(int ii=0; ii<_N_HISTOGRAMS_DNDY_Y; ii++)
        H1D_DNDY_Y[ii]->Scale(1.0 / (nEvents * dX));

  std::ofstream mFile;
  std::stringstream output_filename;

  std::string hadron_name[_N_HISTOGRAMS_DNDY_Y] = {"hpm", "211", "-211", "321", "-321", "2212", "-2212", "3122", "-3122", "313", "-313", "333" };
  int        hadron_index[_N_HISTOGRAMS_DNDY_Y] = {  0,     1 ,     2 ,    3 ,     4 ,     5 ,      6 ,     7 ,      8  ,   9  ,    10 ,   11  };


 for(int ix =0; ix < _N_HISTOGRAMS_DNDY_Y; ix++){
   output_filename.str("");
   output_filename << "results/dndy_y-" << hadron_name[ix];
   output_filename << "_pT_" << pT_min << "_" << pT_max ;
   output_filename << ".dat";

   mFile.open(output_filename.str().c_str(), std::ios::out );
   for(int i=1; i<= H1D_DNDY_Y[hadron_index[ix]]->GetNbinsX(); i++){
       mFile << H1D_DNDY_Y[hadron_index[ix]]->GetBinCenter(i) << "\t" << H1D_DNDY_Y[hadron_index[ix]]->GetBinContent(i) 
              << "\t" << H1D_DNDY_Y[hadron_index[ix]]->GetBinError(i) << std::endl;
   }
   mFile.close();
 }

   for(int ixx=0; ixx < _N_HISTOGRAMS_DNDY_Y ; ixx++){
     H1D_DNDY_Y[ixx]->Clear(); 
   }


 delete H1D_DNDY_Y[_N_HISTOGRAMS_DNDY_Y] ; 

}



void observables::calculate_invariant_yield_vs_pt( int yflag, double Rap_min, double Rap_max ){

  const int _N_HISTOGRAMS_INVYLD_PT_ = 11 ; 

  double pT_min  = 0.05 ; 
  double pT_max  = 3.0  ;
  int    pT_bins = 24   ;
  double dpT = ( pT_max - pT_min ) / pT_bins ; 
  double dY = ( Rap_max - Rap_min )  ; 

  TH1D*                 H1D_INVYLD_PT[_N_HISTOGRAMS_INVYLD_PT_] ; 
  H1D_INVYLD_PT[0] = new TH1D("PT0", "inavriant_yield_pion_plus", pT_bins, pT_min, pT_max);
  H1D_INVYLD_PT[1] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT1");  H1D_INVYLD_PT[1]->SetTitle("pion_minus");
  H1D_INVYLD_PT[2] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT2");  H1D_INVYLD_PT[2]->SetTitle("kaon_plus");
  H1D_INVYLD_PT[3] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT3");  H1D_INVYLD_PT[3]->SetTitle("kaon_minus");
  H1D_INVYLD_PT[4] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT4");  H1D_INVYLD_PT[4]->SetTitle("proton");
  H1D_INVYLD_PT[5] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT5");  H1D_INVYLD_PT[5]->SetTitle("anti_proton");
  H1D_INVYLD_PT[6] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT6");  H1D_INVYLD_PT[6]->SetTitle("Lambda");
  H1D_INVYLD_PT[7] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT7");  H1D_INVYLD_PT[7]->SetTitle("anti_Lambda");
  H1D_INVYLD_PT[8] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT8");  H1D_INVYLD_PT[8]->SetTitle("kstar0");
  H1D_INVYLD_PT[9] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT9");  H1D_INVYLD_PT[9]->SetTitle("kstar0_bar");
  H1D_INVYLD_PT[10] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT10");  H1D_INVYLD_PT[10]->SetTitle("phi");


  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ; 
      double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double W   = Event->get_particle(jj)->get_weight() ; 
      double Rap ; 

      if( fabs(E-Pz) < 1E-10 )
        continue ; 

      if( fabs(P-Pz) < 1E-10 )
        continue ; 


      if (yflag > 0 )
        Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
      else
        Rap = 0.5 * TMath::Log( ( P + Pz ) / (P - Pz) );

      double Pt = sqrt( Px * Px + Py * Py )  ;
      if(Rap > Rap_max || Rap < Rap_min)
        continue ; 

      if(PID == 211){
        H1D_INVYLD_PT[0]->Fill(Pt,W/Pt);
      }
      if(PID == -211){
        H1D_INVYLD_PT[1]->Fill(Pt,W/Pt);
      }
      if(PID == 321){
        H1D_INVYLD_PT[2]->Fill(Pt,W/Pt);
      }
      if(PID == -321){
        H1D_INVYLD_PT[3]->Fill(Pt,W/Pt);
      }
      if(PID == 2212){
        H1D_INVYLD_PT[4]->Fill(Pt,W/Pt);
      }
      if(PID == -2212){
        H1D_INVYLD_PT[5]->Fill(Pt,W/Pt);
      }
      if(PID == 3122){
        H1D_INVYLD_PT[6]->Fill(Pt,W/Pt);
      }
      if(PID == -3122){
        H1D_INVYLD_PT[7]->Fill(Pt,W/Pt);
      }
      if(PID == 313){
        H1D_INVYLD_PT[8]->Fill(Pt,W/Pt);
      }
      if(PID == -313){
        H1D_INVYLD_PT[9]->Fill(Pt,W/Pt);
      }
      if(PID == 333){
        H1D_INVYLD_PT[10]->Fill(Pt,W/Pt);
      }
     } // particle loop
    } // event loop

  for(int iix=0; iix<_N_HISTOGRAMS_INVYLD_PT_; iix++){
    H1D_INVYLD_PT[iix]->Scale( 1.0 / ( nEvents  * 2.0 * TMath::Pi() * dY * dpT ) );
  }

  std::ofstream mFile;
  std::stringstream output_filename;

  std::string hadron_name[_N_HISTOGRAMS_INVYLD_PT_] = { "211", "-211", "321", "-321", "2212", "-2212", "3122", "-3122", "313", "-313", "333" };
  int        hadron_index[_N_HISTOGRAMS_INVYLD_PT_] = {   0,     1 ,     2 ,    3 ,     4 ,     5 ,      6 ,     7 ,      8  ,   9  ,    10  };


 for(int ix =0; ix < _N_HISTOGRAMS_INVYLD_PT_; ix++){
   output_filename.str("");

   if(yflag > 0 ){
     output_filename << "results/dnptdptdy_pt-" << hadron_name[ix];
     output_filename << "_y_" << Rap_min << "_" << Rap_max ;
   }
  else{
     output_filename << "results/dnptdptdy_pt-" << hadron_name[ix];
     output_filename << "_eta_" << Rap_min << "_" << Rap_max ;
  }
   output_filename << ".dat";

   mFile.open(output_filename.str().c_str(), std::ios::out );
   for(int i=1; i<= H1D_INVYLD_PT[hadron_index[ix]]->GetNbinsX(); i++){
       mFile << H1D_INVYLD_PT[hadron_index[ix]]->GetBinCenter(i) << "\t" << H1D_INVYLD_PT[hadron_index[ix]]->GetBinContent(i) 
              << "\t" << H1D_INVYLD_PT[hadron_index[ix]]->GetBinError(i) << std::endl;
   }
   mFile.close();
 }

   for(int ixx=0; ixx < _N_HISTOGRAMS_INVYLD_PT_ ; ixx++){
     H1D_INVYLD_PT[ixx]->Clear(); 
   }


}


void observables::calculate_v1_vs_y_or_eta(int yflag, double psi1,  double pT_min, double pT_max ){

   double Y_min   = -5. ; 
   double Y_max   =  5. ; 
   int    Y_bins  =  40 ; 
 
   const int   _N_HISTOGRAMS_V1_Y_ = 16 ; 
   TProfile*   PROFILE_V1_Y[_N_HISTOGRAMS_V1_Y_] ; 
   PROFILE_V1_Y[0] = new TProfile("PV0", "v1_y_charged_particle", Y_bins, Y_min, Y_max);
   PROFILE_V1_Y[1] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV1");  PROFILE_V1_Y[1]->SetTitle("pion_plus");
   PROFILE_V1_Y[2] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV2");  PROFILE_V1_Y[2]->SetTitle("pion_minus");
   PROFILE_V1_Y[3] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV3");  PROFILE_V1_Y[3]->SetTitle("kaon_plus");
   PROFILE_V1_Y[4] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV4");  PROFILE_V1_Y[4]->SetTitle("kaon_minus");
   PROFILE_V1_Y[5] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV5");  PROFILE_V1_Y[5]->SetTitle("proton");
   PROFILE_V1_Y[6] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV6");  PROFILE_V1_Y[6]->SetTitle("anti_proton");
   PROFILE_V1_Y[7] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV7");  PROFILE_V1_Y[7]->SetTitle("lambda");
   PROFILE_V1_Y[8] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV8");  PROFILE_V1_Y[8]->SetTitle("anti_lambda");
   PROFILE_V1_Y[9] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV9");  PROFILE_V1_Y[9]->SetTitle("kstar0");
   PROFILE_V1_Y[10] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV10");  PROFILE_V1_Y[10]->SetTitle("kstar0_bar");
   PROFILE_V1_Y[11] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV11");  PROFILE_V1_Y[11]->SetTitle("phi");
   PROFILE_V1_Y[12] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV12");  PROFILE_V1_Y[12]->SetTitle("xi");
   PROFILE_V1_Y[13] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV13");  PROFILE_V1_Y[13]->SetTitle("xi_bar");
   PROFILE_V1_Y[14] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV14");  PROFILE_V1_Y[14]->SetTitle("omega");
   PROFILE_V1_Y[15] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV15");  PROFILE_V1_Y[15]->SetTitle("omega_bar");

   int nEvents = rif->get_event_buffer_size() ; 
   for(int ii=0; ii<nEvents; ii++){
     events* Event = rif->get_event(ii) ; 
     int nParticles = Event->get_multiplicity_of_the_event();
     for(int jj=0; jj<nParticles; jj++){
       int    PID = Event->get_particle(jj)->get_pid() ; 
       double Px  = Event->get_particle(jj)->get_px()  ; 
       double Py  = Event->get_particle(jj)->get_py()  ; 
       double Pz  = Event->get_particle(jj)->get_pz()  ; 
       double E   = Event->get_particle(jj)->get_e()   ; 
       double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
       double Pt = sqrt( Px * Px + Py * Py )  ;
       double W   = Event->get_particle(jj)->get_weight() ; 
       
       if(Pt > pT_max || Pt < pT_min)
	 continue ; 
       
       if( fabs(E-Pz) < 1E-10 )
	 continue ; 
       
       if( fabs(P-Pz) < 1E-10 )
	 continue ; 
       
       double Rap ;
              
       if (yflag > 0 )
	 Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
       else
	 Rap = 0.5 * TMath::Log( ( P + Pz ) / (P - Pz) );
       
       double v1 = Px / Pt ; 
       
       if(PID ==  211  || PID == -211  || PID ==  321 || PID == -321  || 
	  PID == 2212  || PID == -2212 || PID == 3122 || PID == -3122 ||
	  PID == 3222  || PID == -3222 || PID == 3112 || PID == -3112 ||
	  PID == 3312  || PID == -3312 || PID == 3334 || PID == -3334  ){
	 PROFILE_V1_Y[0]->Fill(Rap,v1,W);
       }
       
       if(PID == 211){
	 PROFILE_V1_Y[1]->Fill(Rap,v1,W);
       }
       if(PID == -211){
	 PROFILE_V1_Y[2]->Fill(Rap,v1,W);
       }
       if(PID == 321){
	 PROFILE_V1_Y[3]->Fill(Rap,v1,W);
       }
       if(PID == -321){
	 PROFILE_V1_Y[4]->Fill(Rap,v1,W);
       }
       if(PID == 2212){
	 PROFILE_V1_Y[5]->Fill(Rap,v1,W);
       }
       if(PID == -2212){
	 PROFILE_V1_Y[6]->Fill(Rap,v1,W);
       }
       if(PID == 3122){
	 PROFILE_V1_Y[7]->Fill(Rap,v1,W);
       }
       if(PID == -3122){
	 PROFILE_V1_Y[8]->Fill(Rap,v1,W);
       }
       if(PID == 313){
	 PROFILE_V1_Y[9]->Fill(Rap,v1,W);
       }
       if(PID == -313){
	 PROFILE_V1_Y[10]->Fill(Rap,v1,W);
       }
       if(PID == 333){
	 PROFILE_V1_Y[11]->Fill(Rap,v1,W);
       }
       if(PID == 3312){
	 PROFILE_V1_Y[12]->Fill(Rap,v1,W);
       }
       if(PID == -3312){
	 PROFILE_V1_Y[13]->Fill(Rap,v1,W);
       }
       if(PID == 3334){
	 PROFILE_V1_Y[14]->Fill(Rap,v1,W);
       }
       if(PID == -3334){
	 PROFILE_V1_Y[15]->Fill(Rap,v1,W);
       }
       
     } // particle loop
   } // event loop
   
   std::ofstream mFile;
   std::stringstream output_filename;
   
   std::string hadron_name[_N_HISTOGRAMS_V1_Y_] = {"hpm", "211", "-211", "321", "-321", "2212", "-2212", "3122", "-3122", "313", "-313", "333", "3312", "-3312", "3334", "-3334" };
   int  hadron_index[_N_HISTOGRAMS_V1_Y_] = {  0,   1 ,   2 ,  3 ,   4,   5,   6,  7,   8,   9,    10,   11,   12,    13,   14,   15  };
   
   
   for(int ix =0; ix < _N_HISTOGRAMS_V1_Y_; ix++){
     output_filename.str("");
     
     if(yflag > 0 ){
       output_filename << "results/v1_y-" << hadron_name[ix];
       output_filename << "_pT_" << pT_min << "_" << pT_max ;
     }
     else{
       output_filename << "results/v1_eta-" << hadron_name[ix];
       output_filename << "_pT_" << pT_min << "_" << pT_max ;
     }
     output_filename << ".dat";
     
     mFile.open(output_filename.str().c_str(), std::ios::out );
     for(int i=1; i<= PROFILE_V1_Y[hadron_index[ix]]->GetNbinsX(); i++){
       mFile << PROFILE_V1_Y[hadron_index[ix]]->GetBinCenter(i) << "\t" << PROFILE_V1_Y[hadron_index[ix]]->GetBinContent(i) 
	     << "\t" << PROFILE_V1_Y[hadron_index[ix]]->GetBinError(i) << std::endl;
     }
     mFile.close();
   }
   
   
   // Slope calculation //
   int fit_range_iterations = 6 ;
   double fit_range[fit_range_iterations] = { 0.55, 0.65, 0.75, 0.85, 0.95, 1.05 } ;
   
   for(int ifri=0; ifri<fit_range_iterations; ifri++){
     for(int iff=0; iff<2; iff++){
       
       for(int ix =0; ix < _N_HISTOGRAMS_V1_Y_; ix++){
	 output_filename.str("");
	 
	 if(yflag > 0 ){
	   output_filename << "results/v1_y-" << hadron_name[ix];
	   output_filename << "_pT_" << pT_min << "_" << pT_max ;
	 }
	 else{
	   output_filename << "results/v1_eta-" << hadron_name[ix];
	   output_filename << "_pT_" << pT_min << "_" << pT_max ;
	 }
	 
	 output_filename << "_slope_fitrange_";
	 output_filename << fit_range[ifri] ;
	 output_filename << "_fit_func";
	 
	 if(iff == 0 ){
	   output_filename << "_linear" ;
	 }
	 else{
	   output_filename << "_linear_plus_cubic" ;
	 }
	 
	 output_filename << ".dat";
	 
	 double xx_val[20] ;
	 double yy_val[20] ;
	 double yy_up_err[20] ;
	 double yy_low_err[20] ;
	 for(int id=0; id<20; id++){
	   xx_val[id] = 0. ;
	   yy_val[id] = 0. ;
	   yy_up_err[id] = 0. ;
	   yy_low_err[id] = 0. ;
	 }
	 int npoints = 0 ;
	 
	 for(int i=1; i<= PROFILE_V1_Y[hadron_index[ix]]->GetNbinsX(); i++){
	   if( fabs( PROFILE_V1_Y[hadron_index[ix]]->GetBinCenter(i) ) > fit_range[ifri] )
	     continue ;
	   xx_val[npoints] = PROFILE_V1_Y[hadron_index[ix]]->GetBinCenter(i)    ;
	   yy_val[npoints] = PROFILE_V1_Y[hadron_index[ix]]->GetBinContent(i)   ;
	   if( PROFILE_V1_Y[hadron_index[ix]]->GetBinContent(i) > 0 ){
	     yy_up_err[npoints] = PROFILE_V1_Y[hadron_index[ix]]->GetBinContent(i) +
	       PROFILE_V1_Y[hadron_index[ix]]->GetBinError(i)   ;
	     yy_low_err[npoints] = PROFILE_V1_Y[hadron_index[ix]]->GetBinContent(i) -
	       PROFILE_V1_Y[hadron_index[ix]]->GetBinError(i)   ;
	   }
	   else{
	     yy_up_err[npoints] = PROFILE_V1_Y[hadron_index[ix]]->GetBinContent(i) -
	       PROFILE_V1_Y[hadron_index[ix]]->GetBinError(i)   ;
	     yy_low_err[npoints] = PROFILE_V1_Y[hadron_index[ix]]->GetBinContent(i) +
	       PROFILE_V1_Y[hadron_index[ix]]->GetBinError(i)   ;
	   }
	   
	   npoints += 1 ;
	 }
	 
	 
	 double slope_val ;
	 double slope_up_err ;
	 double slope_low_err ;
	 
	 if(iff == 0 ){
	   slope_val = fit_a_straight_line_and_get_slope(npoints,xx_val,yy_val) ;
	   slope_up_err = fit_a_straight_line_and_get_slope(npoints,xx_val,yy_up_err) ;
	   slope_low_err = fit_a_straight_line_and_get_slope(npoints,xx_val,yy_low_err) ;
	 }
	 else{
	   slope_val = fit_a_cubic_plus_linear_func(npoints,xx_val,yy_val) ;
	   slope_up_err = fit_a_cubic_plus_linear_func(npoints,xx_val,yy_up_err) ;
	   slope_low_err = fit_a_cubic_plus_linear_func(npoints,xx_val,yy_low_err) ;
	 }
	 mFile.open(output_filename.str().c_str(), std::ios::out );
	 mFile << npoints << "  " << slope_val << "  " << (slope_up_err - slope_val) << "  " << ( slope_val - slope_low_err ) << endl ;
	 // Actual slope Error = Slope found from  points with ( value + error ) - (value)
	 mFile.close();
	 
       }
     } // iff 
   } // ifri
   // slope calculation end //
   
   
   
   for(int ixx=0; ixx < _N_HISTOGRAMS_V1_Y_ ; ixx++){
     PROFILE_V1_Y[ixx]->Clear(); 
   }
   
   
}





void observables::calculate_v2_vs_y_or_eta(int yflag, double psi2,  double pT_min, double pT_max ){

   double Y_min   = -8. ; 
   double Y_max   =  8. ; 
   int    Y_bins  =  81 ; 
 
   const int   _N_HISTOGRAMS_V2_Y_ = 16 ; 
   TProfile*   PROFILE_V2_Y[_N_HISTOGRAMS_V2_Y_] ; 
   PROFILE_V2_Y[0] = new TProfile("PV0", "v2_y_charged_particle", Y_bins, Y_min, Y_max);
   PROFILE_V2_Y[1] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV1");  PROFILE_V2_Y[1]->SetTitle("pion_plus");
   PROFILE_V2_Y[2] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV2");  PROFILE_V2_Y[2]->SetTitle("pion_minus");
   PROFILE_V2_Y[3] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV3");  PROFILE_V2_Y[3]->SetTitle("kaon_plus");
   PROFILE_V2_Y[4] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV4");  PROFILE_V2_Y[4]->SetTitle("kaon_minus");
   PROFILE_V2_Y[5] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV5");  PROFILE_V2_Y[5]->SetTitle("proton");
   PROFILE_V2_Y[6] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV6");  PROFILE_V2_Y[6]->SetTitle("anti_proton");
   PROFILE_V2_Y[7] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV7");  PROFILE_V2_Y[7]->SetTitle("lambda");
   PROFILE_V2_Y[8] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV8");  PROFILE_V2_Y[8]->SetTitle("anti_lambda");
   PROFILE_V2_Y[9] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV9");  PROFILE_V2_Y[9]->SetTitle("kstar0");
   PROFILE_V2_Y[10] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV10");  PROFILE_V2_Y[10]->SetTitle("kstar0_bar");
   PROFILE_V2_Y[11] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV11");  PROFILE_V2_Y[11]->SetTitle("phi");
   PROFILE_V2_Y[12] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV12");  PROFILE_V2_Y[12]->SetTitle("xi");
   PROFILE_V2_Y[13] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV13");  PROFILE_V2_Y[13]->SetTitle("xi_bar");
   PROFILE_V2_Y[14] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV14");  PROFILE_V2_Y[14]->SetTitle("omega");
   PROFILE_V2_Y[15] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV15");  PROFILE_V2_Y[15]->SetTitle("omega_bar");

  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ; 
      double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double Pt = sqrt( Px * Px + Py * Py )  ;
      double W   = Event->get_particle(jj)->get_weight() ; 

      if(Pt > pT_max || Pt < pT_min)
        continue ; 

      if( fabs(E-Pz) < 1E-10 )
        continue ; 

      if( fabs(P-Pz) < 1E-10 )
        continue ; 

      double Rap ;


      if (yflag > 0 )
        Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
      else
        Rap = 0.5 * TMath::Log( ( P + Pz ) / (P - Pz) );

     double v2 = ( Px * Px - Py * Py ) / ( Pt * Pt ) ; 

     if(PID ==  211  || PID == -211  || PID ==  321 || PID == -321  || 
        PID == 2212  || PID == -2212 || PID == 3122 || PID == -3122 ||
        PID == 3222  || PID == -3222 || PID == 3112 || PID == -3112 ||
        PID == 3312  || PID == -3312 || PID == 3334 || PID == -3334  ){
        PROFILE_V2_Y[0]->Fill(Rap,v2,W);
      }

      if(PID == 211){
        PROFILE_V2_Y[1]->Fill(Rap,v2,W);
      }
      if(PID == -211){
        PROFILE_V2_Y[2]->Fill(Rap,v2,W);
      }
      if(PID == 321){
        PROFILE_V2_Y[3]->Fill(Rap,v2,W);
      }
      if(PID == -321){
        PROFILE_V2_Y[4]->Fill(Rap,v2,W);
      }
      if(PID == 2212){
        PROFILE_V2_Y[5]->Fill(Rap,v2,W);
      }
      if(PID == -2212){
        PROFILE_V2_Y[6]->Fill(Rap,v2,W);
      }
      if(PID == 3122){
        PROFILE_V2_Y[7]->Fill(Rap,v2,W);
      }
      if(PID == -3122){
        PROFILE_V2_Y[8]->Fill(Rap,v2,W);
      }
      if(PID == 313){
        PROFILE_V2_Y[9]->Fill(Rap,v2,W);
      }
      if(PID == -313){
        PROFILE_V2_Y[10]->Fill(Rap,v2,W);
      }
      if(PID == 333){
        PROFILE_V2_Y[11]->Fill(Rap,v2,W);
      }
      if(PID == 3312){
        PROFILE_V2_Y[12]->Fill(Rap,v2,W);
      }
      if(PID == -3312){
        PROFILE_V2_Y[13]->Fill(Rap,v2,W);
      }
      if(PID == 3334){
        PROFILE_V2_Y[14]->Fill(Rap,v2,W);
      }
      if(PID == -3334){
        PROFILE_V2_Y[15]->Fill(Rap,v2,W);
      }
     } // particle loop
    } // event loop

  std::ofstream mFile;
  std::stringstream output_filename;

  std::string hadron_name[_N_HISTOGRAMS_V2_Y_] = {"hpm", "211", "-211", "321", "-321", "2212", "-2212", "3122", "-3122", "313", "-313", "333", "3312", "-3312", "3334", "-3334" };
  int  hadron_index[_N_HISTOGRAMS_V2_Y_] = {  0,   1 ,   2 ,  3 ,   4,   5,   6,  7,   8,   9,    10,   11,   12,    13,   14,   15  };


 for(int ix =0; ix < _N_HISTOGRAMS_V2_Y_; ix++){
   output_filename.str("");

   if(yflag > 0 ){
     output_filename << "results/v2_y-" << hadron_name[ix];
     output_filename << "_pT_" << pT_min << "_" << pT_max ;
   }
  else{
     output_filename << "results/v2_eta-" << hadron_name[ix];
     output_filename << "_pT_" << pT_min << "_" << pT_max ;
  }
   output_filename << ".dat";

   mFile.open(output_filename.str().c_str(), std::ios::out );
   for(int i=1; i<= PROFILE_V2_Y[hadron_index[ix]]->GetNbinsX(); i++){
       mFile << PROFILE_V2_Y[hadron_index[ix]]->GetBinCenter(i) << "\t" << PROFILE_V2_Y[hadron_index[ix]]->GetBinContent(i) 
              << "\t" << PROFILE_V2_Y[hadron_index[ix]]->GetBinError(i) << std::endl;
   }
   mFile.close();
 }

   for(int ixx=0; ixx < _N_HISTOGRAMS_V2_Y_ ; ixx++){
      PROFILE_V2_Y[ixx]->Clear(); 
   }

}



void observables::calculate_v2_pt( int yflag, double Rap_min, double Rap_max ){

  const int _N_HISTOGRAMS_V2_PT = 16 ; 
  double v2_pt_bins[13] = {0.01, 0.1, 0.15, 0.3, 0.5, 0.75,1.0, 1.25, 1.5 , 1.75, 2.1, 2.5, 3.0} ; 

  TProfile*                 PROFILE_V2_PT[_N_HISTOGRAMS_V2_PT] ; 
  PROFILE_V2_PT[0]  = new TProfile("V2PT0", "pion_pT_differential_v2", 12, v2_pt_bins);
  PROFILE_V2_PT[1]  = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT00");  PROFILE_V2_PT[1]->SetTitle("pion_plus");
  PROFILE_V2_PT[2]  = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT1");   PROFILE_V2_PT[2]->SetTitle("pion_minus");
  PROFILE_V2_PT[3]  = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT2");   PROFILE_V2_PT[3]->SetTitle("kaon_plus");
  PROFILE_V2_PT[4]  = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT3");   PROFILE_V2_PT[4]->SetTitle("kaon_minus");
  PROFILE_V2_PT[5]  = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT4");   PROFILE_V2_PT[5]->SetTitle("proton");
  PROFILE_V2_PT[6]  = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT5");   PROFILE_V2_PT[6]->SetTitle("anti_proton");
  PROFILE_V2_PT[7]  = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT6");   PROFILE_V2_PT[7]->SetTitle("Lambda");
  PROFILE_V2_PT[8]  = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT7");   PROFILE_V2_PT[8]->SetTitle("anti_Lambda");
  PROFILE_V2_PT[9]  = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT8");   PROFILE_V2_PT[9]->SetTitle("kstar0");
  PROFILE_V2_PT[10] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT9");   PROFILE_V2_PT[10]->SetTitle("kstar0_bar");
  PROFILE_V2_PT[11] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT10");  PROFILE_V2_PT[11]->SetTitle("phi");
  PROFILE_V2_PT[12] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT11");  PROFILE_V2_PT[12]->SetTitle("xi");
  PROFILE_V2_PT[13] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT12");  PROFILE_V2_PT[13]->SetTitle("xi_bar");
  PROFILE_V2_PT[14] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT13");  PROFILE_V2_PT[14]->SetTitle("omega");
  PROFILE_V2_PT[15] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT14");  PROFILE_V2_PT[15]->SetTitle("omega_bar");

  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ; 
      double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double W   = Event->get_particle(jj)->get_weight() ; 
      double Rap ; 

      if( fabs(E-Pz) < 1E-10 )
        continue ; 

      if( fabs(P-Pz) < 1E-10 )
        continue ; 


      if (yflag > 0 )
        Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
      else
        Rap = 0.5 * TMath::Log( ( P + Pz ) / (P - Pz) );

      double Pt = sqrt( Px * Px + Py * Py )  ;
      if(Rap > Rap_max || Rap < Rap_min)
        continue ; 
      double v2 = ( Px * Px - Py * Py ) / ( Pt * Pt ) ;

     if(PID ==  211  || PID == -211  || PID ==  321 || PID == -321  || 
        PID == 2212  || PID == -2212 || PID == 3122 || PID == -3122 ||
        PID == 3222  || PID == -3222 || PID == 3112 || PID == -3112 ||
        PID == 3312  || PID == -3312 || PID == 3334 || PID == -3334  ){
        PROFILE_V2_PT[0]->Fill(Pt,v2,W);
      }

      if(PID == 211){
        PROFILE_V2_PT[1]->Fill(Pt,v2,W);
      }
      if(PID == -211){
        PROFILE_V2_PT[2]->Fill(Pt,v2,W);
      }
      if(PID == 321){
        PROFILE_V2_PT[3]->Fill(Pt,v2,W);
      }
      if(PID == -321){
        PROFILE_V2_PT[4]->Fill(Pt,v2,W);
      }
      if(PID == 2212){
        PROFILE_V2_PT[5]->Fill(Pt,v2,W);
      }
      if(PID == -2212){
        PROFILE_V2_PT[6]->Fill(Pt,v2,W);
      }
      if(PID == 3122){
        PROFILE_V2_PT[7]->Fill(Pt,v2,W);
      }
      if(PID == -3122){
        PROFILE_V2_PT[8]->Fill(Pt,v2,W);
      }
      if(PID == 313){
        PROFILE_V2_PT[9]->Fill(Pt,v2,W);
      }
      if(PID == -313){
        PROFILE_V2_PT[10]->Fill(Pt,v2,W);
      }
      if(PID == 333){
        PROFILE_V2_PT[11]->Fill(Pt,v2,W);
      }
      if(PID == 3312){
        PROFILE_V2_PT[12]->Fill(Pt,v2,W);
      }
      if(PID == -3312){
        PROFILE_V2_PT[13]->Fill(Pt,v2,W);
      }
      if(PID == 3334){
        PROFILE_V2_PT[14]->Fill(Pt,v2,W);
      }
      if(PID == -3334){
        PROFILE_V2_PT[15]->Fill(Pt,v2,W);
      }

     } // particle loop
    } // event loop


  std::ofstream mFile;
  std::stringstream output_filename;

  std::string hadron_name[_N_HISTOGRAMS_V2_PT] = {"hpm", "211", "-211", "321", "-321", "2212", "-2212", "3122", "-3122", "313", "-313", "333", "3312", "-3312", "3334", "-3334" };
  int        hadron_index[_N_HISTOGRAMS_V2_PT] = {  0,     1 ,     2 ,    3 ,     4 ,     5 ,      6 ,     7 ,      8  ,   9  ,    10 ,   11, 12,  13,  14,  15  };

 for(int ix =0; ix < _N_HISTOGRAMS_V2_PT; ix++){
   output_filename.str("");

   if(yflag > 0 ){
     output_filename << "results/v2_pt-" << hadron_name[ix];
     output_filename << "_y_" << Rap_min << "_" << Rap_max ;
   }
  else{
     output_filename << "results/v2_pt-" << hadron_name[ix];
     output_filename << "_eta_" << Rap_min << "_" << Rap_max ;
  }
   output_filename << ".dat";

   mFile.open(output_filename.str().c_str(), std::ios::out );
   for(int i=1; i<= PROFILE_V2_PT[hadron_index[ix]]->GetNbinsX(); i++){
       mFile << PROFILE_V2_PT[hadron_index[ix]]->GetBinCenter(i) << "\t" << PROFILE_V2_PT[hadron_index[ix]]->GetBinContent(i) 
              << "\t" << PROFILE_V2_PT[hadron_index[ix]]->GetBinError(i) << std::endl;
   }
   mFile.close();
 }

   for(int ixx=0; ixx < _N_HISTOGRAMS_V2_PT ; ixx++){
     PROFILE_V2_PT[ixx]->Clear(); 
   }

}



void observables::calculate_v1_pt( int yflag, double Rap_min, double Rap_max, int reflection_flag ){

  const int _N_HISTOGRAMS_V1_PT = 16 ; 
  double v1_pt_bins[13] = {0.01, 0.1, 0.15, 0.3, 0.5, 0.75,1.0, 1.25, 1.5 , 1.75, 2.1, 2.5, 3.0} ; 

  TProfile*                 PROFILE_V1_PT[_N_HISTOGRAMS_V1_PT] ; 
  PROFILE_V1_PT[0]  = new TProfile("V2PT0", "pion_pT_differential_v2", 12, v1_pt_bins);
  PROFILE_V1_PT[1]  = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT00");  PROFILE_V1_PT[1]->SetTitle("pion_plus");
  PROFILE_V1_PT[2]  = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT1");   PROFILE_V1_PT[2]->SetTitle("pion_minus");
  PROFILE_V1_PT[3]  = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT2");   PROFILE_V1_PT[3]->SetTitle("kaon_plus");
  PROFILE_V1_PT[4]  = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT3");   PROFILE_V1_PT[4]->SetTitle("kaon_minus");
  PROFILE_V1_PT[5]  = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT4");   PROFILE_V1_PT[5]->SetTitle("proton");
  PROFILE_V1_PT[6]  = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT5");   PROFILE_V1_PT[6]->SetTitle("anti_proton");
  PROFILE_V1_PT[7]  = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT6");   PROFILE_V1_PT[7]->SetTitle("Lambda");
  PROFILE_V1_PT[8]  = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT7");   PROFILE_V1_PT[8]->SetTitle("anti_Lambda");
  PROFILE_V1_PT[9]  = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT8");   PROFILE_V1_PT[9]->SetTitle("kstar0");
  PROFILE_V1_PT[10] = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT9");   PROFILE_V1_PT[10]->SetTitle("kstar0_bar");
  PROFILE_V1_PT[11] = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT10");  PROFILE_V1_PT[11]->SetTitle("phi");
  PROFILE_V1_PT[12] = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT11");  PROFILE_V1_PT[12]->SetTitle("xi");
  PROFILE_V1_PT[13] = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT12");  PROFILE_V1_PT[13]->SetTitle("xi_bar");
  PROFILE_V1_PT[14] = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT13");  PROFILE_V1_PT[14]->SetTitle("omega");
  PROFILE_V1_PT[15] = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT14");  PROFILE_V1_PT[15]->SetTitle("omega_bar");

  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ; 
      double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double W   = Event->get_particle(jj)->get_weight() ; 
      double Rap ; 

      if( fabs(E-Pz) < 1E-10 )
        continue ; 

      if( fabs(P-Pz) < 1E-10 )
        continue ; 


      if (yflag > 0 )
        Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
      else
        Rap = 0.5 * TMath::Log( ( P + Pz ) / (P - Pz) );

      double Pt = sqrt( Px * Px + Py * Py )  ;
      double v1 = Px / Pt ; ;

      if(Rap > Rap_max || Rap < Rap_min)
        continue ; 

      if(reflection_flag>0){
        if( fabs(Rap) > Rap_max || fabs(Rap) < Rap_min ){
          continue ;
        } 
        if(Rap < 0){
          v1 *= -1. ;
        }
      }
      else{
        if( Rap > Rap_max || Rap < Rap_min ){
          continue ;
        } 
      }

     if(PID ==  211  || PID == -211  || PID ==  321 || PID == -321  || 
        PID == 2212  || PID == -2212 || PID == 3122 || PID == -3122 ||
        PID == 3222  || PID == -3222 || PID == 3112 || PID == -3112 ||
        PID == 3312  || PID == -3312 || PID == 3334 || PID == -3334  ){
        PROFILE_V1_PT[0]->Fill(Pt,v1,W);
      }

      if(PID == 211){
        PROFILE_V1_PT[1]->Fill(Pt,v1,W);
      }
      if(PID == -211){
        PROFILE_V1_PT[2]->Fill(Pt,v1,W);
      }
      if(PID == 321){
        PROFILE_V1_PT[3]->Fill(Pt,v1,W);
      }
      if(PID == -321){
        PROFILE_V1_PT[4]->Fill(Pt,v1,W);
      }
      if(PID == 2212){
        PROFILE_V1_PT[5]->Fill(Pt,v1,W);
      }
      if(PID == -2212){
        PROFILE_V1_PT[6]->Fill(Pt,v1,W);
      }
      if(PID == 3122){
        PROFILE_V1_PT[7]->Fill(Pt,v1,W);
      }
      if(PID == -3122){
        PROFILE_V1_PT[8]->Fill(Pt,v1,W);
      }
      if(PID == 313){
        PROFILE_V1_PT[9]->Fill(Pt,v1,W);
      }
      if(PID == -313){
        PROFILE_V1_PT[10]->Fill(Pt,v1,W);
      }
      if(PID == 333){
        PROFILE_V1_PT[11]->Fill(Pt,v1,W);
      }
      if(PID == 3312){
        PROFILE_V1_PT[12]->Fill(Pt,v1,W);
      }
      if(PID == -3312){
        PROFILE_V1_PT[13]->Fill(Pt,v1,W);
      }
      if(PID == 3334){
        PROFILE_V1_PT[14]->Fill(Pt,v1,W);
      }
      if(PID == -3334){
        PROFILE_V1_PT[15]->Fill(Pt,v1,W);
      }


     } // particle loop
    } // event loop


  std::ofstream mFile;
  std::stringstream output_filename;

  std::string hadron_name[_N_HISTOGRAMS_V1_PT] = {"hpm", "211", "-211", "321", "-321", "2212", "-2212", "3122", "-3122", "313", "-313", "333", "3312", "-3312", "3334", "-3334" };
  int        hadron_index[_N_HISTOGRAMS_V1_PT] = {  0,     1 ,     2 ,    3 ,     4 ,     5 ,      6 ,     7 ,      8  ,   9  ,    10 ,   11, 12,  13,  14,  15  };


 for(int ix =0; ix < _N_HISTOGRAMS_V1_PT; ix++){
   output_filename.str("");

   if(yflag > 0 ){
     output_filename << "results/v1_pt-" << hadron_name[ix];
     output_filename << "_y_" << Rap_min << "_" << Rap_max ;
   }
  else{
     output_filename << "results/v1_pt-" << hadron_name[ix];
     output_filename << "_eta_" << Rap_min << "_" << Rap_max ;
  }

  if(reflection_flag>0){
   output_filename << "_reflected" ; 
  }
  else{
   output_filename << "_non_reflected" ; 
  }
   output_filename << ".dat";

   mFile.open(output_filename.str().c_str(), std::ios::out );
   for(int i=1; i<= PROFILE_V1_PT[hadron_index[ix]]->GetNbinsX(); i++){
       mFile << PROFILE_V1_PT[hadron_index[ix]]->GetBinCenter(i) << "\t" << PROFILE_V1_PT[hadron_index[ix]]->GetBinContent(i) 
              << "\t" << PROFILE_V1_PT[hadron_index[ix]]->GetBinError(i) << std::endl;
   }
   mFile.close();
 }

   for(int ixx=0; ixx < _N_HISTOGRAMS_V1_PT ; ixx++){
     PROFILE_V1_PT[ixx]->Clear(); 
   }

}


double observables::fit_a_straight_line_and_get_slope(int n, double *x, double *y){
  double sum_x     = 0. ; 
  double sum_x_sqr = 0. ;
  double sum_y     = 0. ; 
  double sum_xy    = 0. ; 

  for(int i=0; i<n; i++){
    sum_x += x[i] ;
    sum_x_sqr += ( x[i] * x[i] ) ; 
    sum_y += y[i] ;
    sum_xy += ( x[i] * y[i] ) ;
  }

  return ( n * sum_xy - sum_x * sum_y ) / ( n * sum_x_sqr - sum_x * sum_x ) ;
}


double observables::fit_a_cubic_plus_linear_func(int n, double *x, double *y){
  // fit function ax+bx^3 ; returns a
  double sum_x_raised_3_y_raised_1     = 0. ;
  double sum_x_raised_4                = 0. ;
  double sum_x_raised_1_y_raised_1     = 0. ;
  double sum_x_raised_6                = 0. ;
  double sum_x_raised_2                = 0. ;

  for(int i=0; i<n; i++){
    sum_x_raised_3_y_raised_1  += ( pow(x[i],3)*pow(y[i],1) ) ;
    sum_x_raised_4             += ( pow(x[i],4) ) ;
    sum_x_raised_1_y_raised_1  += ( pow(x[i],1)*pow(y[i],1) ) ;
    sum_x_raised_6             += ( pow(x[i],6) ) ;
    sum_x_raised_2             += ( pow(x[i],2) ) ;
  }

  return ( sum_x_raised_3_y_raised_1 * sum_x_raised_4 - sum_x_raised_1_y_raised_1 * sum_x_raised_6 ) /
              ( sum_x_raised_4 * sum_x_raised_4 - sum_x_raised_2 * sum_x_raised_6 ) ;
}



void observables::calculate_mean_pt_rap( int yflag, double pTmin, double pTmax){

   double Y_min   = -5. ; 
   double Y_max   =  5. ; 
   int    Y_bins  =  50 ; 
 
   const int   _N_HISTOGRAMS_MEANPT_Y_ = 12 ; 
   TProfile*   PROFILE_MEANPT_Y[_N_HISTOGRAMS_MEANPT_Y_] ; 
   PROFILE_MEANPT_Y[0]  = new TProfile("PV0", "meanpt_y_charged_particle", Y_bins, Y_min, Y_max);
   PROFILE_MEANPT_Y[1]  = (TProfile*) PROFILE_MEANPT_Y[0]->Clone("MPT1");   PROFILE_MEANPT_Y[1]->SetTitle("pion_plus");
   PROFILE_MEANPT_Y[2]  = (TProfile*) PROFILE_MEANPT_Y[0]->Clone("MPT2");   PROFILE_MEANPT_Y[2]->SetTitle("pion_minus");
   PROFILE_MEANPT_Y[3]  = (TProfile*) PROFILE_MEANPT_Y[0]->Clone("MPT3");   PROFILE_MEANPT_Y[3]->SetTitle("kaon_plus");
   PROFILE_MEANPT_Y[4]  = (TProfile*) PROFILE_MEANPT_Y[0]->Clone("MPT4");   PROFILE_MEANPT_Y[4]->SetTitle("kaon_minus");
   PROFILE_MEANPT_Y[5]  = (TProfile*) PROFILE_MEANPT_Y[0]->Clone("MPT5");   PROFILE_MEANPT_Y[5]->SetTitle("proton");
   PROFILE_MEANPT_Y[6]  = (TProfile*) PROFILE_MEANPT_Y[0]->Clone("MPT6");   PROFILE_MEANPT_Y[6]->SetTitle("anti_proton");
   PROFILE_MEANPT_Y[7]  = (TProfile*) PROFILE_MEANPT_Y[0]->Clone("MPT7");   PROFILE_MEANPT_Y[7]->SetTitle("lambda");
   PROFILE_MEANPT_Y[8]  = (TProfile*) PROFILE_MEANPT_Y[0]->Clone("MPT8");   PROFILE_MEANPT_Y[8]->SetTitle("anti_lambda");
   PROFILE_MEANPT_Y[9]  = (TProfile*) PROFILE_MEANPT_Y[0]->Clone("MPT9");   PROFILE_MEANPT_Y[9]->SetTitle("kstar0");
   PROFILE_MEANPT_Y[10] = (TProfile*) PROFILE_MEANPT_Y[0]->Clone("MPT10");  PROFILE_MEANPT_Y[10]->SetTitle("kstar0_bar");
   PROFILE_MEANPT_Y[11] = (TProfile*) PROFILE_MEANPT_Y[0]->Clone("MPT11");  PROFILE_MEANPT_Y[11]->SetTitle("phi");

  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ; 
      double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double W   = Event->get_particle(jj)->get_weight() ; 
      double Rap ; 

      if( fabs(E-Pz) < 1E-10 )
        continue ; 

      if( fabs(P-Pz) < 1E-10 )
        continue ; 


      if (yflag > 0 )
        Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
      else
        Rap = 0.5 * TMath::Log( ( P + Pz ) / (P - Pz) );

      double Pt = sqrt( Px * Px + Py * Py ) ;
      if(Pt > pTmax || Pt < pTmin)
        continue ; 



     if(PID ==  211  || PID == -211  || PID ==  321 || PID == -321  || 
        PID == 2212  || PID == -2212 || PID == 3122 || PID == -3122 ||
        PID == 3222  || PID == -3222 || PID == 3112 || PID == -3112 ||
        PID == 3312  || PID == -3312 || PID == 3334 || PID == -3334  ){
        PROFILE_MEANPT_Y[0]->Fill(Rap,Pt,W);
      }

      if(PID == 211){
        PROFILE_MEANPT_Y[1]->Fill(Rap,Pt,W);
      }
      if(PID == -211){
        PROFILE_MEANPT_Y[2]->Fill(Rap,Pt,W);
      }
      if(PID == 321){
        PROFILE_MEANPT_Y[3]->Fill(Rap,Pt,W);
      }
      if(PID == -321){
        PROFILE_MEANPT_Y[4]->Fill(Rap,Pt,W);
      }
      if(PID == 2212){
        PROFILE_MEANPT_Y[5]->Fill(Rap,Pt,W);
      }
      if(PID == -2212){
        PROFILE_MEANPT_Y[6]->Fill(Rap,Pt,W);
      }
      if(PID == 3122){
        PROFILE_MEANPT_Y[7]->Fill(Rap,Pt,W);
      }
      if(PID == -3122){
        PROFILE_MEANPT_Y[8]->Fill(Rap,Pt,W);
      }
      if(PID == 313){
        PROFILE_MEANPT_Y[9]->Fill(Rap,Pt,W);
      }
      if(PID == -313){
        PROFILE_MEANPT_Y[10]->Fill(Rap,Pt,W);
      }
      if(PID == 333){
        PROFILE_MEANPT_Y[11]->Fill(Rap,Pt,W);
      }

     } // particle loop
    } // event loop


  std::ofstream mFile;
  std::stringstream output_filename;

  std::string hadron_name[_N_HISTOGRAMS_MEANPT_Y_] = {"hpm", "211", "-211", "321", "-321", "2212", "-2212", "3122", "-3122", "313", "-313", "333" };
  int        hadron_index[_N_HISTOGRAMS_MEANPT_Y_] = {  0,     1 ,     2 ,    3 ,     4 ,     5 ,      6 ,     7 ,      8  ,   9  ,    10 ,   11  };

 for(int ix =0; ix < _N_HISTOGRAMS_MEANPT_Y_; ix++){
   output_filename.str("");

   if(yflag > 0 ){
     output_filename << "results/mean_pt_y-" << hadron_name[ix];
     output_filename << "_pT_" << pTmin << "_" << pTmax ;
   }
  else{
     output_filename << "results/mean_pt_eta-" << hadron_name[ix];
     output_filename << "_pT_" << pTmin << "_" << pTmax ;
  }
   output_filename << ".dat";

   mFile.open(output_filename.str().c_str(), std::ios::out );
   for(int i=1; i<= PROFILE_MEANPT_Y[hadron_index[ix]]->GetNbinsX(); i++){
       mFile << PROFILE_MEANPT_Y[hadron_index[ix]]->GetBinCenter(i) << "\t" << PROFILE_MEANPT_Y[hadron_index[ix]]->GetBinContent(i) 
              << "\t" << PROFILE_MEANPT_Y[hadron_index[ix]]->GetBinError(i) << std::endl;
   }
   mFile.close();
 }

   for(int ixx=0; ixx < _N_HISTOGRAMS_MEANPT_Y_ ; ixx++){
     PROFILE_MEANPT_Y[ixx]->Clear(); 
   }

}



void observables::calculate_amn( int ppid, int yflag, double ymin, double ymax, double ptcutmin, double ptcutmax){
  
  std::cout << "Calculating Amn of "  << ppid << " ... " <<  std::endl ; 
  double pT0 = diskptradius ; // in GeV (disc radius)
  int mmax=8 ; 
  int nmax=8 ; 
  double amnRe[mmax][nmax];
  double amnIm[mmax][nmax];
  double temp_amnRe[mmax][nmax];
  double temp_amnIm[mmax][nmax];
  double lambda[mmax][nmax];
  for(int im=0; im<mmax; im++){
    for(int in=0; in<nmax; in++){
      amnRe[im][in] = 0. ; 
      amnIm[im][in] = 0. ; 
      temp_amnRe[im][in] = 0. ; 
      temp_amnIm[im][in] = 0. ; 
      lambda[im][in] = 0. ; 
    }
  }
  double particles_norm = 0 ; 
  for(int im=0; im<mmax; im++){
    for(int in=0; in<nmax; in++){
      if(in==0){
        continue ;
      }
      lambda[im][in] = gsl_sf_bessel_zero_Jnu(fabs(im),in) ; 
    }
  }

  std::ofstream mFile;
  std::stringstream output_filename;
  output_filename.str("");
  if(yflag > 0 ){
    output_filename << "results/Amn-" << ppid ;
    output_filename << "_y_" << ymin << "_" << ymax ;
  }
  else{
    output_filename << "results/Amn-" << ppid;
    output_filename << "_eta_" << ymin << "_" << ymax ;
  }
  output_filename << "_pt_" << ptcutmin << "_" << ptcutmax ;
  output_filename << "_pt0_" << pT0 ;
  output_filename << ".dat";
  mFile.open(output_filename.str().c_str(), std::ios::out );
  mFile << "# m   n   AmnRe  AmnIm  |Amn| " << std::endl ; 
  
  
  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    if(ii%500 == 0){
      std::cout << "Analysing event " << ii << std::endl ; 
    }

    // set tempAmn for each event to zero at starting of event //
    for(int im=0; im<mmax; im++){
      for(int in=0; in<nmax; in++){
        temp_amnRe[im][in] = 0. ; 
        temp_amnIm[im][in] = 0. ; 
      }
    }

    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    double temp_particles_norm = 0 ; 
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ; 
      double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double W   = Event->get_particle(jj)->get_weight() ; 
      double Rap ; 
      
      if( fabs(PID-ppid) > 0.1 ){
	continue ; 
      }
      

      if( fabs(E-Pz) < 1E-10 )
	continue ; 
      
      if( fabs(P-Pz) < 1E-10 )
	continue ; 
      
      
      if (yflag > 0 )
	Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
      else
	Rap = 0.5 * TMath::Log( ( P + Pz ) / (P - Pz) );
      
      if(Rap > ymax || Rap < ymin ){
	continue ; 
      }
      
      double Pt = sqrt( Px * Px + Py * Py ) ;
      double Phi = atan2(Py,Px) ;
      if(Pt < ptcutmin || Pt > ptcutmax){
	continue ; 
      }
      
      temp_particles_norm += 1. ; 
      double lambda_mn ; 
      for(int im=0; im<mmax; im++){
	for(int in=0; in<nmax; in++){
	  
	  if(in==0){
	    continue ; 
	  }
	  
	  // get n-th zero of J_m(x) // 
	  lambda_mn =  lambda[im][in] ; 

	  temp_amnRe[im][in] +=  ( 1./gsl_sf_bessel_Jn( fabs(im)+1, lambda_mn ) * gsl_sf_bessel_Jn(im,  Pt / pT0 * lambda_mn ) * cos(im*Phi) )  ; 
	  temp_amnIm[im][in] +=  ( 1./gsl_sf_bessel_Jn( fabs(im)+1, lambda_mn ) * gsl_sf_bessel_Jn(im,  Pt / pT0 * lambda_mn ) * sin(im*Phi) * (-1.0) )  ; 
	} // in loop
      } // im loop
      
    } // particle loop

  for(int im=0; im<mmax; im++){
    for(int in=0; in<nmax; in++){
      amnRe[im][in] += temp_amnRe[im][in] ; 
      amnIm[im][in] += temp_amnIm[im][in] ; 
    }
  }
   particles_norm += temp_particles_norm ; 

  } // event loop
  
  
  for(int im=0; im<mmax; im++){
    for(int in=0; in<nmax; in++){
      amnRe[im][in] /= ( particles_norm * M_PI * pT0 * pT0 ) ; 
      amnIm[im][in] /= ( particles_norm * M_PI * pT0 * pT0 ) ; 
      mFile << im << "  " << in << "  " << amnRe[im][in] << "  " << amnIm[im][in] 
	    << "  " << sqrt(pow(amnRe[im][in],2)+pow(amnIm[im][in],2)) << std::endl ; 
    }
  }
  
  mFile.close();
}



void observables::calculate_amn_of_charged_hadrons( int yflag, double ymin, double ymax, double ptcutmin, double ptcutmax){
  
  std::cout << "Calculating Amn of charged hadrons ... "  << std::endl ; 
  double pT0 = diskptradius ; // in GeV (disc radius)
  int mmax=8 ; 
  int nmax=8 ; 
  double amnRe[mmax][nmax];
  double amnIm[mmax][nmax];
  double temp_amnRe[mmax][nmax];
  double temp_amnIm[mmax][nmax];
  double lambda[mmax][nmax];
  for(int im=0; im<mmax; im++){
    for(int in=0; in<nmax; in++){
      amnRe[im][in] = 0. ; 
      amnIm[im][in] = 0. ; 
      temp_amnRe[im][in] = 0. ; 
      temp_amnIm[im][in] = 0. ; 
      lambda[im][in] = 0. ; 
    }
  }
  double particles_norm = 0 ; 

  for(int im=0; im<mmax; im++){
    for(int in=0; in<nmax; in++){
      if(in==0){
        continue ;
      }
      lambda[im][in] = gsl_sf_bessel_zero_Jnu(fabs(im),in) ; 
    }
  }


  std::ofstream mFile;
  std::stringstream output_filename;
  output_filename.str("");
  if(yflag > 0 ){
    output_filename << "results/Amnch" ;
    output_filename << "_y_" << ymin << "_" << ymax ;
  }
  else{
    output_filename << "results/Amnch";
    output_filename << "_eta_" << ymin << "_" << ymax ;
  }
  output_filename << "_pt_" << ptcutmin << "_" << ptcutmax ;
  output_filename << "_pt0_" << pT0 ;
  output_filename << ".dat";
  mFile.open(output_filename.str().c_str(), std::ios::out );
  mFile << "# m   n   AmnRe  AmnIm  |Amn| " << std::endl ; 
  
  
  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    if(ii%500 == 0){
      std::cout << "Analysing event " << ii << std::endl ; 
    }

    // set tempAmn for each event to zero at starting of event //
    for(int im=0; im<mmax; im++){
      for(int in=0; in<nmax; in++){
        temp_amnRe[im][in] = 0. ; 
        temp_amnIm[im][in] = 0. ; 
      }
    }

    int ppid = 0 ; 


    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    double temp_particles_norm = 0 ; 
    for(int jj=0; jj<nParticles; jj++){
      ppid = 0 ; 
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ; 
      double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double W   = Event->get_particle(jj)->get_weight() ; 
      double Rap ; 
      
      if( fabs(E-Pz) < 1E-10 )
	continue ; 
      
      if( fabs(P-Pz) < 1E-10 )
	continue ; 
      
      if (yflag > 0 )
	Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
      else
	Rap = 0.5 * TMath::Log( ( P + Pz ) / (P - Pz) );
      
      if(Rap > ymax || Rap < ymin ){
	continue ; 
      }
      
      double Pt = sqrt( Px * Px + Py * Py ) ;
      double Phi = atan2(Py,Px) ;
      if(Pt < ptcutmin || Pt > ptcutmax){
	continue ; 
      }


     if(PID ==  211  || PID == -211  || PID ==  321 || PID == -321  ||
        PID == 2212  || PID == -2212 || PID == 3122 || PID == -3122 ||
        PID == 3222  || PID == -3222 || PID == 3112 || PID == -3112 ||
        PID == 3312  || PID == -3312 || PID == 3334 || PID == -3334  ){
        ppid = 1 ;
      }

      if( ppid < 1 ){
        continue ;
      }

      temp_particles_norm += 1. ; 
      double lambda_mn ; 
      for(int im=0; im<mmax; im++){
	for(int in=0; in<nmax; in++){
	  
	  if(in==0){
	    continue ; 
	  }
	  
	  // get n-th zero of J_m(x) // 
	  lambda_mn =  lambda[im][in] ; 

	  temp_amnRe[im][in] +=  ( 1./gsl_sf_bessel_Jn( fabs(im)+1, lambda_mn ) * gsl_sf_bessel_Jn(im,  Pt / pT0 * lambda_mn ) * cos(im*Phi) )  ; 
	  temp_amnIm[im][in] +=  ( 1./gsl_sf_bessel_Jn( fabs(im)+1, lambda_mn ) * gsl_sf_bessel_Jn(im,  Pt / pT0 * lambda_mn ) * sin(im*Phi) * (-1.0) )  ; 
	} // in loop
      } // im loop
      
    } // particle loop

  for(int im=0; im<mmax; im++){
    for(int in=0; in<nmax; in++){
      amnRe[im][in] += temp_amnRe[im][in] ; 
      amnIm[im][in] += temp_amnIm[im][in] ; 
    }
  }
   particles_norm += temp_particles_norm ; 

  } // event loop
  
  
  for(int im=0; im<mmax; im++){
    for(int in=0; in<nmax; in++){
      amnRe[im][in] /= ( particles_norm * M_PI * pT0 * pT0 ) ; 
      amnIm[im][in] /= ( particles_norm * M_PI * pT0 * pT0 ) ; 
      mFile << im << "  " << in << "  " << amnRe[im][in] << "  " << amnIm[im][in] 
	    << "  " << sqrt(pow(amnRe[im][in],2)+pow(amnIm[im][in],2)) << std::endl ; 
    }
  }
  
  mFile.close();
}



void observables::calculate_amn_from_smeared_grid( int ppid, int yflag, double ymin, double ymax, double ptcutmin, double ptcutmax){
  
  double ptmin = 0.0 ; // in GeV
  double ptmax = 2.0 ; // in GeV
  int npt = 100 ; 
  double dpt = ( ptmax - ptmin ) / npt ; 
  
  double phimin = 0.0 ; // 
  double phimax = 2.0 * M_PI ; // 
  int nphi = 50 ; 
  double dphi = ( phimax - phimin ) / nphi ; 
  
  double fptphi[npt][nphi] ; 
  for(int ii=0; ii<npt; ii++){
    for(int jj=0; jj<nphi; jj++){
      fptphi[ii][jj] = 0. ; 
    }
  }
  
  
  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ; 
      double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double W   = Event->get_particle(jj)->get_weight() ; 
      double Rap ; 
      
      if( fabs(PID-ppid) > 0.1 ){
	continue ; 
      }
      
      if( fabs(E-Pz) < 1E-10 )
	continue ; 
      
      if( fabs(P-Pz) < 1E-10 )
	continue ; 
      
      if (yflag > 0 )
	Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
      else
	Rap = 0.5 * TMath::Log( ( P + Pz ) / (P - Pz) );
      
      if(Rap > ymax || Rap < ymin ){
	continue ; 
      }
      
      double Pt = sqrt( Px * Px + Py * Py ) ;
      double Phi = atan2(Py,Px) ;
      if(Pt < ptcutmin || Pt > ptcutmax){
	continue ; 
      }
      
      int ptbin = floor((Pt - ptmin) / dpt );
      int phibin = floor((Phi - phimin) / dphi ); 
      fptphi[ptbin][phibin] += (1./Pt) ; 
      
    } // particle loop
  } // event loop
  
  
  // Calculate dN / pT dpT dphi // 
  for(int ii=0; ii<npt; ii++){
    for(int jj=0; jj<nphi; jj++){
      fptphi[ii][jj] /= ( nEvents * dpt * dphi )  ; 
    }
  }
  
  // calculate the normalisation constant
  double intval = 0. ; 
  for(int ii=0; ii<npt; ii++){
    for(int jj=0; jj<nphi; jj++){
      double temp_pt = ptmin + (ii * dpt) + (dpt/2.);
      intval += fptphi[ii][jj] * temp_pt * dpt * dphi ; 
    }
  }
  
  // Normalized dN / pT dpT dphi // 
  for(int ii=0; ii<npt; ii++){
    for(int jj=0; jj<nphi; jj++){
      fptphi[ii][jj] /= intval  ; 
    }
  }
  
  
  // Now proceed to Amn calculation //
  double pT0 = diskptradius ; // in GeV (disc radius)
  int mmax=8 ; 
  int nmax=8 ; 
  double amnRe[mmax][nmax];
  double amnIm[mmax][nmax];
  for(int im=0; im<mmax; im++){
    for(int in=0; in<nmax; in++){
      amnRe[im][in] = 0. ; 
      amnIm[im][in] = 0. ; 
    }
  }
  
  
  double lambda_mn ; 
  for(int im=0; im<mmax; im++){
    for(int in=0; in<nmax; in++){
      
      if(in==0){
	continue ; 
      }
      
      // get n-th zero of J_m(x) // 
      lambda_mn = gsl_sf_bessel_zero_Jnu(fabs(im),in) ; 
      for(int ii=0; ii<npt; ii++){
	for(int jj=0; jj<nphi; jj++){
	  double temp_pt = ptmin + (ii * dpt) + (dpt/2.);
	  double temp_phi = phimin + (jj * dphi) + (dphi/2.);
	  amnRe[im][in] +=  (1./M_PI) * (1./pow(pT0,2)) *
                     ( 1./gsl_sf_bessel_Jn( fabs(im)+1, lambda_mn ) * 
			gsl_sf_bessel_Jn(im,  temp_pt / pT0 * lambda_mn ) * 
                             cos(im*temp_phi) * temp_pt * dpt * dphi * fptphi[ii][jj] )  ; 
	  amnIm[im][in] +=  (1./M_PI) * (1./pow(pT0,2)) * 
                      ( 1./gsl_sf_bessel_Jn( fabs(im)+1, lambda_mn ) * 
                          gsl_sf_bessel_Jn(im,  temp_pt / pT0 * lambda_mn ) * 
                              sin(im*temp_phi) * (-1.0) * temp_pt * dpt * dphi * fptphi[ii][jj] )  ; 
	}
      }
      
    } // in loop
  } // im loop
  
  
  
  std::ofstream mFile;
  std::stringstream output_filename;
  output_filename.str("");
  if(yflag > 0 ){
    output_filename << "results/Amn-" << ppid ;
    output_filename << "_y_" << ymin << "_" << ymax ;
  }
  else{
    output_filename << "results/Amn-" << ppid;
    output_filename << "_eta_" << ymin << "_" << ymax ;
  }
  output_filename << "_pt_" << ptcutmin << "_" << ptcutmax ;
  output_filename << "_pt0_" << pT0 ;
  output_filename << "_grid.dat";
  mFile.open(output_filename.str().c_str(), std::ios::out );
  mFile << "# m   n   AmnRe  AmnIm  |Amn| " << std::endl ; 
  
  for(int im=0; im<mmax; im++){
    for(int in=0; in<nmax; in++){
      mFile << im << "  " << in << "  " << amnRe[im][in] << "  " << amnIm[im][in] 
	    << "  " << sqrt(pow(amnRe[im][in],2)+pow(amnIm[im][in],2)) << std::endl ; 
    }
  }
  
  mFile.close();
  
  
}



void observables::calculate_amn_vs_rapidity( int ppid, int yflag, double ptcutmin, double ptcutmax){
  
  std::cout << "Calculating Amn vs rapidity of "  << ppid << " ... " <<  std::endl ; 
  double pT0 = diskptradius ; // in GeV (disc radius)
  int mmax=5 ; 
  int nmax=5 ; 


  double ymin = -3.0 ; // in GeV
  double ymax = 3.0 ; // in GeV
  int ny = 31 ; 
  double dy = ( ymax - ymin ) / ny ;


  double amnRe[ny][mmax][nmax];
  double amnIm[ny][mmax][nmax];
  double temp_amnRe[ny][mmax][nmax];
  double temp_amnIm[ny][mmax][nmax];
  double lambda[mmax][nmax];
  for(int iy=0; iy<ny; iy++){
   for(int im=0; im<mmax; im++){
    for(int in=0; in<nmax; in++){
      amnRe[iy][im][in] = 0. ; 
      amnIm[iy][im][in] = 0. ; 
      temp_amnRe[iy][im][in] = 0. ; 
      temp_amnIm[iy][im][in] = 0. ; 
      lambda[im][in] = 0. ; 
    }
   }
  }

  double particles_norm[ny] ; 
  for(int iy=0; iy<ny; iy++){
    particles_norm[iy] = 0 ; 
  }

  // calculate zeros of Bessel function
  for(int im=0; im<mmax; im++){
    for(int in=0; in<nmax; in++){
      if(in==0){
        continue ;
      }
      lambda[im][in] = gsl_sf_bessel_zero_Jnu(fabs(im),in) ; 
    }
  }

  
  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    if(ii%500 == 0){
      std::cout << "Analysing event " << ii << std::endl ; 
    }

    // set tempAmn for each event to zero at starting of event //
    for(int iy=0; iy<ny; iy++){
     for(int im=0; im<mmax; im++){
      for(int in=0; in<nmax; in++){
        temp_amnRe[iy][im][in] = 0. ; 
        temp_amnIm[iy][im][in] = 0. ; 
      }
     }
    }

    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    double temp_particles_norm[ny] ; 
    for(int iy=0; iy<ny; iy++){
      temp_particles_norm[iy] = 0 ; 
    }
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      if( fabs(PID-ppid) > 0.1 ){
	continue ; 
      }
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ; 
      double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double W   = Event->get_particle(jj)->get_weight() ; 
      double Rap ; 
            
      if( fabs(E-Pz) < 1E-10 )
	continue ; 
      
      if( fabs(P-Pz) < 1E-10 )
	continue ; 
         
      if (yflag > 0 )
	Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
      else
	Rap = 0.5 * TMath::Log( ( P + Pz ) / (P - Pz) );
      
      if(Rap >= ymax || Rap <= ymin ){
	continue ; 
      }

      int ybin = floor((Rap - ymin) / dy ); 
      
      double Pt = sqrt( Px * Px + Py * Py ) ;
      double Phi = atan2(Py,Px) ;
      if(Pt < ptcutmin || Pt > ptcutmax){
	continue ; 
      }
      
      temp_particles_norm[ybin] += 1. ; 
      double lambda_mn ; 
      for(int im=0; im<mmax; im++){
	for(int in=0; in<nmax; in++){
	  if(in==0){
	    continue ; 
	  }
	  
	  // get n-th zero of J_m(x) // 
	  lambda_mn =  lambda[im][in] ; 

	  temp_amnRe[ybin][im][in] +=  ( 1./gsl_sf_bessel_Jn( fabs(im)+1, lambda_mn ) * gsl_sf_bessel_Jn(im,  Pt / pT0 * lambda_mn ) * cos(im*Phi) )  ; 
	  temp_amnIm[ybin][im][in] +=  ( 1./gsl_sf_bessel_Jn( fabs(im)+1, lambda_mn ) * gsl_sf_bessel_Jn(im,  Pt / pT0 * lambda_mn ) * sin(im*Phi) * (-1.0) )  ; 
	} // in loop
      } // im loop
      
    } // particle loop

    for(int iy=0; iy<ny; iy++){
     for(int im=0; im<mmax; im++){
      for(int in=0; in<nmax; in++){
        if(in==0){
          continue ; 
        }
        amnRe[iy][im][in] += temp_amnRe[iy][im][in] ; 
        amnIm[iy][im][in] += temp_amnIm[iy][im][in] ; 
      }
     }
    }

   for(int iy=0; iy<ny; iy++){
     particles_norm[iy] += temp_particles_norm[iy] ;
   } 

  } // event loop
  
  
  for(int iy=0; iy<ny; iy++){
   for(int im=0; im<mmax; im++){
    for(int in=0; in<nmax; in++){
      if(in==0){
        continue ; 
      }
      amnRe[iy][im][in] /= ( particles_norm[iy] * M_PI * pT0 * pT0 ) ; 
      amnIm[iy][im][in] /= ( particles_norm[iy] * M_PI * pT0 * pT0 ) ; 
    }
   }
  }
  
  std::ofstream mFile;
  std::stringstream output_filename;

  for(int im=0; im<mmax; im++){
   for(int in=0; in<nmax; in++){
     if(in==0){
      continue ; 
     }
     output_filename.str("");
     if(yflag > 0 ){
      output_filename << "results/A" << im << in << "_" << ppid ;
     output_filename << "_with_y" ;
    }
    else{
     output_filename << "results/Amn-" << ppid;
     output_filename << "_with_eta" ;
    }
    output_filename << "_pt_" << ptcutmin << "_" << ptcutmax ;
    output_filename << "_pt0_" << pT0 ;
    output_filename << ".dat";
    mFile.open(output_filename.str().c_str(), std::ios::out );
    mFile << "#(m= " << im << "  n= " << in << ")  rap  AmnRe  AmnIm  |Amn| " << std::endl ; 
    for(int iy=0; iy<ny; iy++){  
      mFile << ymin + (iy * dy) + (0.5 * dy) << "  " << amnRe[iy][im][in] << "  " << amnIm[iy][im][in] 
	    << "  " << sqrt(pow(amnRe[iy][im][in],2)+pow(amnIm[iy][im][in],2)) << std::endl ; 
   }
  mFile.close(); 
  }
 }

}









