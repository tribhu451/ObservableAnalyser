#include "observables.h"

observables::observables(input_paramters &iparam_, read_input_file* rif_) : iparam(iparam_), rif(rif_) {
}


void observables::calculate_dnchdeta_eta( double pT_min, double pT_max ){

  const int _N_HISTOGRAMS_DNCHDETA_ETA = 9 ; 

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

  int nEvents = iparam.nEvents ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      //double E   = Event->get_particle(jj)->get_e()   ; 
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
        H1D_DNCHDETA_ETA[1]->Fill(Eta,1.);
      }
      if(PID == -211){
        H1D_DNCHDETA_ETA[2]->Fill(Eta,1.);
      }
      if(PID == 321){
        H1D_DNCHDETA_ETA[3]->Fill(Eta,1.);
      }
      if(PID == -321){
        H1D_DNCHDETA_ETA[4]->Fill(Eta,1.);
      }
      if(PID == 2212){
        H1D_DNCHDETA_ETA[5]->Fill(Eta,1.);
      }
      if(PID == -2212){
        H1D_DNCHDETA_ETA[6]->Fill(Eta,1.);
      }
      if(PID == 3122){
        H1D_DNCHDETA_ETA[7]->Fill(Eta,1.);
      }
      if(PID == -3122){
        H1D_DNCHDETA_ETA[8]->Fill(Eta,1.);
      }
     } // particle loop
    } // event loop


  double dX = ( Eta_max - Eta_min ) / Eta_bins ; 
  for(int ii=0; ii<_N_HISTOGRAMS_DNCHDETA_ETA; ii++)
        H1D_DNCHDETA_ETA[ii]->Scale(1.0 / (nEvents * dX));

  std::ofstream mFile;
  std::stringstream output_filename;

  std::string hadron_name[_N_HISTOGRAMS_DNCHDETA_ETA] = {"hpm", "211", "-211", "321", "-321", "2212", "-2212", "3122", "-3122" };
  int        hadron_index[_N_HISTOGRAMS_DNCHDETA_ETA] = {  0,     1 ,     2 ,    3 ,     4 ,     5 ,      6 ,     7 ,      8   };


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

  const int _N_HISTOGRAMS_DNDY_Y = 9 ; 

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

  int nEvents = iparam.nEvents ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ; 
      //double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
      double Pt = sqrt( Px * Px + Py * Py )  ;
      if(Pt > pT_max || Pt < pT_min)
        continue ; 

     if(PID ==  211  || PID == -211  || PID ==  321 || PID == -321  || 
        PID == 2212  || PID == -2212 || PID == 3122 || PID == -3122 ||
        PID == 3222  || PID == -3222 || PID == 3112 || PID == -3112 ||
        PID == 3312  || PID == -3312 || PID == 3334 || PID == -3334  ){
        H1D_DNDY_Y[0]->Fill(Rap,1.);
      }

      if(PID == 211){
        H1D_DNDY_Y[1]->Fill(Rap,1.);
      }
      if(PID == -211){
        H1D_DNDY_Y[2]->Fill(Rap,1.);
      }
      if(PID == 321){
        H1D_DNDY_Y[3]->Fill(Rap,1.);
      }
      if(PID == -321){
        H1D_DNDY_Y[4]->Fill(Rap,1.);
      }
      if(PID == 2212){
        H1D_DNDY_Y[5]->Fill(Rap,1.);
      }
      if(PID == -2212){
        H1D_DNDY_Y[6]->Fill(Rap,1.);
      }
      if(PID == 3122){
        H1D_DNDY_Y[7]->Fill(Rap,1.);
      }
      if(PID == -3122){
        H1D_DNDY_Y[8]->Fill(Rap,1.);
      }
     } // particle loop
    } // event loop


  double dX = ( Y_max - Y_min ) / Y_bins ; 
  for(int ii=0; ii<_N_HISTOGRAMS_DNDY_Y; ii++)
        H1D_DNDY_Y[ii]->Scale(1.0 / (nEvents * dX));

  std::ofstream mFile;
  std::stringstream output_filename;

  std::string hadron_name[_N_HISTOGRAMS_DNDY_Y] = {"hpm", "211", "-211", "321", "-321", "2212", "-2212", "3122", "-3122" };
  int        hadron_index[_N_HISTOGRAMS_DNDY_Y] = {  0,     1 ,     2 ,    3 ,     4 ,     5 ,      6 ,     7 ,      8   };


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

  const int _N_HISTOGRAMS_INVYLD_PT_ = 8 ; 

  double pT_min  = 0.01 ; 
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


  int nEvents = iparam.nEvents ; 
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
        H1D_INVYLD_PT[0]->Fill(Pt,1.0/Pt);
      }
      if(PID == -211){
        H1D_INVYLD_PT[1]->Fill(Pt,1.0/Pt);
      }
      if(PID == 321){
        H1D_INVYLD_PT[2]->Fill(Pt,1.0/Pt);
      }
      if(PID == -321){
        H1D_INVYLD_PT[3]->Fill(Pt,1.0/Pt);
      }
      if(PID == 2212){
        H1D_INVYLD_PT[4]->Fill(Pt,1.0/Pt);
      }
      if(PID == -2212){
        H1D_INVYLD_PT[5]->Fill(Pt,1.0/Pt);
      }
      if(PID == 3122){
        H1D_INVYLD_PT[6]->Fill(Pt,1.0/Pt);
      }
      if(PID == -3122){
        H1D_INVYLD_PT[7]->Fill(Pt,1.0/Pt);
      }
     } // particle loop
    } // event loop


  H1D_INVYLD_PT[0]->Scale( 1.0 / ( nEvents  * 2.0 * TMath::Pi() * dY * dpT ) );
  H1D_INVYLD_PT[1]->Scale( 1.0 / ( nEvents  * 2.0 * TMath::Pi() * dY * dpT ) );
  H1D_INVYLD_PT[2]->Scale( 1.0 / ( nEvents  * 2.0 * TMath::Pi() * dY * dpT ) );
  H1D_INVYLD_PT[3]->Scale( 1.0 / ( nEvents  * 2.0 * TMath::Pi() * dY * dpT ) );
  H1D_INVYLD_PT[4]->Scale( 1.0 / ( nEvents  * 2.0 * TMath::Pi() * dY * dpT ) );
  H1D_INVYLD_PT[5]->Scale( 1.0 / ( nEvents  * 2.0 * TMath::Pi() * dY * dpT ) );

  std::ofstream mFile;
  std::stringstream output_filename;

  std::string hadron_name[_N_HISTOGRAMS_INVYLD_PT_] = { "211", "-211", "321", "-321", "2212", "-2212", "3122", "-3122" };
  int        hadron_index[_N_HISTOGRAMS_INVYLD_PT_] = {   0,      1 ,     2 ,    3 ,     4 ,     5 ,      6 ,     7    };


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

   double Y_min   = -8. ; 
   double Y_max   =  8. ; 
   int    Y_bins  =  48 ; 
 
   const int   _N_HISTOGRAMS_V1_Y_ = 9 ; 
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

  int nEvents = iparam.nEvents ; 
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
        PROFILE_V1_Y[0]->Fill(Rap,v1);
      }

      if(PID == 211){
        PROFILE_V1_Y[1]->Fill(Rap,v1);
      }
      if(PID == -211){
        PROFILE_V1_Y[2]->Fill(Rap,v1);
      }
      if(PID == 321){
        PROFILE_V1_Y[3]->Fill(Rap,v1);
      }
      if(PID == -321){
        PROFILE_V1_Y[4]->Fill(Rap,v1);
      }
      if(PID == 2212){
        PROFILE_V1_Y[5]->Fill(Rap,v1);
      }
      if(PID == -2212){
        PROFILE_V1_Y[6]->Fill(Rap,v1);
      }
      if(PID == 3122){
        PROFILE_V1_Y[7]->Fill(Rap,v1);
      }
      if(PID == -3122){
        PROFILE_V1_Y[8]->Fill(Rap,v1);
      }

     } // particle loop
    } // event loop

  std::ofstream mFile;
  std::stringstream output_filename;

  std::string hadron_name[_N_HISTOGRAMS_V1_Y_] = {"hpm", "211", "-211", "321", "-321", "2212", "-2212", "3122", "-3122" };
  int        hadron_index[_N_HISTOGRAMS_V1_Y_] = {  0,     1 ,     2 ,    3 ,     4 ,     5 ,      6 ,     7 ,      8   };


 for(int ix =0; ix < _N_HISTOGRAMS_V1_Y_; ix++){
   output_filename.str("");

   if(yflag > 0 ){
     output_filename << "results/v1_y-" << hadron_name[ix];
     output_filename << "_pT_" << pT_min << "_" << pT_max ;
   }
  else{
     output_filename << "results/v1_y-" << hadron_name[ix];
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

   for(int ixx=0; ixx < _N_HISTOGRAMS_V1_Y_ ; ixx++){
      PROFILE_V1_Y[ixx]->Clear(); 
   }


}



void observables::calculate_v2_pt( int yflag, double Rap_min, double Rap_max ){

  const int _N_HISTOGRAMS_V2_PT = 9 ; 
  double v2_pt_bins[13] = {0.01, 0.1, 0.15, 0.3, 0.5, 0.75,1.0, 1.25, 1.5 , 1.75, 2.1, 2.5, 3.0} ; 

  TProfile*                 PROFILE_V2_PT[_N_HISTOGRAMS_V2_PT] ; 
  PROFILE_V2_PT[0] = new TProfile("V2PT0", "pion_pT_differential_v2", 12, v2_pt_bins);
  PROFILE_V2_PT[1] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT0");  PROFILE_V2_PT[1]->SetTitle("pion_plus");
  PROFILE_V2_PT[2] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT1");  PROFILE_V2_PT[2]->SetTitle("pion_minus");
  PROFILE_V2_PT[3] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT2");  PROFILE_V2_PT[3]->SetTitle("kaon_plus");
  PROFILE_V2_PT[4] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT3");  PROFILE_V2_PT[4]->SetTitle("kaon_minus");
  PROFILE_V2_PT[5] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT4");  PROFILE_V2_PT[5]->SetTitle("proton");
  PROFILE_V2_PT[6] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT5");  PROFILE_V2_PT[6]->SetTitle("anti_proton");
  PROFILE_V2_PT[7] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT6");  PROFILE_V2_PT[7]->SetTitle("Lambda");
  PROFILE_V2_PT[8] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT7");  PROFILE_V2_PT[8]->SetTitle("anti_Lambda");

  int nEvents = iparam.nEvents ; 
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
        PROFILE_V2_PT[0]->Fill(Pt,v2);
      }

      if(PID == 211){
        PROFILE_V2_PT[1]->Fill(Pt,v2);
      }
      if(PID == -211){
        PROFILE_V2_PT[2]->Fill(Pt,v2);
      }
      if(PID == 321){
        PROFILE_V2_PT[3]->Fill(Pt,v2);
      }
      if(PID == -321){
        PROFILE_V2_PT[4]->Fill(Pt,v2);
      }
      if(PID == 2212){
        PROFILE_V2_PT[5]->Fill(Pt,v2);
      }
      if(PID == -2212){
        PROFILE_V2_PT[6]->Fill(Pt,v2);
      }
      if(PID == 3122){
        PROFILE_V2_PT[7]->Fill(Pt,v2);
      }
      if(PID == -3122){
        PROFILE_V2_PT[8]->Fill(Pt,v2);
      }
     } // particle loop
    } // event loop


  std::ofstream mFile;
  std::stringstream output_filename;

  std::string hadron_name[_N_HISTOGRAMS_V2_PT] = {"hpm", "211", "-211", "321", "-321", "2212", "-2212", "3122", "-3122" };
  int        hadron_index[_N_HISTOGRAMS_V2_PT] = {  0,     1 ,     2 ,    3 ,     4 ,     5 ,      6 ,     7 ,      8   };


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

























