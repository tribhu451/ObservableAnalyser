#include "observables.h"

observables::observables(){
   eta_min = -8. ; 
   eta_max = 8. ; 
   eta_bins = 48 ; 

   pT_min  = 0 ; 
   pT_max  = 3.0 ; 
   pT_bins = 16 ; 

   y_min = -6.;
   y_max = 6.;
   y_bins = 36;

   double v2_pt_bins[13] = {0.01, 0.1, 0.15, 0.3, 0.5, 0.75,1.0, 1.25, 1.5 , 1.75, 2.1, 2.5, 3.0} ; 

   H1D_DNCHDETA_ETA[0] = new TH1D("H0", "pion_eta_differential_yield", eta_bins, eta_min, eta_max);
   H1D_DNCHDETA_ETA[1] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H1");  H1D_DNCHDETA_ETA[1]->SetTitle("kaon");
   H1D_DNCHDETA_ETA[2] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H2");  H1D_DNCHDETA_ETA[2]->SetTitle("proton");
   H1D_DNCHDETA_ETA[3] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H3");  H1D_DNCHDETA_ETA[3]->SetTitle("anti_proton");
   H1D_DNCHDETA_ETA[4] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H4");  H1D_DNCHDETA_ETA[4]->SetTitle("charged_particle");
   H1D_DNCHDETA_ETA[5] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H5");  H1D_DNCHDETA_ETA[5]->SetTitle("Lambda");
   H1D_DNCHDETA_ETA[6] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H6");  H1D_DNCHDETA_ETA[6]->SetTitle("anti_Lambda");

   PROFILE_V2_PT[0] = new TProfile("PH0", "pion_pT_differential_v2", 12, v2_pt_bins);
   PROFILE_V2_PT[1] = (TProfile*) PROFILE_V2_PT[0]->Clone("PH1");  PROFILE_V2_PT[1]->SetTitle("kaon");
   PROFILE_V2_PT[2] = (TProfile*) PROFILE_V2_PT[0]->Clone("PH2");  PROFILE_V2_PT[2]->SetTitle("proton");
   PROFILE_V2_PT[3] = (TProfile*) PROFILE_V2_PT[0]->Clone("PH3");  PROFILE_V2_PT[3]->SetTitle("anti_proton");
   PROFILE_V2_PT[4] = (TProfile*) PROFILE_V2_PT[0]->Clone("PH4");  PROFILE_V2_PT[4]->SetTitle("charged_particle");

   PROFILE_V1_Y[0] = new TProfile("PV0", "pion_plus_v1_y", y_bins, y_min, y_max);
   PROFILE_V1_Y[1] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV1");  PROFILE_V1_Y[1]->SetTitle("pion_minus");
   PROFILE_V1_Y[2] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV2");  PROFILE_V1_Y[2]->SetTitle("proton");
   PROFILE_V1_Y[3] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV3");  PROFILE_V1_Y[3]->SetTitle("anti_proton");
   PROFILE_V1_Y[4] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV4");  PROFILE_V1_Y[4]->SetTitle("charged_particle");

   PROFILE_V1_ETA[0] = new TProfile("PX0", "charged_particle_v1_eta", eta_bins, eta_min, eta_max);

   H1D_INVYLD_PT[0] = new TH1D("PT0", "pion_pT_differential_v2", pT_bins, pT_min, pT_max);
   H1D_INVYLD_PT[1] = (TProfile*) H1D_INVYLD_PT[0]->Clone("PT1");  H1D_INVYLD_PT[1]->SetTitle("kaon");
   H1D_INVYLD_PT[2] = (TProfile*) H1D_INVYLD_PT[0]->Clone("PT2");  H1D_INVYLD_PT[2]->SetTitle("proton");
   H1D_INVYLD_PT[3] = (TProfile*) H1D_INVYLD_PT[0]->Clone("PT3");  H1D_INVYLD_PT[3]->SetTitle("anti_proton");
   H1D_INVYLD_PT[4] = (TProfile*) H1D_INVYLD_PT[0]->Clone("PT4");  H1D_INVYLD_PT[4]->SetTitle("Lambda");
   H1D_INVYLD_PT[5] = (TProfile*) H1D_INVYLD_PT[0]->Clone("PT5");  H1D_INVYLD_PT[5]->SetTitle("anti_Lambda");

}


observables::~observables(){
  delete[] H1D_DNCHDETA_ETA  ; 
  delete[] PROFILE_V2_PT  ; 
  delete[] PROFILE_V1_Y  ; 
  delete[] PROFILE_V1_ETA  ; 
  delete[] H1D_INVYLD_PT ; 
}

void observables::fill_histogram_of_dnchdeta_eta(events* Event, double pT_min, double pT_max ){

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
      if(PID == 211){
        H1D_DNCHDETA_ETA[0]->Fill(Eta,1.);
        H1D_DNCHDETA_ETA[4]->Fill(Eta,1.);
      }
      if(PID == 321){
        H1D_DNCHDETA_ETA[1]->Fill(Eta,1.);
        H1D_DNCHDETA_ETA[4]->Fill(Eta,1.);
      }
      if(PID == 2212){
        H1D_DNCHDETA_ETA[2]->Fill(Eta,1.);
        H1D_DNCHDETA_ETA[4]->Fill(Eta,1.);
      }
      if(PID == -2212){
        H1D_DNCHDETA_ETA[3]->Fill(Eta,1.);
        H1D_DNCHDETA_ETA[4]->Fill(Eta,1.);
      }
      if(PID == 3122){
        H1D_DNCHDETA_ETA[5]->Fill(Eta,1.);
        H1D_DNCHDETA_ETA[4]->Fill(Eta,1.);
      }
      if(PID == -3122){
        H1D_DNCHDETA_ETA[6]->Fill(Eta,1.);
        H1D_DNCHDETA_ETA[4]->Fill(Eta,1.);
      }
     if(PID == -211 || PID == -321 || 
        PID == 3222 || PID == -3222 || PID == 3112 || PID == -3112 ||
        PID == 3312 || PID == -3312 || PID == 3334 || PID == -3334  ){
        H1D_DNCHDETA_ETA[4]->Fill(Eta,1.);
      }

  
    }
}


void observables::calculate_dnchdeta_eta(int nEvents){
  double dX = ( eta_max - eta_min ) / eta_bins ; 
  for(int ii=0; ii<_N_HISTOGRAMS_DNCHDETA_ETA; ii++)
        H1D_DNCHDETA_ETA[ii]->Scale(1.0 / (nEvents * dX));

  std::ofstream mFile;
  mFile.open("results/Dnchdeta_eta-211.dat");
  for(int i=1; i< H1D_DNCHDETA_ETA[0]->GetNbinsX(); i++){
      mFile << H1D_DNCHDETA_ETA[0]->GetBinCenter(i) << "\t" << H1D_DNCHDETA_ETA[0]->GetBinContent(i) 
             << "\t" << H1D_DNCHDETA_ETA[0]->GetBinError(i) << std::endl;
  }
  mFile.close();

  mFile.open("results/Dnchdeta_eta-321.dat");
  for(int i=1; i< H1D_DNCHDETA_ETA[1]->GetNbinsX(); i++){
      mFile << H1D_DNCHDETA_ETA[1]->GetBinCenter(i) << "\t" << H1D_DNCHDETA_ETA[1]->GetBinContent(i) 
             << "\t" << H1D_DNCHDETA_ETA[1]->GetBinError(i) << std::endl;
  }
  mFile.close();

  mFile.open("results/Dnchdeta_eta-2212.dat");
  for(int i=1; i< H1D_DNCHDETA_ETA[2]->GetNbinsX(); i++){
      mFile << H1D_DNCHDETA_ETA[2]->GetBinCenter(i) << "\t" << H1D_DNCHDETA_ETA[2]->GetBinContent(i) 
             << "\t" << H1D_DNCHDETA_ETA[2]->GetBinError(i) << std::endl;
  }
  mFile.close();

  mFile.open("results/Dnchdeta_eta--2212.dat");
  for(int i=1; i< H1D_DNCHDETA_ETA[3]->GetNbinsX(); i++){
      mFile << H1D_DNCHDETA_ETA[3]->GetBinCenter(i) << "\t" << H1D_DNCHDETA_ETA[3]->GetBinContent(i) 
             << "\t" << H1D_DNCHDETA_ETA[3]->GetBinError(i) << std::endl;
  }
  mFile.close();

  mFile.open("results/Dnchdeta_eta-charged_particles.dat");
  for(int i=1; i< H1D_DNCHDETA_ETA[4]->GetNbinsX(); i++){
      mFile << H1D_DNCHDETA_ETA[4]->GetBinCenter(i) << "\t" << H1D_DNCHDETA_ETA[4]->GetBinContent(i) 
             << "\t" << H1D_DNCHDETA_ETA[4]->GetBinError(i) << std::endl;
  }
  mFile.close();


  mFile.open("results/Dnchdeta_eta-3122.dat");
  for(int i=1; i< H1D_DNCHDETA_ETA[5]->GetNbinsX(); i++){
      mFile << H1D_DNCHDETA_ETA[5]->GetBinCenter(i) << "\t" << H1D_DNCHDETA_ETA[5]->GetBinContent(i) 
             << "\t" << H1D_DNCHDETA_ETA[5]->GetBinError(i) << std::endl;
  }
  mFile.close();


  mFile.open("results/Dnchdeta_eta--3122.dat");
  for(int i=1; i< H1D_DNCHDETA_ETA[6]->GetNbinsX(); i++){
      mFile << H1D_DNCHDETA_ETA[6]->GetBinCenter(i) << "\t" << H1D_DNCHDETA_ETA[6]->GetBinContent(i) 
             << "\t" << H1D_DNCHDETA_ETA[6]->GetBinError(i) << std::endl;
  }
  mFile.close();

}


void observables::fill_histogram_of_v2_pt(events* Event, double Rap_min, double Rap_max ){

  int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ; 
      double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
      double Pt = sqrt( Px * Px + Py * Py )  ;
      double v2 = ( Px * Px - Py * Py ) / ( Pt * Pt ) ; 
      if(Rap > Rap_max || Rap < Rap_min)
        continue ; 
      if(PID == 211){
        PROFILE_V2_PT[0]->Fill(Pt,v2);
      }
      if(PID == 321){
        PROFILE_V2_PT[1]->Fill(Pt,v2);
      }
      if(PID == 2212){
        PROFILE_V2_PT[2]->Fill(Pt,v2);
      }
      if(PID == -2212){
        PROFILE_V2_PT[3]->Fill(Pt,v2);
      }
      if( PID == 211   || PID == -211  || PID ==  321  || PID == -321  || 
          PID ==  2212 || PID == -2212 || PID ==  3122 || PID == -3122 ||
          PID == 3222  || PID == -3222 || PID == 3112  || PID == -3112 ||
          PID == 3312  || PID == -3312 || PID == 3334  || PID == -3334    ){
            PROFILE_V2_PT[4]->Fill(Pt,v2);
      }
  
    } // particle loop
}

void observables::calculate_v2_pt(){

  std::ofstream mFile;
  mFile.open("results/v2_pT-211.dat");
  for(int i=1; i< PROFILE_V2_PT[0]->GetNbinsX(); i++){
      mFile << PROFILE_V2_PT[0]->GetBinCenter(i) << "\t" << PROFILE_V2_PT[0]->GetBinContent(i) 
             << "\t" << PROFILE_V2_PT[0]->GetBinError(i) << std::endl;
  }
  mFile.close();

  mFile.open("results/v2_pT-321.dat");
  for(int i=1; i< PROFILE_V2_PT[1]->GetNbinsX(); i++){
      mFile << PROFILE_V2_PT[1]->GetBinCenter(i) << "\t" << PROFILE_V2_PT[1]->GetBinContent(i) 
             << "\t" << PROFILE_V2_PT[1]->GetBinError(i) << std::endl;
  }
  mFile.close();


  mFile.open("results/v2_pT-2212.dat");
  for(int i=1; i< PROFILE_V2_PT[2]->GetNbinsX(); i++){
      mFile << PROFILE_V2_PT[2]->GetBinCenter(i) << "\t" << PROFILE_V2_PT[2]->GetBinContent(i) 
             << "\t" << PROFILE_V2_PT[2]->GetBinError(i) << std::endl;
  }
  mFile.close();

  mFile.open("results/v2_pT--2212.dat");
  for(int i=1; i< PROFILE_V2_PT[3]->GetNbinsX(); i++){
      mFile << PROFILE_V2_PT[3]->GetBinCenter(i) << "\t" << PROFILE_V2_PT[3]->GetBinContent(i) 
             << "\t" << PROFILE_V2_PT[3]->GetBinError(i) << std::endl;
  }
  mFile.close();

  mFile.open("results/v2_pT-charged-particle.dat");
  for(int i=1; i< PROFILE_V2_PT[4]->GetNbinsX(); i++){
      mFile << PROFILE_V2_PT[4]->GetBinCenter(i) << "\t" << PROFILE_V2_PT[4]->GetBinContent(i) 
             << "\t" << PROFILE_V2_PT[4]->GetBinError(i) << std::endl;
  }
  mFile.close();
}






void observables::fill_histogram_of_v1_y(events* Event, double pT_min, double pT_max ){

  int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ; 
      double Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
      double Pt = sqrt( Px * Px + Py * Py )  ;
      if(Pt > pT_max || Pt < pT_min)
        continue ; 

      double v1 = Px / Pt ; 
      if(PID == 211){
        PROFILE_V1_Y[0]->Fill(Rap,v1);
      }
      if(PID == -211){
        PROFILE_V1_Y[1]->Fill(Rap,v1);
      }
      if(PID == 2212){
        PROFILE_V1_Y[2]->Fill(Rap,v1);
      }
      if(PID == -2212){
        PROFILE_V1_Y[3]->Fill(Rap,v1);
      }
      if( PID == 211   || PID == -211  || PID ==  321  || PID == -321  || 
          PID ==  2212 || PID == -2212 || PID ==  3122 || PID == -3122 ||
          PID == 3222  || PID == -3222 || PID == 3112  || PID == -3112 ||
          PID == 3312  || PID == -3312 || PID == 3334  || PID == -3334    ){
           PROFILE_V1_Y[4]->Fill(Rap,v1);
      }
  
    } // particle loop
}


void observables::calculate_v1_y(){
  std::ofstream mFile;
  mFile.open("results/v1_y-211.dat");
  for(int i=1; i< PROFILE_V1_Y[0]->GetNbinsX(); i++){
      mFile << PROFILE_V1_Y[0]->GetBinCenter(i) << "\t" << PROFILE_V1_Y[0]->GetBinContent(i) 
             << "\t" << PROFILE_V1_Y[0]->GetBinError(i) << std::endl;
  }
  mFile.close();

  mFile.open("results/v1_y--211.dat");
  for(int i=1; i< PROFILE_V1_Y[1]->GetNbinsX(); i++){
      mFile << PROFILE_V1_Y[1]->GetBinCenter(i) << "\t" << PROFILE_V1_Y[1]->GetBinContent(i) 
             << "\t" << PROFILE_V1_Y[1]->GetBinError(i) << std::endl;
  }
  mFile.close();

  mFile.open("results/v1_y-2212.dat");
  for(int i=1; i< PROFILE_V1_Y[2]->GetNbinsX(); i++){
      mFile << PROFILE_V1_Y[2]->GetBinCenter(i) << "\t" << PROFILE_V1_Y[2]->GetBinContent(i) 
             << "\t" << PROFILE_V1_Y[2]->GetBinError(i) << std::endl;
  }
  mFile.close();

  mFile.open("results/v1_y--2212.dat");
  for(int i=1; i< PROFILE_V1_Y[3]->GetNbinsX(); i++){
      mFile << PROFILE_V1_Y[3]->GetBinCenter(i) << "\t" << PROFILE_V1_Y[3]->GetBinContent(i) 
             << "\t" << PROFILE_V1_Y[3]->GetBinError(i) << std::endl;
  }
  mFile.close();

  mFile.open("results/v1_y--charged_particle.dat");
  for(int i=1; i< PROFILE_V1_Y[4]->GetNbinsX(); i++){
      mFile << PROFILE_V1_Y[4]->GetBinCenter(i) << "\t" << PROFILE_V1_Y[4]->GetBinContent(i) 
             << "\t" << PROFILE_V1_Y[4]->GetBinError(i) << std::endl;
  }
  mFile.close();

}



void observables::fill_histogram_of_v1_eta(events* Event, double pT_min, double pT_max ){

  int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double Rap = 0.5 * TMath::Log( ( P + Pz ) / (P - Pz) );
      double Pt = sqrt( Px * Px + Py * Py )  ;
      if(Pt > pT_max || Pt < pT_min)
        continue ; 

      double v1 = Px / Pt ; 
      if( PID == 211   || PID == -211  || PID ==  321  || PID == -321  || 
          PID ==  2212 || PID == -2212 || PID ==  3122 || PID == -3122 ||
          PID == 3222  || PID == -3222 || PID == 3112  || PID == -3112 ||
          PID == 3312  || PID == -3312 || PID == 3334  || PID == -3334    ){
           PROFILE_V1_ETA[0]->Fill(Rap,v1);
      }
  
    } // particle loop
}

void observables::calculate_v1_eta(){
  std::ofstream mFile;
  mFile.open("results/v1_eta-charged_particle.dat");
  for(int i=1; i< PROFILE_V1_ETA[0]->GetNbinsX(); i++){
      mFile << PROFILE_V1_ETA[0]->GetBinCenter(i) << "\t" << PROFILE_V1_ETA[0]->GetBinContent(i) 
             << "\t" << PROFILE_V1_ETA[0]->GetBinError(i) << std::endl;
  }
  mFile.close();

}



void observables::fill_histogram_of_invariant_yield_pt(events* Event, double Rap_min, double Rap_max ){

  int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ; 
      double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
      double Pt = sqrt( Px * Px + Py * Py )  ;
      if(Rap > Rap_max || Rap < Rap_min)
        continue ; 
      if(PID == 211){
        H1D_INVYLD_PT[0]->Fill(Pt,1.0/Pt);
      }
      if(PID == 321){
        H1D_INVYLD_PT[1]->Fill(Pt,1.0/Pt);
      }
      if(PID == 2212){
        H1D_INVYLD_PT[2]->Fill(Pt,1.0/Pt);
      }
      if(PID == -2212){
        H1D_INVYLD_PT[3]->Fill(Pt,1.0/Pt);
      }
      if(PID == 3122){
        H1D_INVYLD_PT[4]->Fill(Pt,1.0/Pt);
      }
      if(PID == -3122){
        H1D_INVYLD_PT[5]->Fill(Pt,1.0/Pt);
      }
    } // particle loop
}

void observables::calculate_invariant_yield_pt(int Nevents, double Rap_min, double Rap_max){

  double dpT = ( pT_max - pT_min ) / pT_bins ; 
  double dX0 = ( eta_max - eta_min )  ; 
  H1D_INVYLD_PT[0]->Scale( 1.0 / ( Nevents  * 2.0 * TMath::Pi() * dX0 * dpT ) );
  H1D_INVYLD_PT[1]->Scale( 1.0 / ( Nevents  * 2.0 * TMath::Pi() * dX0 * dpT ) );
  H1D_INVYLD_PT[2]->Scale( 1.0 / ( Nevents  * 2.0 * TMath::Pi() * dX0 * dpT ) );
  H1D_INVYLD_PT[3]->Scale( 1.0 / ( Nevents  * 2.0 * TMath::Pi() * dX0 * dpT ) );
  H1D_INVYLD_PT[4]->Scale( 1.0 / ( Nevents  * 2.0 * TMath::Pi() * dX0 * dpT ) );
  H1D_INVYLD_PT[5]->Scale( 1.0 / ( Nevents  * 2.0 * TMath::Pi() * dX0 * dpT ) );

  std::ofstream mFile;
  mFile.open("results/invariant_yield_pt-211.dat");
  for(int i=1; i< H1D_INVYLD_PT[0]->GetNbinsX(); i++){
      mFile << H1D_INVYLD_PT[0]->GetBinCenter(i) << "\t" << H1D_INVYLD_PT[0]->GetBinContent(i) 
             << "\t" << H1D_INVYLD_PT[0]->GetBinError(i) << std::endl;
  }
  mFile.close();

  mFile.open("results/invariant_yield_pt-321.dat");
  for(int i=1; i< H1D_INVYLD_PT[1]->GetNbinsX(); i++){
      mFile << H1D_INVYLD_PT[1]->GetBinCenter(i) << "\t" << H1D_INVYLD_PT[1]->GetBinContent(i) 
             << "\t" << H1D_INVYLD_PT[1]->GetBinError(i) << std::endl;
  }
  mFile.close();

  mFile.open("results/invariant_yield_pt-2212.dat");
  for(int i=1; i< H1D_INVYLD_PT[2]->GetNbinsX(); i++){
      mFile << H1D_INVYLD_PT[2]->GetBinCenter(i) << "\t" << H1D_INVYLD_PT[2]->GetBinContent(i) 
             << "\t" << H1D_INVYLD_PT[2]->GetBinError(i) << std::endl;
  }
  mFile.close();

  mFile.open("results/invariant_yield_pt--2212.dat");
  for(int i=1; i< H1D_INVYLD_PT[3]->GetNbinsX(); i++){
      mFile << H1D_INVYLD_PT[3]->GetBinCenter(i) << "\t" << H1D_INVYLD_PT[3]->GetBinContent(i) 
             << "\t" << H1D_INVYLD_PT[3]->GetBinError(i) << std::endl;
  }
  mFile.close();

  mFile.open("results/invariant_yield_pt-3122.dat");
  for(int i=1; i< H1D_INVYLD_PT[4]->GetNbinsX(); i++){
      mFile << H1D_INVYLD_PT[4]->GetBinCenter(i) << "\t" << H1D_INVYLD_PT[4]->GetBinContent(i) 
             << "\t" << H1D_INVYLD_PT[4]->GetBinError(i) << std::endl;
  }
  mFile.close();

  mFile.open("results/invariant_yield_pt--3122.dat");
  for(int i=1; i< H1D_INVYLD_PT[5]->GetNbinsX(); i++){
      mFile << H1D_INVYLD_PT[5]->GetBinCenter(i) << "\t" << H1D_INVYLD_PT[5]->GetBinContent(i) 
             << "\t" << H1D_INVYLD_PT[5]->GetBinError(i) << std::endl;
  }
  mFile.close();




}


























