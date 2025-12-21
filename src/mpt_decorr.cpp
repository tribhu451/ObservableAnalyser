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
  deta = ( etamax - etamin ) / Neta ; 
  eta_bin_centers = new float[Neta];
  for (int ii = 0; ii < Neta; ii++) {
     eta_bin_centers[ii] = etamin + (ii+0.5) * deta ;
     std::cout << etamin + (ii) * deta  << "   " << eta_bin_centers[ii] 
     << "   " <<  etamin + (ii+1) * deta  << std::endl ; 
   }

  // storing [pt] in event x EtaBin matrix 
  IptI = new float[Nevents*Neta]; // (i-event, j-etabin)

  double temp_mpt[Neta];
  double temp_mpt_count[Neta]; 

  // calculate meanpt of each event at each rapidity window and store it.
  for(int ii=0; ii<Nevents; ii++){
    events* Event = rif->get_event(ii) ;
    for(int ieta=0; ieta<Neta; ieta++){
      temp_mpt[ieta] = 0. ; 
      temp_mpt_count[ieta] = 0. ; 
    } 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ; 
      double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double Pt = sqrt( Px * Px + Py * Py ) ;
      double Rap ; 

      if( fabs(E-Pz) < 1E-10 )
        continue ; 

      if( fabs(P-Pz) < 1E-10 )
        continue ; 

      if (yflag > 0 )
        Rap = 0.5 * log( ( E + Pz ) / (E - Pz) );
      else
        Rap = 0.5 * log( ( P + Pz ) / (P - Pz) );

      if(Pt > ptmax || Pt < ptmin)
        continue ; 

      if(Rap > etamax || Rap < etamin)
        continue ; 

      int RapBin = get_rap_bin_index(Rap);

      if(part==0){ // if charged particle
       if( PID==211 || PID==-211 || PID==321 || PID==-321 || PID==2212 || PID==-2212 ){
         temp_mpt[RapBin] += Pt;
         temp_mpt_count[RapBin] += 1. ; 
       }
       else{
         continue ; 
       }
      }
      else if(part==PID){
        temp_mpt[RapBin] += Pt;
        temp_mpt_count[RapBin] += 1. ; 
      }
      else{
       continue ; 
      }

    } // particle loop

    for(int ieta=0; ieta<Neta; ieta++){
      temp_mpt[ieta] /= temp_mpt_count[ieta] ; 
      set_IptI(ii, ieta, temp_mpt[ieta]);
    }

  } // event loop

}


mpt_decorr::~mpt_decorr(){
  delete[] eta_bin_centers ; 
  delete[] IptI ; 
}

std::vector<int> mpt_decorr::get_an_event_ensemble(){
  std::vector<int> event_ID_ens;
  for(int ii=0; ii<Nevents; ii++){
    int evID =  Nevents * rand->rand_uniform() ;
    event_ID_ens.push_back(evID);
  }
  return event_ID_ens;
}

void mpt_decorr::calculate_mean_and_variance_of_pt_of_one_ensemble(int ieta,
      std::vector<int> event_ID_ens, double&  Mpt, double&  var_Mpt){
  double sum_mpt = 0. ; 
  double sum_mpt_sq = 0. ; 
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ;
    sum_mpt += get_IptI(eventID,ieta); 
    sum_mpt_sq += get_IptI(eventID,ieta)*get_IptI(eventID,ieta); 
  }
  sum_mpt /= event_ID_ens.size() ; 
  sum_mpt_sq /= event_ID_ens.size() ; 
  Mpt = sum_mpt ; 
  var_Mpt = sum_mpt_sq - sum_mpt * sum_mpt ; 
}

void mpt_decorr::write_mean_and_variance_of_pt_with_eta(){

  std::vector<int> event_ID_ens;
  double sumMean[Neta];
  double sumMeanSq[Neta];
  double sumVar[Neta];
  double sumVarSq[Neta];
  double temp_Mpt, temp_varMpt ; 

  for(int ieta=0; ieta<Neta; ieta++){
   sumMean[ieta]   = 0. ; 
   sumMeanSq[ieta] = 0. ; 
   sumVar[ieta]    = 0. ; 
   sumVarSq[ieta]  = 0. ; 
  }

  for(int ievents=0; ievents<Nevents; ievents++){
    event_ID_ens = get_an_event_ensemble();
    for(int ieta=0; ieta<Neta; ieta++){
      calculate_mean_and_variance_of_pt_of_one_ensemble(ieta,event_ID_ens,temp_Mpt,temp_varMpt);
      sumMean[ieta]   += temp_Mpt ;
      sumMeanSq[ieta] += temp_Mpt * temp_Mpt ;
      sumVar[ieta]    += temp_varMpt ;
      sumVarSq[ieta]  += temp_varMpt * temp_varMpt ;
    } // loop over eta
  } // ievents loop

  for(int ieta=0; ieta<Neta; ieta++){
   sumMean[ieta]   /= Nevents ; 
   sumMeanSq[ieta] /= Nevents ; 
   sumVar[ieta]    /= Nevents ; 
   sumVarSq[ieta]  /= Nevents ; 
  }

  std::ofstream mFile;
  std::stringstream output_filename;
  // write to file
  output_filename.str("");
  output_filename << "results/meanpt_and_variance_of_meanpt";
  output_filename << "_pt_";
  output_filename << ptmin << "_" << ptmax << "_with" ;
  if(yflag==1){
   output_filename << "_y" ;
  }
  else{
   output_filename << "_eta" ;
  }
  output_filename << "_" << part ;
  output_filename << ".dat";
  mFile.open(output_filename.str().c_str(), std::ios::out );
  mFile << "eta   <pt>   error   Var(<pt>)  error" << std::endl ;
  for(int ieta=0; ieta<Neta; ieta++){
   mFile << eta_bin_centers[ieta] << "   " << sumMean[ieta] << "   " 
   << sqrt(sumMeanSq[ieta] - sumMean[ieta]*sumMean[ieta]) << "   "
   << sumVar[ieta] << "  "
   << sqrt(sumVarSq[ieta] - sumVar[ieta]*sumVar[ieta]) << std::endl ; 
  }
  mFile.close();

}




void mpt_decorr::calculate_covariance_of_mean_pt_of_one_ensemble(int ieta1, int ieta2, 
std::vector<int> event_ID_ens, double& cov, double& RMpt){
  double sum_mpt_1 = 0. ; 
  double sum_mpt_2 = 0. ; 
  double sum_mpt_1_sq = 0. ; 
  double sum_mpt_2_sq = 0. ; 
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ;
    sum_mpt_1 += get_IptI(eventID,ieta1); 
    sum_mpt_2 += get_IptI(eventID,ieta2); 
    sum_mpt_1_sq += get_IptI(eventID,ieta1) * get_IptI(eventID,ieta1); 
    sum_mpt_2_sq += get_IptI(eventID,ieta2) * get_IptI(eventID,ieta2); 

  }
  sum_mpt_1    /= event_ID_ens.size() ; 
  sum_mpt_2    /= event_ID_ens.size() ; 
  sum_mpt_1_sq /= event_ID_ens.size() ; 
  sum_mpt_2_sq /= event_ID_ens.size() ; 

  double sum_cov=0. ;
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ;
    sum_cov +=  (get_IptI(eventID,ieta1) - sum_mpt_1) * (get_IptI(eventID,ieta2) - sum_mpt_2) ; 
  }
  cov = sum_cov / event_ID_ens.size() ; 
  if( (sum_mpt_1_sq - sum_mpt_1 * sum_mpt_1 ) < 0 || (sum_mpt_2_sq - sum_mpt_2 * sum_mpt_2 ) < 0 ){
    std::cerr << "Cov_<pT> / R_<pT> calculation error !!!" << std::endl ; exit(-1); 
  }
  else{ 
    RMpt = cov / ( sqrt(sum_mpt_1_sq - sum_mpt_1 * sum_mpt_1 ) * sqrt(sum_mpt_2_sq - sum_mpt_2 * sum_mpt_2 ) );
  }
}



void mpt_decorr::calculate_r_mean_pt_of_one_ensemble(int ieta1, int ieta2, 
std::vector<int> event_ID_ens, double& rMpt){
 double cov_num ;
 double cov_den ; 
 double dummy ;
 // the bin at a given rapidity window -eta in the other side of pseudorapidity is (Neta-1-ieta2)
 calculate_covariance_of_mean_pt_of_one_ensemble(ieta1, (Neta-1-ieta2), event_ID_ens, cov_num, dummy);
 calculate_covariance_of_mean_pt_of_one_ensemble(ieta1, ieta2, event_ID_ens, cov_den, dummy);
 rMpt = cov_num / cov_den ; 
}



void mpt_decorr::write_covariance_of_meanpt(){
  std::ofstream mFile;
  std::stringstream output_filename;
  std::vector<int> event_ID_ens;
  double sumMean1[Neta];
  double sumMeanSq1[Neta];
  double sumMean2[Neta];
  double sumMeanSq2[Neta];
  double sumMean3[Neta];
  double sumMeanSq3[Neta];
  double temp1, temp2, temp3 ; 

  for(int ieta=0; ieta<Neta; ieta++){
   sumMean1[ieta]   = 0. ; 
   sumMeanSq1[ieta] = 0. ; 
   sumMean2[ieta]   = 0. ; 
   sumMeanSq2[ieta] = 0. ; 
   sumMean3[ieta]   = 0. ; 
   sumMeanSq3[ieta] = 0. ; 
  }

  for(int ieta1=0; ieta1<Neta; ieta1++){
    for(int ieta2=0; ieta2<Neta; ieta2++){
      sumMean1[ieta2]   = 0. ; 
      sumMeanSq1[ieta2] = 0. ; 
      sumMean2[ieta2]   = 0. ; 
      sumMeanSq2[ieta2] = 0. ; 
      sumMean3[ieta2]   = 0. ; 
      sumMeanSq3[ieta2] = 0. ; 
      for(int ievents=0; ievents<Nevents; ievents++){
        event_ID_ens = get_an_event_ensemble();
        calculate_covariance_of_mean_pt_of_one_ensemble(ieta1, ieta2, event_ID_ens, temp1, temp2);
        calculate_r_mean_pt_of_one_ensemble(ieta1, ieta2, event_ID_ens, temp3);
        sumMean1[ieta2]   += temp1 ;
        sumMeanSq1[ieta2] += temp1 * temp1 ;
        sumMean2[ieta2]   += temp2 ;
        sumMeanSq2[ieta2] += temp2 * temp2 ;
        sumMean3[ieta2]   += temp3 ;
        sumMeanSq3[ieta2] += temp3 * temp3 ;
      } // ievents loop

      sumMean1[ieta2]   /= Nevents ; 
      sumMeanSq1[ieta2] /= Nevents ; 
      sumMean2[ieta2]   /= Nevents ; 
      sumMeanSq2[ieta2] /= Nevents ; 
      sumMean3[ieta2]   /= Nevents ; 
      sumMeanSq3[ieta2] /= Nevents ;  
     
    } // loop over eta2


      // write to file
      output_filename.str("");
      output_filename << "results/Covariance_of_meanpt";
      output_filename << "_pt_";
      output_filename << ptmin << "_" << ptmax << "_with" ;
      if(yflag==1){
       output_filename << "_y" ;
      }
      else{
       output_filename << "_eta" ;
      }
      output_filename << "_" << part ;
      output_filename << "_etaref_" << eta_bin_centers[ieta1] ;
      output_filename << ".dat";
      mFile.open(output_filename.str().c_str(), std::ios::out );
      mFile << "eta1  eta2   cov(eta1,eta2)   error   Rpt(eta1,eta2)  error  rpt(eta1,eta2)" << std::endl ;
      for(int ieta2=0; ieta2<Neta; ieta2++){
        mFile << eta_bin_centers[ieta1] << "   " << eta_bin_centers[ieta2] << "   " <<  sumMean1[ieta2] << "  "
         << sqrt(sumMeanSq1[ieta2] - sumMean1[ieta2]*sumMean1[ieta2]) << "   " << sumMean2[ieta2] << "  "
         << sqrt(sumMeanSq2[ieta2] - sumMean2[ieta2]*sumMean2[ieta2]) << "   " << sumMean3[ieta2] << "  "
         << sqrt(sumMeanSq3[ieta2] - sumMean3[ieta2]*sumMean3[ieta2]) << "   " << std::endl ; 
      }
      mFile.close();

  } // loop over eta1

}



