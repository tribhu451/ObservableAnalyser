void convert_to_root_tree(){

  std::string dummy ;
  int nch ;
  int index, pid ;
  double  x, y, e, px, py, pz,mass, t, z;
  int temp_pid, temp_iso ;
  int Event, Multiplicity ; 
  int TotalEvents = 1 ; 

  istringstream* iss;
  char   buff[400];

  const int Nmult = 5000 ;  
  double T[Nmult] ; 
  double X[Nmult] ; 
  double Y[Nmult] ; 
  double Z[Nmult] ; 
  double E[Nmult] ; 
  double PX[Nmult] ; 
  double PY[Nmult] ; 
  double PZ[Nmult] ; 
  double MASS[Nmult] ; 
  double PID[Nmult] ; 

  for(int ii = 0; ii < Nmult; ii++){
	  T[ii] = 0 ; 
  }



  // =================================  //
  TFile* File2= new TFile ("particle_list.root","recreate");
  TTree* tree1= new TTree("TreePart","Tree for UrQMD generated particles");
  tree1->Branch("Event",&Event,"Event/I");
  tree1->Branch("Multiplicity",&Multiplicity,"Multiplicity/I");
  tree1->Branch("PID",&PID,"PID[Nmult]/I");
  tree1->Branch("t",&T,"T[Nmult]/F");
  tree1->Branch("x",&X,"X[Nmult]/D");
  tree1->Branch("y",&Y,"Y[Nmult]/D");
  tree1->Branch("x",&Z,"Z[Nmult]/D");
  tree1->Branch("e",&E,"E[Nmult]/D");
  tree1->Branch("px",&PX,"PX[Nmult]/D");
  tree1->Branch("py",&PY, "PY[Nmult]/D");
  tree1->Branch("pz",&PZ, "PZ[Nmult]/D");
  tree1->Branch("mass",&MASS, "MASS[Nmult]/D");
  // =================================  //
  tree1->SetAutoSave(10737418240);
  File2->SetCompressionLevel(9);

  std::cout << "reading particle_list.dat from urqmd ..."
              << std::endl ;

  std::ifstream file;
  file.open("particle_list.dat");
  if(!file){
    std::cout << "file not found."
                << std::endl;
    exit(1);
  }

  for(int ii=0; ii<TotalEvents; ii++){
    
    Event = ii ; 

    for(int ij=0; ij<17; ij++){
      file.getline(buff,200) ;
    }

    file.getline(buff,100);
    iss = new istringstream(buff);
    *iss >> nch >> dummy ;
    delete iss;
    std::cout << ii << "\t" << nch << std::endl ;
    Multiplicity = nch ;
    file.getline(buff,200) ; // read an useless line

    for(int i=0; i<nch ; i++){
      file.getline(buff,320);
      iss = new istringstream(buff);
      *iss       >> dummy >> dummy >> dummy >> dummy
                 >> dummy >> dummy >> dummy >> dummy
                 >> mass >> pid >> temp_iso
                 >> dummy >> dummy >> dummy >> dummy
                 >> t >> x >> y >> z
                 >> e >> px >> py >> pz ;

    MASS[i] = mass;
    PID[i] = pid;
    T[i] = t;
    X[i] = x;
    Y[i] = y;
    Z[i] = z;
    E[i] = e;
    PX[i] = px;
    PY[i] = py;
    PZ[i] = pz ;
   
    tree1->Fill() ;  
    delete iss;
    
    } // particle loop

  }

tree1->Write();
File2->Close();



}
