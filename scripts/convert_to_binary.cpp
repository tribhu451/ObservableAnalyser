#include <iostream>
#include <fstream>
#include <sstream>

using std::istringstream ;

int main(){

  int EventsInEachFile = 1000 ;  
  int start_ifile = 0 ; 
  int end_ifile = 1 ; 

  std::string dummy ;
  int nch ;
  int index, pid ;
  double  x, y, e, px, py, pz,mass, t, z;
  int temp_pid, temp_iso ;
  int TotalEvents = 0 ; 
  istringstream* iss;
  char   buff[400];

  std::cout << "reading particle_list.dat from urqmd ..."
              << std::endl ;

  std::ofstream outbin;
  std::stringstream output_filename1;
  output_filename1 <<"particle_list.bin";
  outbin.open(output_filename1.str().c_str(), std::ios::out | std::ios::binary);

  std::ifstream file;
  
  for(int ifile=start_ifile; ifile<end_ifile; ifile++ ){

  std::stringstream input_filename;
  input_filename.str(std::string());
  input_filename << "FILE"<<ifile <<"/particle_list.dat";
  file.open(input_filename.str().c_str(), std::ios::in);
  std::cout << "reading file : " << input_filename.str().c_str() << "... " << std::endl ; 
  if(!file){
    std::cout << "file not found."
                << std::endl;
    exit(1);
  }


  for(int ii=0; ii<EventsInEachFile; ii++){
    

    for(int ij=0; ij<17; ij++){
      file.getline(buff,200) ;
    }

    file.getline(buff,100);
    iss = new istringstream(buff);
    *iss >> nch >> dummy ;
    delete iss;
    std::cout << TotalEvents << "  " << ii << "  " << nch << std::endl ;
    file.getline(buff,200) ; // read an useless line
    outbin.write(reinterpret_cast<char *>(&nch),
                         sizeof(int));
    for(int i=0; i<nch ; i++){
      file.getline(buff,320);
      iss = new istringstream(buff);
      *iss       >> dummy >> dummy >> dummy >> dummy
                 >> dummy >> dummy >> dummy >> dummy
                 >> mass >> pid >> temp_iso
                 >> dummy >> dummy >> dummy >> dummy
                 >> t >> x >> y >> z
                 >> e >> px >> py >> pz ;
      float particle_array[] = {          static_cast<float>(mass),
                                          static_cast<float>(pid),
                                          static_cast<float>(temp_iso),
                                          static_cast<float>(t),
                                          static_cast<float>(x),
                                          static_cast<float>(y),
                                          static_cast<float>(z),
                                          static_cast<float>(e),
                                          static_cast<float>(px),
                                          static_cast<float>(py),
                                          static_cast<float>(pz)};
                for (int ii = 0; ii < 11 ; ii++) {
                    outbin.write(
                        reinterpret_cast<char *>(&(particle_array[ii])),
                        sizeof(float));
                }
    delete iss;
    
    } // particle loop
   TotalEvents += 1 ; 
  } // event loop
  
  file.close() ; 
  } // ifile loop
  
  std::cout << "Total Events : " << TotalEvents << std::endl ; 
  outbin.close() ; 

}
