#include "inparams.h"

read_parameters::read_parameters(){
}

read_parameters::~read_parameters(){
}

// This functions reads the input file and sets the input parameters in a structure for global use in the code. 
void read_parameters::read_parameters_from_file(input_paramters &iparam, string input_file_name)
{
  string param_name; char param_value[50];
  
  istringstream* iss;
  char   buff[200];
  
  fstream File0;
  File0.open(input_file_name,ios::in);
  if(!File0){
     cout<<"No input parameter file found.\nexit(1)";
     exit(1);
  }
  
  int number = 0;
  while(!File0.eof())
    {
      File0.getline(buff,200);
      
      if (!(*buff) || (*buff == '#')){
	number ++; 
	continue;
      }
      
      iss = new istringstream(buff);
      *iss >> param_name >> param_value ;
      delete iss ; 
      
      if(param_name == "nEvents"  ){
	iparam.nEvents = atof(param_value);
      }
      
      if(param_name == "input_read_mode"  ){
	iparam.input_read_mode = atof(param_value);
      }

      if(param_name == "include_decay"  ){
	iparam.include_weak_decay = atof(param_value);
      }

      if(param_name == "reconstruct_phi_flag"  ){
	iparam.reconstruct_phi_flag = atof(param_value);
      }

      if(param_name == "reconstruct_kstar0_flag"  ){
	iparam.reconstruct_kstar0_flag = atof(param_value);
      }

      if(param_name == "reconstruct_kstar0_bar_flag"  ){
	iparam.reconstruct_kstar0_bar_flag = atof(param_value);
      }

      if(param_name == "reconst_phi_mass"  ){
	iparam.reconst_phi_mass = atof(param_value);
      }

      if(param_name == "reconst_phi_width"  ){
	iparam.reconst_phi_width = atof(param_value);
      }
   
      if(param_name == "reconst_kstar0_mass"  ){
	iparam.reconst_kstar0_mass = atof(param_value);
      }

      if(param_name == "reconst_kstar0_width"  ){
	iparam.reconst_kstar0_width = atof(param_value);
      }

      if(param_name == "reconst_kstar0_bar_mass"  ){
	iparam.reconst_kstar0_bar_mass = atof(param_value);
      }

      if(param_name == "reconst_kstar0_bar_width"  ){
	iparam.reconst_kstar0_bar_width = atof(param_value);
      }

      if(param_name == "pdg_type"  ){
	iparam.pdg_type = atof(param_value);
      }


    number++; 
    
    } 
  File0.close();
}











