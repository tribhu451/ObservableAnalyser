#include "random.h"



random_gen::random_gen() : rand_uniform_dist(0.,1.) { 
  seed = ran_dev();
  //seed = 12345 ;  // avoid the randomness by turning on this one (only for check purpose) ... :)
  ran_generator = std::unique_ptr<std::mt19937>(new std::mt19937(seed));
}


random_gen::~random_gen(){
}


double random_gen::rand_uniform(){
  return(rand_uniform_dist(*ran_generator));
}


int random_gen::get_seed() {
  return(seed);
}






