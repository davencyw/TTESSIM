#include <climits>

#include <omp.h>
#include <iostream>
#include <string>

#include "simenv.hh"
#include "cmdparser.hh"
#include "pdesolver/pdesolver.hh"
#include "tstorageunit/tstorageunit.hh"

int main(int argc, char const *argv[]) {

  // parse command line arguments
  SimEnv simenv;
  cmdpars(argc, argv, simenv);

  // create terminal output
  std::cout << "\n\n\n\n\n"
            <<"_______ _______ _______ _______ _______     _______\n"
            <<"   |       |    |______ |______ |______  |  |  |  |\n"
            <<"   |       |    |______ ______| ______|  |  |  |  |\n"
            << "\n\n\n\n"
            << "author: david schmidig [david@davencyw.net]\n"
            << "        davencyw code  [davencyw.net]\n"
            << "        ETH Zurich\n\n"
            <<"\033[1;31mT E S T I N G T E S T I N G T E S T I N T E S T I N G\033[0m\n\n\n";


  // start testing
  Pdesolver pdesolver(&simenv);
  pdesolver.testing();

  return 0;
}