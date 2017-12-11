#include <climits>

#include <omp.h>
#include <iostream>
#include <string>

#include "cmdparser.hh"
#include "pdesolver/pdesolver.hh"
#include "simenv.hh"
#include "tstorageunit/tstorageunit.hh"

int main(int argc, char const *argv[]) {
  // parse command line arguments
  SimEnv simenv;
  cmdpars(argc, argv, simenv);

  // terminal output
  std::cout << "\n\n\n\n\n"
            << "_______ _______ _______ _______ _______     _______\n"
            << "   |       |    |______ |______ |______  |  |  |  |\n"
            << "   |       |    |______ ______| ______|  |  |  |  |\n"
            << "\n\n\n\n"
            << "author: david schmidig [david@davencyw.net]\n"
            << "        davencyw code  [davencyw.net]\n"
            << "        ETH Zurich\n\n"
            << "\033[1;31m"
            << "T E S T I N G T E S T I N G T E S T I N T E S T I N G\n"
            << "T E S T I N G T E S T I N G T E S T I N T E S T I N G\n"
            << "T E S T I N G T E S T I N G T E S T I N T E S T I N G"
            << "\033[0m\n\n\n";

  // start testing

  const unsigned int incremental(10);
  const unsigned int start(10);
  array_t heights(4);
  heights << 0.1, 1.0, 20.0, 100.0;

  // for different peclet numbers
  for (int j = 0; j < 4; ++j) {
    // for different meshsizes
    for (int i = 1; i < 20; ++i) {
      simenv._numcells = start + i * incremental;
      simenv._storage_height = heights(j);
      Pdesolver pdesolver(&simenv);
      pdesolver.testing();
    }
  }

  return 0;
}