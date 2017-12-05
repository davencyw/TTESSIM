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
  std::cout
      << "\n\n\n\n\n"
      << "_______ _______ _______ _______ _______     _______\n"
      << "   |       |    |______ |______ |______  |  |  |  |\n"
      << "   |       |    |______ ______| ______|  |  |  |  |\n"
      << "\n\n\n\n"
      << "author: david schmidig [david@davencyw.net]\n"
      << "        davencyw code  [davencyw.net]\n"
      << "        ETH Zurich\n\n"
      << "\033[1;31mT E S T I N G T E S T I N G T E S T I N T E S T I N G\n"
      << "T E S T I N G T E S T I N G T E S T I N T E S T I N G\n"
      << "T E S T I N G T E S T I N G T E S T I N T E S T I N G\033[0m\n\n\n";

  // start testing

  unsigned int incremental(100);

  for (int i = 1; i < 20; ++i) {
    simenv._numcells = i * incremental;
    Pdesolver pdesolver(&simenv);
    pdesolver.testing();
  }

  // TODO(dave): OVS -> change simenv and start testing again, file descriptor
  // is still the same

  return 0;
}