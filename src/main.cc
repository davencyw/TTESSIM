#include <climits>

#include <omp.h>
#include <iostream>
#include <string>

#include "simenv.hh"
#include "cmdparser.hh"
#include "tstorageunit/tstorageunit.hh"

int main(int argc, char const *argv[]) {

  // parse command line arguments
  SimEnv simenv;
  cmdpars(argc, argv, simenv);

  // create terminal output
  std::cout << "\n\n\n\n\n"
            <<"_______ _______ _______ _______ _______ _____ _______\n"
            <<"   |       |    |______ |______ |______   |   |  |  |\n"
            <<"   |       |    |______ ______| ______| __|__ |  |  |\n"
            << "\n\n\n\n"
            << "author: david schmidig [david@davencyw.net]\n"
            << "        davencyw code  [davencyw.net]\n"
            << "        ETH Zurich\n\n\n";

  // start main program
  Tstorageunit tsunit(simenv);
  tsunit.writetocsv();

  return 0;
}