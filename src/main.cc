#include <chrono>
#include <climits>

#include <omp.h>
#include <iostream>
#include <string>

#include "cmdparser.hh"
#include "simenv.hh"
#include "tstorageunit/tstorageunit.hh"

int main(int argc, char const *argv[]) {
  // parse command line arguments
  SimEnv simenv;
  cmdpars(argc, argv, simenv);

  auto now = std::chrono::system_clock::now();
  auto starttime = std::chrono::system_clock::to_time_t(now);

  // create terminal output
  std::cout << "\n\n\n\n\n"
            << "_______ _______ _______ _______ _______     _______\n"
            << "   |       |    |______ |______ |______  |  |  |  |\n"
            << "   |       |    |______ ______| ______|  |  |  |  |\n"
            << "\n\n\n\n"
            << "author: david schmidig [david@davencyw.net]\n"
            << "        davencyw code  [davencyw.net]\n"
            << "        ETH Zurich\n\n\n"
            << "\n\n\n\n"
            << "started at: " << std::ctime(&starttime) << "\n\n\n\n";

  omp_set_num_threads(simenv._nThreads);

  // start main program
  Tstorageunit tsunit(simenv);
  tsunit.run(simenv._numcycles);

  return 0;
}