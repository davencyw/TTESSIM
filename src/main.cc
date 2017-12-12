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

  auto start = std::chrono::system_clock::now();
  auto starttime = std::chrono::system_clock::to_time_t(start);

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
            << "started at: " << std::ctime(&starttime) << "\n";

  omp_set_num_threads(simenv._nThreads);

  // start main program
  Tstorageunit tsunit(simenv);
  tsunit.run(simenv._numcycles);

  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
  std::cout << "elapsed: " << static_cast<float>(elapsed.count()) / 60.0f
            << " minutes\n\n\n";

  return 0;
}