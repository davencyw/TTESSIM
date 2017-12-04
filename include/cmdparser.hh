/*__DECLARATION__
 *
 *
 *      TTESSIM
 *      simulation of operating a thermocline thermal energy storage unit at
 * high temperatures
 *
 *
 *
 *      author: david schmidig [     david@davencyw.net   ]
 *      ETH Zurich             [ davschmi@student.ethz.ch ]
 *      DAVENCYW CODE          [        davencyw.net      ]
 */
#ifndef __CMDPARSER_HH__
#define __CMDPARSER_HH__

#include <cstdlib>
#include <ctime>
#include <string>
#include "boost/program_options.hpp"

#include "global.hh"
#include "simenv.hh"

void cmdpars(int argc, char const* argv[], SimEnv& simenv) {
  // set default simenv
  simenv._nThreads = __P_DEF_NTHREADS;
  simenv._scheduling = __P_DEF_SCHEDULE;
  simenv._cuda = __P_DEF_CUDA;
  simenv._outfolder = __P_DEF_OUTFOLDER;

  // BOOST PRORGAM OPTIONS
  namespace po = boost::program_options;
  po::options_description desc("Parameters");
  desc.add_options()("help", "Print help messages")(
      ",N", po::value<int>(&(simenv._numcells))->required(), "number of cells")(
      ",h", po::value<precision_t>(&(simenv._storage_height))->required(),
      "storage height")(
      ",d", po::value<precision_t>(&(simenv._storage_diameter))->required(),
      "storage diameter")(
      ",t", po::value<precision_t>(&(simenv._fluid_initemp))->required(),
      "initial temperature of fluid")(
      "tds0", po::value<precision_t>(&(simenv._timedurstate0))->required(),
      "time duration for state 0")(
      "tds1", po::value<precision_t>(&(simenv._timedurstate1))->required(),
      "time duration for state 1")(
      "tds2", po::value<precision_t>(&(simenv._timedurstate2))->required(),
      "time duration for state 2")(
      "tds3", po::value<precision_t>(&(simenv._timedurstate3))->required(),
      "time duration for state 3")(
      "kf", po::value<precision_t>(&(simenv._kf))->required(), "kf")(
      "ks", po::value<precision_t>(&(simenv._ks))->required(), "ks")(
      "hf", po::value<precision_t>(&(simenv._hf))->required(), "hf")(
      "hs", po::value<precision_t>(&(simenv._hs))->required(), "hs")(
      "rhof", po::value<precision_t>(&(simenv._rhof))->required(), "rhof")(
      "rhos", po::value<precision_t>(&(simenv._rhos))->required(), "rhos")(
      "cf", po::value<precision_t>(&(simenv._cf))->required(), "cf")(
      "cs", po::value<precision_t>(&(simenv._cs))->required(), "cs")(
      "epsilon", po::value<precision_t>(&(simenv._epsilon))->required(),
      "epsilon")("uf", po::value<precision_t>(&(simenv._uf))->required(), "uf")(
      "dt", po::value<precision_t>(&(simenv._deltat))->required(),
      "timestepsize")("numcycles,c",
                      po::value<int>(&(simenv._numcycles))->required(),
                      "number of cycles")(
      "tstbc,b", po::value<int>(&(simenv._tsteppercycle))->required(),
      "timesteps between cycless")("folder",
                                   po::value<std::string>(&(simenv._outfolder)),
                                   "output folder [optional]")(
      "nthreads,n", po::value<int>(&(simenv._nThreads)),
      "number of threads [optional]")(
      "schedule,s", po::value<int>(&(simenv._scheduling)),
      "omp scheduling [optional]")("cuda", "enable cuda support [optional]");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);  // can throw

    /** --help option
     */
    if (vm.count("help")) {
      std::cout << __P_NAME << std::endl << desc << std::endl;
      exit(0);
    }

    if (vm.count("cuda")) simenv._cuda = true;

    po::notify(vm);  // throws on error, so do after help in case
                     // there are any problems
  } catch (po::error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }

  // sanity check input
  if (simenv._nThreads > 256) {
    std::cerr << "MISCONFIG: nThreads = " << simenv._nThreads << std::endl;
    exit(1);
  }

  if (simenv._numcells > 10e10) {
    std::cerr << "MISCONFIG: numcells = " << simenv._numcells << std::endl;
    exit(1);
  }

  // TODO(dave): sanitize (folder) input
  // runhash
  simenv._runhash = std::time(0);
  simenv._fullpath_fluid =
      simenv._outfolder + "r_" + std::to_string(simenv._runhash) + "_f.csv";
  simenv._fullpath_solid =
      simenv._outfolder + "r_" + std::to_string(simenv._runhash) + "_s.csv";
}

#endif  //__CMDPARSER_HH__