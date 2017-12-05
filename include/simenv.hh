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
#ifndef __SIMENV_HH__
#define __SIMENV_HH__

#include <fstream>
#include <string>

#include "global.hh"

struct SimEnv {
  SimEnv() {
    _fs_fluid = new std::ofstream;
    _fs_solid = new std::ofstream;
  }

  // physical environment
  //____________________

  // number of cells
  int _numcells;
  // height of storage unit (cylindrical)
  precision_t _storage_height;
  // diameter of storage unit (cylindrical)
  precision_t _storage_diameter;
  // initial temperature of fluid
  precision_t _fluid_initemp;
  // timestepsize
  precision_t _deltat;

  // times for states
  precision_t _timedurstate0;
  precision_t _timedurstate1;
  precision_t _timedurstate2;
  precision_t _timedurstate3;

  // fluid properties
  precision_t _kf;
  precision_t _ks;
  precision_t _rhof;
  precision_t _rhos;
  precision_t _cf;
  precision_t _cs;
  precision_t _epsilon;
  precision_t _uf;
  precision_t _hf;
  precision_t _hs;

  // computing environment
  //____________________

  // number of threads
  int _nThreads;
  // omp scheduling
  int _scheduling;
  // cuda enabled
  bool _cuda;
  // number of cycles
  int _numcycles;
  // tsteps per cycle;
  int _tsteppercycle;
  // every _ops writetocsv
  int _ops;

  // data environment
  //____________________

  std::string _outfolder;
  int _runhash;
  std::string _fullpath_solid;
  std::string _fullpath_fluid;
  std::ofstream* _fs_fluid;
  std::ofstream* _fs_solid;
};

#endif  //__SIMENV_HH__