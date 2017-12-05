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
#ifndef __TSTORAGEUNIT_HH__
#define __TSTORAGEUNIT_HH__

#include "const.hh"
#include "global.hh"
#include "pdesolver/pdesolver.hh"
#include "simenv.hh"

#include <immintrin.h>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>

class Tstorageunit {
 public:
  Tstorageunit(const SimEnv& simenv) : _simenv(simenv) {
    _state = 0;
    _fluid_temperature =
        (precision_t*)_mm_malloc(sizeof(precision_t) * _simenv._numcells, 32);
    _solid_temperature =
        (precision_t*)_mm_malloc(sizeof(precision_t) * _simenv._numcells, 32);
    _fluid_temperature_o =
        (precision_t*)_mm_malloc(sizeof(precision_t) * _simenv._numcells, 32);
    _solid_temperature_o =
        (precision_t*)_mm_malloc(sizeof(precision_t) * _simenv._numcells, 32);

    _pdesolver = Pdesolver(&simenv);
    precision_t inittemp(_simenv._fluid_initemp);

    for (int i = 0; i < _simenv._numcells; ++i) {
      _fluid_temperature[i] = _simenv._fluid_initemp;
      _solid_temperature[i] = _simenv._fluid_initemp;
      _fluid_temperature_o[i] = _simenv._fluid_initemp;
      _solid_temperature_o[i] = _simenv._fluid_initemp;
      _solid_temperature[i] = _simenv._fluid_initemp;
    }

    // TODO(dave): add outputstep such that physical timesteps are known
    // initialize files with metadata

    _simenv._fs_fluid->open(_simenv._fullpath_fluid,
                            std::ofstream::out | std::ofstream::app);
    *_simenv._fs_fluid << _simenv._numcells << ";" << _simenv._deltat << ";"
                       << "OUTPUTSTEP MISSING;" << _simenv._storage_height
                       << "\n";

    _simenv._fs_solid->open(_simenv._fullpath_solid,
                            std::ofstream::out | std::ofstream::trunc);
    *_simenv._fs_solid << _simenv._numcells << ";" << _simenv._deltat << ";"
                       << "OUTPUTSTEP MISSING;" << _simenv._storage_height
                       << "\n";
  };

  ~Tstorageunit() {
    _mm_free(_fluid_temperature);
    _mm_free(_fluid_temperature_o);
    _mm_free(_solid_temperature);
    _mm_free(_solid_temperature_o);

    // Close fs
    _simenv._fs_fluid->close();
    _simenv._fs_solid->close();
  }

  // step-by-step simulation
  bool simstep();
  bool simsteps(const int);
  bool simsteps(const int, const int);

 private:
  const int getstate();
  const bool isidle();

  // output functions
  const bool writetocsv(precision_t*, int, std::ofstream* stream);

  // simulation environment
  const SimEnv _simenv;

  // physical properties of system
  precision_t _total_time;
  precision_t* _fluid_temperature;
  precision_t* _solid_temperature;
  precision_t* _fluid_temperature_o;
  precision_t* _solid_temperature_o;
  // state
  // 0 = charging (+)
  // 1 = discharging (-)
  // 2 = idle between charging and discharging (+-)
  // 3 = idle between discharging and charging (-+)
  int _state;

  Pdesolver _pdesolver;
};

#endif  //__TSTORAGEUNIT_HH__