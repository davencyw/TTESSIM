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

#include <eigen3/Eigen/Dense>

class Tstorageunit {
 public:
  Tstorageunit(const SimEnv& simenv) : _simenv(simenv) {
    _fluid_temperature =
        array_t::Constant(_simenv._numcells, _simenv._fluid_initemp);
    _solid_temperature =
        array_t::Constant(_simenv._numcells, _simenv._fluid_initemp);
    _fluid_temperature_o =
        array_t::Constant(_simenv._numcells, _simenv._fluid_initemp);
    _solid_temperature_o =
        array_t::Constant(_simenv._numcells, _simenv._fluid_initemp);

    _fluid_temperature_ptr = &_fluid_temperature;
    _solid_temperature_ptr = &_solid_temperature;
    _fluid_temperature_o_ptr = &_fluid_temperature_o;
    _solid_temperature_o_ptr = &_solid_temperature_o;

    _pdesolver = Pdesolver(&simenv);
    precision_t inittemp(_simenv._fluid_initemp);

    _dx = _simenv._storage_height / static_cast<precision_t>(_simenv._numcells);
    _uf = _simenv._uf;
    _total_time = 0.0;

    precision_t alphaf =
        _simenv._kf / (_simenv._epsilon * _simenv._rhof * _simenv._cf);
    precision_t alphas =
        _simenv._ks / ((1 - _simenv._epsilon) * _simenv._rhos * _simenv._cs);

    _solid_diffusionnumber = _simenv._deltat / (_dx * _dx) * alphas;
    _fluid_diffusionnumber = _simenv._deltat / (_dx * _dx) * alphaf;
    updatecfl(_uf);
    // TODO(dave): Check!
    _boundary_temperature = _simenv._fluid_initemp;

    // initialize files with metadata

    _simenv._fs_fluid->open(_simenv._fullpath_fluid,
                            std::ofstream::out | std::ofstream::app);
    *_simenv._fs_fluid << _simenv._numcells << ";" << _simenv._deltat << ";"
                       << _simenv._ops << ";" << _simenv._storage_height
                       << "\n";

    _simenv._fs_solid->open(_simenv._fullpath_solid,
                            std::ofstream::out | std::ofstream::trunc);
    *_simenv._fs_solid << _simenv._numcells << ";" << _simenv._deltat << ";"
                       << _simenv._ops << ";" << _simenv._storage_height
                       << "\n";

    _time_per_cycle = _simenv._timedurstate0 + _simenv._timedurstate1 +
                      _simenv._timedurstate2 + _simenv._timedurstate3;

    _exergy_flux_array = array_t::Zero(4);
  };

  ~Tstorageunit() {
    // Close fs
    _simenv._fs_fluid->close();
    _simenv._fs_solid->close();
  }

  // step-by-step simulation
  bool simstep();
  bool simsteps(const int);
  bool simsteps(const int, const int);

  // full sim
  bool run(unsigned int cycles);

 private:
  // state
  // 0 = charging (+)
  // 1 = discharging (-)
  // 2 = idle between charging and discharging (+-)
  // 3 = idle between discharging and charging (-+)
  const int getstate();
  void updatecfl(precision_t uf);

  void computeefficiency();

  // output functions
  const bool writetocsv(array_t*, int, std::ofstream* stream);

  // simulation environment
  const SimEnv _simenv;

  // physical properties of system
  precision_t _total_time;
  array_t _fluid_temperature;
  array_t _solid_temperature;
  array_t _fluid_temperature_o;
  array_t _solid_temperature_o;

  array_t* _fluid_temperature_ptr;
  array_t* _solid_temperature_ptr;
  array_t* _fluid_temperature_o_ptr;
  array_t* _solid_temperature_o_ptr;

  precision_t _solid_diffusionnumber;
  precision_t _fluid_diffusionnumber;
  precision_t _cflnumber;
  precision_t _boundary_temperature;
  precision_t _dx;
  precision_t _uf;

  precision_t _time_per_cycle;

  // efficiency
  precision_t _exergy_flux;
  array_t _exergy_flux_array;
  precision_t _capacity_factor;

  // pde solver
  Pdesolver _pdesolver;
};

#endif  //__TSTORAGEUNIT_HH__