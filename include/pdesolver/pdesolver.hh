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
#ifndef __PDESOLVER_HH__
#define __PDESOLVER_HH__

#include <functional>
#include "global.hh"
#include "simenv.hh"

class Pdesolver {
 public:
  Pdesolver() = default;
  Pdesolver(const SimEnv* simenv) : _simenv(simenv) {
    _dt = _simenv->_deltat;
    _dx =
        _simenv->_storage_height / static_cast<precision_t>(_simenv->_numcells);
    _dx2 = _dx * _dx;
    _idx = 1.0 / _dx;
    _idx2 = 1.0 / _dx2;
    _kf = _simenv->_kf;
    _ks = _simenv->_ks;
    _epsilon = _simenv->_epsilon;
    _rhof = _simenv->_rhof;
    _rhos = _simenv->_rhos;
    _cf = _simenv->_cf;
    _cs = _simenv->_cs;
    _alphaf = _kf / (_epsilon * _rhof * _cf);
    _alphas = _ks / ((1 - _epsilon) * _rhos * _cs);
    _alphafidx2 = _alphaf * _idx2;
    _alphasidx2 = _alphas * _idx2;
    _alphafidx2dt = _alphafidx2 * _dt;
    _alphasidx2dt = _alphasidx2 * _dt;
    _uf = _simenv->_uf;
    _numcells = _simenv->_numcells;
    _hf = _simenv->_hf;
    _hs = _simenv->_hs;
  };

  void solvediffusion(array_t* temperature, array_t* temperature_old,
                      precision_t diffusionnumber);
  void solveadvection(array_t* temperature, array_t* temperature_old,
                      precision_t cflnumber, precision_t boundary_temperature);
  void solvefluid(array_t** temperature, array_t** temperature_old,
                  precision_t cflnumber, precision_t diffusionnumber,
                  precision_t boundary_temperature);
  void solvesolid(array_t** temperature, array_t** temperature_old,
                  precision_t diffusionnumber);
  void solvecoupling(array_t* temperature_solid, array_t* temperature_fluid);

#ifdef TESTING
  void testing();
#endif

 private:
#ifdef TESTING
  void verify(precision_t* errorf, precision_t* errors, precision_t* iterf,
              precision_t* iters);
#endif

  const SimEnv* _simenv;

  precision_t _dt;
  precision_t _dx;
  precision_t _dx2;
  precision_t _idx;
  precision_t _idx2;
  precision_t _kf;
  precision_t _ks;
  precision_t _epsilon;
  precision_t _rhof;
  precision_t _rhos;
  precision_t _cf;
  precision_t _cs;
  precision_t _alphaf;
  precision_t _alphas;
  precision_t _alphafidx2;
  precision_t _alphasidx2;
  precision_t _alphafidx2dt;
  precision_t _alphasidx2dt;
  precision_t _uf;
  precision_t _hf;
  precision_t _hs;
  int _numcells;

// DEBUG
#ifdef TESTING
  precision_t _k;
  static constexpr precision_t k_tol = 1e-10;
  static constexpr int k_maxiterations = 500000;
  array_t _source_fluid;
  array_t _source_solid;
#endif
};

#endif  //__PDESOLVER_HH__