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
  };

  void solvefluid(precision_t* fluid_temperature,
                  precision_t* fluid_temperature_o, precision_t boundary);
  void solvesolid(precision_t* solid_temperature,
                  precision_t* solid_temperature_o, precision_t boundary);

#ifdef TESTING
  void testing();
#endif

 private:
#ifdef TESTING
  // These functions solve the governing fluid and
  // solid equations for a MMS given the slack term as
  // default or passed lambda expression.
  bool verifyfluid(precision_t* error);
  bool verfiysolid(precision_t* error);
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
  int _numcells;

// DEBUG
#ifdef TESTING
  int _n;
  precision_t _k;
  static constexpr precision_t _tol = 1e-8;
  static constexpr int _maxiterations = 10000;
  precision_t* _source_fluid;
  precision_t* _source_solid;
#endif
};

#endif  //__PDESOLVER_HH__