#include "pdesolver/pdesolver.hh"
#include "const.hh"
#include "global.hh"
#include "simenv.hh"

#include <stdio.h>
#include <cassert>
#include <cmath>

#include <algorithm>
#include <fstream>
#include <iostream>

// TODO(dave): Optimize arithmetic!! (reordering)
// TODO(dave): change switch to generic generator
void Pdesolver::solvefluid(precision_t** ft, precision_t** fto,
                           precision_t boundary, const unsigned int state) {
  precision_t* fluid_temperature = *ft;
  precision_t* fluid_temperature_o = *fto;

  assert(state < 4);

  const int N(_simenv->_numcells - 1);
  const int Nm1(_simenv->_numcells - 2);

  switch (0) {
    case 0: {
#ifdef __INTEL_COMPILER
#pragma ivdep
#elif __GNUC__
#pragma GCC ivdep
#endif
      for (int i = 1; i < _simenv->_numcells - 1; ++i) {
        const precision_t tfi = fluid_temperature[i];
        const precision_t tfim1 = fluid_temperature[i - 1];
        const precision_t tfip1 = fluid_temperature[i + 1];

        fluid_temperature_o[i] = tfi - _dt * _uf * _idx * (tfi - tfim1) +
                                 _alphafidx2dt * (tfim1 - 2.0 * tfi + tfip1);
      }
      // boundary cells
      fluid_temperature_o[0] =
          fluid_temperature[0] -
          _dt * _uf * _idx * (fluid_temperature[0] - boundary) +
          _alphafidx2dt *
              (fluid_temperature[1] - 2.0 * fluid_temperature[0] + boundary);
      fluid_temperature_o[N] =
          fluid_temperature[N] -
          _dt * _uf * _idx * (fluid_temperature[N] - fluid_temperature[Nm1]) +
          _alphafidx2dt * (fluid_temperature[Nm1] - fluid_temperature[N]);
    } break;
    case 1: {
#ifdef __INTEL_COMPILER
#pragma ivdep
#elif __GNUC__
#pragma GCC ivdep
#endif
      for (int i = 1; i < _simenv->_numcells - 1; ++i) {
        const precision_t tfi = fluid_temperature[i];
        const precision_t tfim1 = fluid_temperature[i - 1];
        const precision_t tfip1 = fluid_temperature[i + 1];

        fluid_temperature_o[i] =
            tfi + _alphafidx2dt * (tfim1 - 2.0 * tfi + tfip1);
      }
      // boundary cells
      fluid_temperature_o[0] =
          fluid_temperature[0] +
          _alphafidx2dt * (fluid_temperature[1] - fluid_temperature[0]);
      fluid_temperature_o[N] =
          fluid_temperature[N] -
          +_alphafidx2dt * (fluid_temperature[Nm1] - fluid_temperature[N]);
    } break;
    case 2: {
#ifdef __INTEL_COMPILER
#pragma ivdep
#elif __GNUC__
#pragma GCC ivdep
#endif
      for (int i = 1; i < _simenv->_numcells - 1; ++i) {
        const precision_t tfi = fluid_temperature[i];
        const precision_t tfim1 = fluid_temperature[i - 1];
        const precision_t tfip1 = fluid_temperature[i + 1];

        fluid_temperature_o[i] = tfi - _dt * _uf * _idx * (tfip1 - tfi) +
                                 _alphafidx2dt * (tfim1 - 2 * tfi + tfip1);
      }
      // boundary cells
      fluid_temperature_o[0] =
          fluid_temperature[0] -
          _dt * _uf * _idx * (fluid_temperature[1] - fluid_temperature[0]) +
          _alphafidx2dt * (fluid_temperature[1] - fluid_temperature[0]);
      fluid_temperature_o[N] =
          fluid_temperature[N] -
          _dt * _uf * _idx * (boundary - fluid_temperature[N]) +
          _alphafidx2dt *
              (boundary - 2.0 * fluid_temperature[N] + fluid_temperature[Nm1]);
    } break;
    case 3: {
#ifdef __INTEL_COMPILER
#pragma ivdep
#elif __GNUC__
#pragma GCC ivdep
#endif
      for (int i = 1; i < _simenv->_numcells - 1; ++i) {
        const precision_t tfi = fluid_temperature[i];
        const precision_t tfim1 = fluid_temperature[i - 1];
        const precision_t tfip1 = fluid_temperature[i + 1];

        fluid_temperature_o[i] =
            tfi + _alphafidx2dt * (tfim1 - 2 * tfi + tfip1);
      }
      // boundary cells
      fluid_temperature_o[0] =
          fluid_temperature[0] +
          _alphafidx2dt * (fluid_temperature[1] - fluid_temperature[0]);
      fluid_temperature_o[N] =
          fluid_temperature[N] -
          +_alphafidx2dt * (fluid_temperature[Nm1] - fluid_temperature[N]);
    } break;

    default: { std::cout << "E R R O R INVALID STATE\n"; }
  }

// mms
#ifdef TESTING
  assert(state == 0);
#ifdef __INTEL_COMPILER
#pragma ivdep
#elif __GNUC__
#pragma GCC ivdep
#endif
  for (int i = 0; i < _simenv->_numcells; ++i) {
    fluid_temperature_o[i] += _source_fluid[i];
  }
#endif

  std::swap(*ft, *fto);
  ;
}

void Pdesolver::solvesolid(precision_t** st, precision_t** sto) {
  // Loop over inner N-2 cells

  precision_t* solid_temperature = *st;
  precision_t* solid_temperature_o = *sto;

#ifdef __INTEL_COMPILER
#pragma ivdep
#elif __GNUC__
#pragma GCC ivdep
#endif
  for (int i = 1; i < _simenv->_numcells - 1; ++i) {
    const precision_t tsi = solid_temperature[i];
    const precision_t tsim1 = solid_temperature[i - 1];
    const precision_t tsip1 = solid_temperature[i + 1];
    solid_temperature_o[i] = tsi + _alphasidx2dt * (tsip1 - 2 * tsi + tsim1);
  }

  const int N(_simenv->_numcells - 1);
  const int Nm1(_simenv->_numcells - 2);
  // Boundary cells
  solid_temperature_o[0] =
      solid_temperature[0] +
      _alphasidx2dt * (solid_temperature[1] - solid_temperature[0]);
  solid_temperature_o[N] =
      solid_temperature[N] +
      _alphasidx2dt * (solid_temperature[Nm1] - solid_temperature[N]);

// mms
#ifdef TESTING
#ifdef __INTEL_COMPILER
#pragma ivdep
#elif __GNUC__
#pragma GCC ivdep
#endif
  for (int i = 0; i < _simenv->_numcells; ++i) {
    solid_temperature_o[i] += _source_solid[i];
  }
#endif

  std::swap(*st, *sto);
}

void Pdesolver::solvecoupling(precision_t* tf, precision_t* ts) {
  // matrix inversion constants
  precision_t a(1.0 + _hf * _dt);
  precision_t b(-_hf * _dt);
  precision_t c(-_hs * _dt);
  precision_t d(1.0 + _hs * _dt);
  precision_t denominator(d * a - b * c);

  precision_t tmp_tf(0.0);

// TODO(dave): check aliasing
#ifdef __INTEL_COMPILER
#pragma ivdep
#elif __GNUC__
#pragma GCC ivdep
#endif
  for (int i = 0; i < _numcells; ++i) {
    tmp_tf = (tf[i] * d - ts[i] * b) / denominator;
    ts[i] = (ts[i] * a - tf[i] * c) / denominator;
    tf[i] = tmp_tf;
  }
}

#ifdef TESTING
void Pdesolver::testing() {
  const int n = 1;
  _k = 2 * __SC_PI * n / _simenv->_storage_height;
  _n = n;

  precision_t* errorf = new precision_t;
  precision_t* errors = new precision_t;
  precision_t* iters = new precision_t;
  precision_t* iterf = new precision_t;
  verify(errorf, errors, iterf, iters);

  std::cout << "ERRF: " << *errorf << "\t\tITERF: " << *iterf << std::endl;
  std::cout << "ERRS: " << *errors << "\t\tITERS: " << *iters << std::endl;

  // output
  std::string filename("testing_OVS_r_" + std::to_string(_simenv->_runhash) +
                       "_f.csv");
  std::string fullpath(_simenv->_outfolder + filename);
  std::ofstream fs;
  fs.open(fullpath, std::ofstream::out | std::ofstream::app);
  fs << _simenv->_numcells << ";" << *errorf << "\n";
  fs.close();
  filename = "testing_OVS_r_" + std::to_string(_simenv->_runhash) + "_s.csv";
  fullpath = _simenv->_outfolder + filename;
  fs.open(fullpath, std::ofstream::out | std::ofstream::app);
  fs << _simenv->_numcells << ";" << *errors << "\n";
  fs.close();
}

bool Pdesolver::verify(precision_t* errorf, precision_t* errors,
                       precision_t* iterf, precision_t* iters) {
  const auto sol = [this](precision_t x) { return std::cos(_k * x); };

  precision_t* solution =
      (precision_t*)_mm_malloc(sizeof(precision_t) * _numcells, 32);
  precision_t* fluid_temperature =
      (precision_t*)_mm_malloc(sizeof(precision_t) * _numcells, 32);
  precision_t* solid_temperature =
      (precision_t*)_mm_malloc(sizeof(precision_t) * _numcells, 32);
  precision_t* fluid_temperature_o =
      (precision_t*)_mm_malloc(sizeof(precision_t) * _numcells, 32);
  precision_t* solid_temperature_o =
      (precision_t*)_mm_malloc(sizeof(precision_t) * _numcells, 32);
  _source_fluid = (precision_t*)_mm_malloc(sizeof(precision_t) * _numcells, 32);
  _source_solid = (precision_t*)_mm_malloc(sizeof(precision_t) * _numcells, 32);

  // compute source terms
  // fill solution and initial values to arrays
  precision_t x(_dx * 0.5);
  precision_t leftboundary(sol(x - 0.5 * _dx));
  for (int i = 0; i < _numcells; ++i) {
    solution[i] = sol(x);
    fluid_temperature[i] = sol(x) + 0.1337;
    solid_temperature[i] = sol(x) + 0.1337;
    _source_solid[i] = _dt * (_alphas * _k * std::sin(_k * x));
    _source_fluid[i] = _dt * (-_uf * _k * std::sin(_k * x) +
                              _alphaf * _k * _k * std::cos(_k * x));
    x += _dx;
  }

  // loop while solution not converged
  precision_t diff(1.0);
  int i(0);
  for (; i < _maxiterations && diff > _tol; ++i) {
    solvefluid(&fluid_temperature, &fluid_temperature_o, leftboundary, 0);
    diff = 0.0;
    for (int j = 0; j < _numcells; ++j) {
      diff += std::abs(fluid_temperature[j] - fluid_temperature_o[j]);
    }
    diff = diff / static_cast<precision_t>(_numcells);

    //*DEBUG*/ std::cout << "DIFF: " << diff << "\n";
  }

  *iterf = i;

  // loop while solution not converged
  diff = 1.0;
  i = 0;
  for (; i < _maxiterations && diff > _tol; ++i) {
    solvesolid(&solid_temperature, &solid_temperature_o);
    diff = 0.0;
    for (int j = 0; j < _numcells; ++j) {
      diff += std::abs(solid_temperature[j] - solid_temperature_o[j]);
    }
    diff = diff / static_cast<precision_t>(_numcells);

    //*DEBUG*/ std::cout << "DIFF: " << diff << "\n";
  }

  *iters = i;

  // compute difference to solution
  precision_t difff_solution(0.0);
  precision_t diffs_solution(0.0);
  for (int i = 0; i < _numcells; ++i) {
    difff_solution +=
        std::abs((fluid_temperature[i] - solution[i]) / solution[i]);
    diffs_solution +=
        std::abs((solid_temperature[i] - solution[i]) / solution[i]);
  }

  *errorf = difff_solution / static_cast<precision_t>(_numcells);
  *errors = diffs_solution / static_cast<precision_t>(_numcells);

  _mm_free(fluid_temperature);
  _mm_free(solid_temperature);
  _mm_free(fluid_temperature_o);
  _mm_free(solid_temperature_o);
  _mm_free(solution);
  _mm_free(_source_fluid);
  _mm_free(_source_solid);
};

#endif