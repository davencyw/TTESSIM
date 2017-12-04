#include "pdesolver/pdesolver.hh"
#include "const.hh"
#include "global.hh"
#include "simenv.hh"

#include <stdio.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

// TODO(dave): Optimize arithmetic!! (reordering)
void Pdesolver::solvefluid(precision_t** ft,
                           precision_t** fto,
                           precision_t boundary, int state) {
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
                             _alphafidx2dt * (tfim1 - 2 * tfi + tfip1);
  }

  const int N(_simenv->_numcells - 1);
  const int Nm1(_simenv->_numcells - 2);
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

// mms
#ifdef TESTING
#ifdef __INTEL_COMPILER
#pragma ivdep
#elif __GNUC__
#pragma GCC ivdep
#endif
  for (int i = 0; i < _simenv->_numcells; ++i) {
    fluid_temperature_o[i] += _source_fluid[i];
  }
#endif

  // TODO(dave): swap pointers
}

void Pdesolver::solvesolid(precision_t* __restrict__ solid_temperature,
                           precision_t* __restrict__ solid_temperature_o,
                           precision_t boundary) {
// Loop over inner N-2 cells

#ifdef __INTEL_COMPILER
#pragma ivdep
#elif __GNUC__
#pragma GCC ivdep
#endif
  for (int i = 0; i < _simenv->_numcells - 1; ++i) {
    const precision_t tsi = solid_temperature[i];
    const precision_t tsim1 = solid_temperature[i - 1];
    const precision_t tsip1 = solid_temperature[i + 1];
    // TODO(dave): Optimize, verify
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
}

#ifdef TESTING

void Pdesolver::testing() {
  const int n = 1;
  _k = 2 * __SC_PI * n / _simenv->_storage_height;
  _n = n;

  precision_t* errorf = new precision_t;
  precision_t* errors = new precision_t;
  verify(errorf, errors);
  std::cout << "ERRF: " << *errorf << std::endl;
  std::cout << "ERRS: " << *errors << std::endl;

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

bool Pdesolver::verify(precision_t* errorf, precision_t* errors) {
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
  precision_t leftboundary(sol(0));
  precision_t x(_dx * 0.5);
  for (int i = 0; i < _numcells; ++i) {
    solution[i] = sol(x);
    fluid_temperature[i] = sol(x);
    solid_temperature[i] = sol(x);
    _source_solid[i] = -_dt * (2.0 * _alphas * _k * _k * std::cos(2 * _k * x));
    _source_fluid[i] = _dt * (-_uf * _k * std::sin(_k * x) +
                              _alphaf * _k * _k * std::cos(_k * x));
    x += _dx;
  }

  // loop while solution not converged
  precision_t diff(1.0);
  int i(0);
  for (; i < _maxiterations && diff > _tol; ++i) {
    solvefluid(fluid_temperature, fluid_temperature_o, leftboundary, 0);
    diff = 0.0;
    for (int j = 0; j < _numcells; ++j) {
      diff += std::abs(fluid_temperature[j] - fluid_temperature_o[j]);
    }
    diff = diff / static_cast<precision_t>(_numcells);

    // TODO(dave): Check pointerswap with solvefluid
    std::swap(fluid_temperature, fluid_temperature_o);
    //*DEBUG*/ std::cout << "DIFF: " << diff << "\n";
  }

  /*DEBUG*/ std::cout << "ITERF: " << i << std::endl;

  // loop while solution not converged
  diff = 1.0;
  i = 0;
  for (; i < _maxiterations && diff > _tol; ++i) {
    solvesolid(solid_temperature, solid_temperature_o, leftboundary);
    diff = 0.0;
    for (int j = 0; j < _numcells; ++j) {
      diff += std::abs(solid_temperature[j] - solid_temperature_o[j]);
    }
    diff = diff / static_cast<precision_t>(_numcells);

    // TODO(dave): Check pointerswap with solvefluid
    std::swap(solid_temperature, solid_temperature_o);
    //*DEBUG*/ std::cout << "DIFF: " << diff << "\n";
  }

  /*DEBUG*/ std::cout << "ITERS: " << i << std::endl;

  // compute difference to solution
  precision_t difff_solution(0.0);
  precision_t diffs_solution(0.0);
  for (int i = 0; i < _numcells; ++i) {
    difff_solution +=
        std::abs((fluid_temperature_o[i] - solution[i]) / solution[i]);
    diffs_solution +=
        std::abs((solid_temperature_o[i] - solution[i]) / solution[i]);
  }

  *errorf = difff_solution / static_cast<precision_t>(_numcells);
  *errors = diffs_solution / static_cast<precision_t>(_numcells);

  _mm_free(fluid_temperature);
  _mm_free(solid_temperature);
  _mm_free(fluid_temperature_o);
  _mm_free(solid_temperature_o);
  _mm_free(solution);
};

#endif