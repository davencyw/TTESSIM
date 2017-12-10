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

// TODO(dave): Implement alpha
// TODO(dave): write documentation
// TODO(dave): write gridview/gridclass!
// URGENT_TODO(dave): write memory view for segments!! -^^^^

void Pdesolver::solvediffusion(array_t* temperature, array_t* temperature_old,
                               precision_t diffusionnumber) {
  assert(temperature_old->size() == temperature->size());
  assert(temperature_old->size() == _numcells);

  array_t bot = temperature->head(_numcells - 2);
  array_t mid = temperature->segment(1, _numcells - 2);
  array_t top = temperature->tail(_numcells - 2);

  temperature_old->segment(1, _numcells - 2) =
      mid + diffusionnumber * (top - 2.0 * mid + bot);

  (*temperature_old)(_numcells - 1) =
      (*temperature)(_numcells - 1) +
      diffusionnumber *
          ((*temperature)(_numcells - 2) - (*temperature)(_numcells - 1));
  (*temperature_old)(0) =
      (*temperature)(0) +
      diffusionnumber * ((*temperature)(1) - (*temperature)(0));
}

void Pdesolver::solveadvection(array_t* temperature, array_t* temperature_old,
                               precision_t cflnumber,
                               precision_t boundary_temperature) {
  array_t bot = temperature->head(_numcells - 2);
  array_t mid = temperature->segment(1, _numcells - 2);

  if (cflnumber > 0) {
    (*temperature_old)(0) = boundary_temperature;
    temperature_old->segment(1, _numcells - 2) = mid - cflnumber * (mid - bot);

  } else if (cflnumber < 0) {
    (*temperature_old)(_numcells - 1) = boundary_temperature;
    temperature_old->head(_numcells - 1) = bot - cflnumber * (mid - bot);
  } else {
    return;
  }
}

void Pdesolver::solvefluid(array_t** temperature, array_t** temperature_old,
                           precision_t cflnumber, precision_t diffusionnumber,
                           precision_t boundary_temperature) {
  assert((*temperature_old)->size() == (*temperature)->size());
  assert((*temperature_old)->size() == _numcells);

  // diffusion part
  solvediffusion(*temperature, *temperature_old, diffusionnumber);
  // advection part
  solveadvection(*temperature, *temperature_old, cflnumber,
                 boundary_temperature);

  if (cflnumber > 0) {
  } else if (cflnumber < 0) {
  }

  std::swap(temperature, temperature_old);
}

void Pdesolver::solvesolid(array_t** temperature, array_t** temperature_old,
                           precision_t diffusionnumber) {
  solvediffusion(*temperature, *temperature_old, diffusionnumber);
  //(**temperature_old) = (**temperature_old) + (**temperature);

  // (*temperature)->head(_numcells) =
  //     (*temperature_old)->head(_numcells) + (*temperature)->head(_numcells);

  std::swap(*temperature, *temperature_old);
#ifdef TESTING
  **temperature += _source_solid;
#endif
}

#ifdef TESTING

void Pdesolver::verify(precision_t* errorf, precision_t* errors,
                       precision_t* iterf, precision_t* iters) {
  const precision_t alpha(0.1);

  // diffusionnumber < 0.5!
  const precision_t diffusionnumber(_dt / (_dx * _dx) * alpha);
  assert(diffusionnumber <= 0.5);

  // create verification environment (grid, solution and source terms)
  array_t grid = array_t::LinSpaced(_numcells, _dx / 2.0,
                                    _simenv->_storage_height - _dx / 2.0);
  array_t gridfaces =
      array_t::LinSpaced(_numcells + 1, 0, _simenv->_storage_height);
  array_t lowergridfaces = gridfaces.head(_numcells);
  array_t uppergridfaces = gridfaces.tail(_numcells);

  _source_solid = (_k * _dt * alpha / _dx) *
                  ((_k * uppergridfaces).sin() - (_k * lowergridfaces).sin());
  array_t solution = (_k * grid).cos();
  array_t temperature = solution;
  array_t temperature_o = array_t::Zero(_numcells);

  array_t* temperature_ptr = &temperature;
  array_t* temperature_o_ptr = &temperature_o;

  // verify solid equations
  precision_t diff(0.0);
  unsigned int iter(0);
  do {
    solvesolid(&temperature_ptr, &temperature_o_ptr, diffusionnumber);

    diff = (temperature - temperature_o).abs().sum();
    diff /= static_cast<precision_t>(_numcells);

    ++iter;
  } while (diff > k_tol && iter < k_maxiterations);

  // computer error:
  *errors = (temperature - solution).abs().sum() /
            static_cast<precision_t>(_numcells);
  *iters = iter;
}

void Pdesolver::testing() {
  const int n = 1;
  _k = 2.0 * __SC_PI * static_cast<precision_t>(n) / _simenv->_storage_height;
  // _n = n;

  precision_t* errorf = new precision_t;
  precision_t* errors = new precision_t;
  precision_t* iters = new precision_t;
  precision_t* iterf = new precision_t;
  verify(errorf, errors, iterf, iters);

  std::cout << "   N: " << _simenv->_numcells << "\n";
  // std::cout << "ERRF: " << *errorf << "\t\tITERF: " << *iterf << "\n";
  std::cout << "ERRS: " << *errors << "\t\tITERS: " << *iters << "\n\n";

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

#endif