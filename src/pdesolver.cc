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
      diffusionnumber * (top - 2.0 * mid + bot);

  (*temperature_old)(_numcells - 1) =
      diffusionnumber *
      ((*temperature)(_numcells - 2) - (*temperature)(_numcells - 1));

  (*temperature_old)(0) =
      diffusionnumber * ((*temperature)(1) - (*temperature)(0));
}

void Pdesolver::solveadvection(array_t* temperature, array_t* temperature_old,
                               precision_t cflnumber,
                               precision_t boundary_temperature) {
  array_t tail = temperature->tail(_numcells - 1);
  array_t head = temperature->head(_numcells - 1);
  array_t diff = head - tail;

  if (cflnumber > 0) {
    (*temperature_old)(0) = boundary_temperature;
    temperature_old->tail(_numcells - 1) = cflnumber * diff;

  } else if (cflnumber < 0) {
    (*temperature_old)(_numcells - 1) = boundary_temperature;
    temperature_old->head(_numcells - 1) = cflnumber * diff;
  } else {
    return;
  }
}

void Pdesolver::solvecoupling(array_t* temperature_solid,
                              array_t* temperature_fluid) {
  precision_t coeff_a(1.0 + _hf * _dt);
  precision_t coeff_b(-(-_hf * _dt));
  precision_t coeff_c(-(-_hs * _dt));
  precision_t coeff_d(1 + _hs * _dt);
  precision_t coeff_determinant(1.0 / (coeff_a * coeff_d - coeff_b * coeff_c));

  *temperature_fluid =
      (*temperature_fluid * coeff_d + *temperature_solid * coeff_b) *
      coeff_determinant;
  *temperature_solid =
      (*temperature_fluid * coeff_c + *temperature_solid * coeff_a) *
      coeff_determinant;
}

void Pdesolver::solvefluid(array_t** temperature, array_t** temperature_old,
                           precision_t cflnumber, precision_t diffusionnumber,
                           precision_t boundary_temperature) {
  assert((*temperature_old)->size() == (*temperature)->size());
  assert((*temperature_old)->size() == _numcells);

  // TODO(dave): add another permanent memory allocator for diffusion_tmp
  // diffusion part
  solvediffusion(*temperature, *temperature_old, diffusionnumber);
  array_t diffusion_tmp = **temperature_old;
  // advection part
  solveadvection(*temperature, *temperature_old, cflnumber,
                 boundary_temperature);

  if (cflnumber > 0) {
    (*temperature_old)->tail(_numcells - 1) +=
        (*temperature)->tail(_numcells - 1);

  } else if (cflnumber < 0) {
    (*temperature_old)->head(_numcells - 1) +=
        (*temperature)->head(_numcells - 1);
  }

  **temperature_old += diffusion_tmp;

  std::swap(*temperature, *temperature_old);
#ifdef TESTING
  **temperature += _source_fluid;
#endif
}

void Pdesolver::solvesolid(array_t** temperature, array_t** temperature_old,
                           precision_t diffusionnumber) {
  solvediffusion(*temperature, *temperature_old, diffusionnumber);
  **temperature_old += **temperature;

  std::swap(*temperature, *temperature_old);
#ifdef TESTING
  **temperature += _source_solid;
#endif
}

void Pdesolver::updateuf(precision_t uf) { _uf = uf; }

#ifdef TESTING

void Pdesolver::verify(precision_t* errorf, precision_t* errors,
                       precision_t* iterf, precision_t* iters) {
  // diffusionnumber < 0.5!
  const precision_t diffusionnumber(_dt / (_dx * _dx) * _alphas);
  assert(diffusionnumber <= 0.5);

  const precision_t cflnumber(_uf * _dt / _dx);

  // create verification environment (grid, solution and source terms)
  array_t grid = array_t::LinSpaced(_numcells, _dx / 2.0,
                                    _simenv->_storage_height - _dx / 2.0);
  array_t gridfaces =
      array_t::LinSpaced(_numcells + 1, 0, _simenv->_storage_height);
  array_t lowergridfaces = gridfaces.head(_numcells);
  array_t uppergridfaces = gridfaces.tail(_numcells);

  _source_solid = (_k * _dt * _alphas / _dx) *
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

    diff = (temperature - temperature_o).abs().sum() /
           static_cast<precision_t>(_numcells);

    ++iter;
  } while (diff > k_tol && iter < k_maxiterations);

  // computer error:
  *errors = (temperature - solution).abs().sum() /
            static_cast<precision_t>(_numcells);
  *iters = iter;

  // verify fluid equations
  _source_fluid =
      _source_solid +
      cflnumber * ((_k * uppergridfaces).cos() - (_k * lowergridfaces).cos());
  temperature = solution;
  temperature_o = array_t::Zero(_numcells);
  temperature_ptr = &temperature;
  temperature_o_ptr = &temperature_o;
  assert(temperature_ptr == &temperature);
  assert(temperature_o_ptr == &temperature_o);

  diff = 0.0;
  iter = 0;
  do {
    solvefluid(&temperature_ptr, &temperature_o_ptr, cflnumber, diffusionnumber,
               1);
    diff = (temperature - temperature_o).abs().sum() /
           static_cast<precision_t>(_numcells);

    ++iter;
  } while (diff > k_tol && iter < k_maxiterations);

  *errorf = (temperature - solution).abs().sum() /
            static_cast<precision_t>(_numcells);
  *iterf = iter;
}

void Pdesolver::testing() {
  const int n = 1;
  _k = 2.0 * __SC_PI * static_cast<precision_t>(n) / _simenv->_storage_height;

  precision_t* errorf = new precision_t;
  precision_t* errors = new precision_t;
  precision_t* iters = new precision_t;
  precision_t* iterf = new precision_t;
  verify(errorf, errors, iterf, iters);

  std::cout.precision(10);
  std::cout << "   N: " << _simenv->_numcells << "\n";
  std::cout << "ERRF: " << *errorf << "\t\tITERF: " << *iterf << "\n";
  std::cout << "ERRS: " << *errors << "\t\tITERS: " << *iters << "\n\n";

  precision_t pecletnumber(_uf * _simenv->_storage_height / _alphas);

  // output
  std::string filename("testing_OVS_r_" + std::to_string(_simenv->_runhash) +
                       "_f.csv");
  std::string fullpath(_simenv->_outfolder + filename);
  std::ofstream fs;
  fs.open(fullpath, std::ofstream::out | std::ofstream::app);
  fs << _simenv->_numcells << ";" << pecletnumber << ";" << *errorf << "\n";
  fs.close();
  filename = "testing_OVS_r_" + std::to_string(_simenv->_runhash) + "_s.csv";
  fullpath = _simenv->_outfolder + filename;
  fs.open(fullpath, std::ofstream::out | std::ofstream::app);
  fs << _simenv->_numcells << ";" << pecletnumber << ";" << *errors << "\n";
  fs.close();
}

#endif