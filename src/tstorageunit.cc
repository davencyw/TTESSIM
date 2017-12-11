#include "tstorageunit/tstorageunit.hh"
#include "const.hh"
#include "global.hh"
#include "simenv.hh"

#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>

bool Tstorageunit::run(precision_t time) {}

bool Tstorageunit::simstep() {
  // TODO(dave): implement

  _pdesolver.solvesolid(&_solid_temperature_ptr, &_solid_temperature_o_ptr,
                        _solid_diffusionnumber);
  _pdesolver.solvefluid(&_fluid_temperature_ptr, &_fluid_temperature_o_ptr,
                        _cflnumber, _fluid_diffusionnumber,
                        _boundary_temperature);
}

bool Tstorageunit::simsteps(const int steps) { return simsteps(steps, 0); }

// writes output every outputnstep-step. If outputnstep is 0, nothing will be
// written
bool Tstorageunit::simsteps(const int steps, const int outputnstep) {
  int opns(outputnstep);

  if (!outputnstep) opns = steps;

  // avoid if with double loop
  for (int stepi(0); stepi < steps; stepi += opns) {
    for (int stepj(0); stepj < opns; ++stepj) {
      simstep();
    }
    // TODO(dave): Make asynchronous
    // output every outputnstep steps
    writetocsv(&_fluid_temperature, _simenv._numcells, _simenv._fs_fluid);
    writetocsv(&_solid_temperature, _simenv._numcells, _simenv._fs_solid);
  }
}

const int Tstorageunit::getstate() { return _state; }

const bool Tstorageunit::isidle() { return _state > 1; }

const bool Tstorageunit::writetocsv(array_t* data, int size,
                                    std::ofstream* stream) {
  // TODO(dave): Optimize

  for (int i = 0; i < size - 1; ++i) {
    *stream << (*data)(i) << ";";
  }
  *stream << (*data)(size - 1) << "\n";
};