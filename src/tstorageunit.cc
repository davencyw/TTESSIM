#include "tstorageunit/tstorageunit.hh"
#include "const.hh"
#include "global.hh"
#include "simenv.hh"

#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>

bool Tstorageunit::run(unsigned int cycles) {
  for (int cycle_i = 0; cycle_i < cycles; ++cycle_i) {
    const int state(getstate());

    unsigned int steps(0);
    if (state == 0) {
      updatecfl(std::abs(_uf));
      steps =
          static_cast<unsigned int>(_simenv._timedurstate0 / _simenv._deltat);
    } else if (state == 1) {
      updatecfl(0.0);
      steps =
          static_cast<unsigned int>(_simenv._timedurstate1 / _simenv._deltat);
    } else if (state == 2) {
      updatecfl(-std::abs(_uf));
      steps =
          static_cast<unsigned int>(_simenv._timedurstate2 / _simenv._deltat);
    } else if (state == 3) {
      updatecfl(0.0);
      steps =
          static_cast<unsigned int>(_simenv._timedurstate3 / _simenv._deltat);
    }

    // do stepping in state
    simsteps(steps, _simenv._ops);
  }
}

bool Tstorageunit::simstep() {
  // TODO(dave): implement

  _pdesolver.solvesolid(&_solid_temperature_ptr, &_solid_temperature_o_ptr,
                        _solid_diffusionnumber);
  _pdesolver.solvefluid(&_fluid_temperature_ptr, &_fluid_temperature_o_ptr,
                        _cflnumber, _fluid_diffusionnumber,
                        _boundary_temperature);
  _pdesolver.solvecoupling(_solid_temperature_ptr, _fluid_temperature_ptr);

  _total_time += _simenv._deltat;
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

const int Tstorageunit::getstate() {
  unsigned int cycles_passed(std::floor(_total_time / _time_per_cycle));
  precision_t time_in_cycle(_total_time - cycles_passed * _time_per_cycle);

  assert(time_in_cycle <= _time_per_cycle);

  if (time_in_cycle < _simenv._timedurstate0) {
    return 0;
  } else if (time_in_cycle < _simenv._timedurstate1) {
    return 1;
  } else if (time_in_cycle < _simenv._timedurstate2) {
    return 2;
  } else {
    return 3;
  }
}

void Tstorageunit::updatecfl(precision_t uf) {
  _cflnumber = _uf * _simenv._deltat / _dx;
  _pdesolver.updateuf(_uf);
}

const bool Tstorageunit::writetocsv(array_t* data, int size,
                                    std::ofstream* stream) {
  // TODO(dave): Optimize

  for (int i = 0; i < size - 1; ++i) {
    *stream << (*data)(i) << ";";
  }
  *stream << (*data)(size - 1) << "\n";
};