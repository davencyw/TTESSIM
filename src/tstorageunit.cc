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
    for (int phase_i = 0; phase_i < 4; ++phase_i) {
      // TODO(dave): change deltat in cmdparser!!!
      const int state(getstate());
      unsigned int steps(0);
      if (state == 0) {
        // compute before charging
        computecapacityfactor();
        updatecfl(std::abs(_uf));
        _boundary_temperature = _simenv._fluid_temp_charge;
        steps =
            static_cast<unsigned int>(_simenv._timedurstate0 / _simenv._deltat);
      } else if (state == 1) {
        updatecfl(0.0);
        steps =
            static_cast<unsigned int>(_simenv._timedurstate1 / _simenv._deltat);
      } else if (state == 2) {
        // compute before discharging
        computecapacityfactor();
        updatecfl(-std::abs(_uf));
        _boundary_temperature = _simenv._fluid_temp_discharge;
        steps =
            static_cast<unsigned int>(_simenv._timedurstate2 / _simenv._deltat);
      } else if (state == 3) {
        updatecfl(0.0);
        steps =
            static_cast<unsigned int>(_simenv._timedurstate3 / _simenv._deltat);
      }

      // do stepping in statephase
      simsteps(steps, _simenv._ops);
    }
    std::cout << "completed cycle " << cycle_i + 1 << "\n";
  }
}

bool Tstorageunit::simstep() {
  _pdesolver.solvesolid(&_solid_temperature_ptr, &_solid_temperature_o_ptr,
                        _solid_diffusionnumber);
  _pdesolver.solvefluid(&_fluid_temperature_ptr, &_fluid_temperature_o_ptr,
                        _cflnumber, _fluid_diffusionnumber,
                        _boundary_temperature);
  _pdesolver.solvecoupling(_solid_temperature_ptr, _fluid_temperature_ptr);

  _total_time += _simenv._deltat;

  computeefficiency();
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

void Tstorageunit::computeefficiency() {
  int state(getstate());

  if (state == 1 || state == 3) {
    return;
  }

  unsigned int inindex(0);
  unsigned int outindex(_simenv._numcells - 1);

  // discharge hhase has index 1 in _exergy_flux_array and in/out locations are
  // swapped
  if (state == 2) {
    state = 1;
    std::swap(inindex, outindex);
  };

  // in
  _exergy_flux_array(state) +=
      _simenv._deltat *
      (10 * _simenv._cf *
       (_fluid_temperature[inindex] - _simenv._fluid_initemp -
        _simenv._fluid_initemp *
            std::log(_fluid_temperature[inindex] / _simenv._fluid_initemp)));
  // out
  _exergy_flux_array(state + 1) +=
      _simenv._deltat *
      (10 * _simenv._cf *
       (_fluid_temperature[outindex] - _simenv._fluid_initemp -
        _simenv._fluid_initemp *
            std::log(_fluid_temperature[outindex] / _simenv._fluid_initemp)));

  // update exergy efficiency
  _exergy_flux = (_exergy_flux_array(1) - _exergy_flux_array(0)) /
                 (_exergy_flux_array(2) - _exergy_flux_array(3));
}

void Tstorageunit::computecapacityfactor() {
  const int state(getstate());
  assert(state == 0 || state == 2);

  const precision_t prefactor(__SC_PI / 4.0 * _simenv._storage_diameter *
                              _simenv._storage_diameter);
  const precision_t fluid_difference =
      ((_fluid_temperature - _simenv._fluid_temp_charge) * _dx).sum();
  const precision_t solid_difference =
      ((_solid_temperature - _simenv._fluid_temp_discharge) * _dx).sum();
  const precision_t thermal_energy =
      prefactor +
      (_simenv._epsilon * _simenv._rhof * _simenv._cf * fluid_difference +
       (1.0 - _simenv._epsilon) * _simenv._rhos * _simenv._cs *
           solid_difference);

  if (state == 0) {
    // before charging
    _thermal_energy_array(0) = thermal_energy;

  } else {
    // before discharging
    _thermal_energy_array(1) = thermal_energy;
  }

  _capacity_factor = (_thermal_energy_array(1) - _thermal_energy_array(0)) /
                     _max_thermal_energy;
}

const bool Tstorageunit::writetocsv(array_t* data, int size,
                                    std::ofstream* stream) {
  // TODO(dave): Optimize
  for (int i = 0; i < size - 1; ++i) {
    *stream << (*data)(i) << ";";
  }
  *stream << (*data)(size - 1) << "\n";
};