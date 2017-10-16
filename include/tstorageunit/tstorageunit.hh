/*__DECLARATION__
 *
 * 
 *      TTESSIM
 *      simulation of operating a thermocline thermal energy storage unit at high temperatures
 *      
 *
 * 
 *      author: david schmidig [     david@davencyw.net   ]
 *      ETH Zurich             [ davschmi@student.ethz.ch ]
 *      DAVENCYW CODE          [        davencyw.net      ]
 */
#ifndef __TSTORAGEUNIT_HH__
#define __TSTORAGEUNIT_HH__

#include "global.hh"
#include "const.hh"
#include "simenv.hh"
#include "pdesolver/pdesolver.hh"

#include <string>
#include <functional>
#include <immintrin.h>

class Tstorageunit
{
public:
	Tstorageunit(const SimEnv& simenv) : _simenv(simenv){

	_state = 0;
  	_fluid_temperature 		= (precision_t *) _mm_malloc(sizeof(precision_t)*_simenv._numcells, 32);
  	_solid_temperature 		= (precision_t *) _mm_malloc(sizeof(precision_t)*_simenv._numcells, 32);
  	_fluid_temperature_o 	= (precision_t *) _mm_malloc(sizeof(precision_t)*_simenv._numcells, 32);
  	_solid_temperature_o 	= (precision_t *) _mm_malloc(sizeof(precision_t)*_simenv._numcells, 32);

  	 _pdesolver = Pdesolver(&simenv);

	};

	~Tstorageunit(){

		_mm_free(_fluid_temperature);
		_mm_free(_fluid_temperature_o);
		_mm_free(_solid_temperature);
		_mm_free(_solid_temperature_o);
	}

	//step-by-step simulation
	bool simstep();
	bool simsteps(const int);
	bool simsteps(const int, const int);

	const int getstate();
	const bool isidle();

	//output functions
	const bool writetocsv();

private:
	//simulation environment
	const SimEnv _simenv;

	//physical properties of system
	precision_t _total_time;
	precision_t* _fluid_temperature;
	precision_t* _solid_temperature;
	precision_t* _fluid_temperature_o;
	precision_t* _solid_temperature_o;
	//state
	//0 = charging (+)
	//1 = discharging (-)
	//2 = idle between charging and discharging (+-)
	//3 = idle between discharging and charging (-+)
	int _state;

	Pdesolver _pdesolver;

};


#endif //__TSTORAGEUNIT_HH__