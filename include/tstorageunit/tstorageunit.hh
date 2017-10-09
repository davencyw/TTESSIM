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

#include <string>
#include <functional>

class Tstorageunit
{
public:
	Tstorageunit(SimEnv& simenv) : _simenv(simenv){

		//TODO(dave): Initialization

	_state = 0;


	};

	//step-by-step simulation
	bool simstep();
	bool simsteps(const int);

	const int getstate();

	//output functions
	const bool writetocsv();

private:
	//simulation environment
	SimEnv _simenv;

	//physical properties of system
	double _total_time;
	//TODO(dave): smart pointers!
	double* _fluid_temperature;
	double* _solid_temperature;
	//state
	//0 = charging (+)
	//1 = discharging (-)
	//2 = idle between charging and discharging (+-)
	//3 = idle between discharging and charging (-+)
	int _state;
};


#endif //__TSTORAGEUNIT_HH__