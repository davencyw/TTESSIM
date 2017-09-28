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
#include <ctime>
#include <functional>

class Tstorageunit
{
public:
	Tstorageunit(SimEnv& simenv) : _simenv(simenv){

		//TODO(dave): Initialization

	_simenv._runhash = std::time(0);
	};

	//step-by-step simulation
	bool simstep();
	bool simsteps(int);

	//output functions
	bool writetocsv();

private:
	//simulation environment
	SimEnv _simenv;

	//physical properties of system
	double _total_time;
	//TODO(dave): smart pointers!
	double* _fluid_temperature;
	double* _solid_temperature;
};


#endif //__TSTORAGEUNIT_HH__