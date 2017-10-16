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
#ifndef __SIMENV_HH__
#define __SIMENV_HH__

#include <string>

#include "global.hh"

struct SimEnv
{
		
	//physical environment
	//____________________

	//number of cells
	int _numcells;
	// height of storage unit (cylindrical)
	double _storage_height;
	//diameter of storage unit (cylindrical)
	double _storage_diameter;
	//initial temperature of fluid
	double _fluid_initemp;

	//times for states
	double _timedurstate0;
	double _timedurstate1;
	double _timedurstate2;
	double _timedurstate3;

	//fluid properties
	double _kf;
	double _ks;
	double _rhof;
	double _rhos;
	double _cf;
	double _cs;
	double _epsilon;


	//computing environment
	//____________________
	
	//number of threads
	int _nThreads;
	//omp scheduling
	int _scheduling;
	//cuda enabled
	bool _cuda;
	//number of cycles
	int _numcycles;
	//tsteps per cycle;
	int _tsteppercycle;

	//data environment
	//____________________

	std::string _outfolder;
	int _runhash;

};

#endif //__SIMENV_HH__