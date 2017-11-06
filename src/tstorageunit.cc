#include "tstorageunit/tstorageunit.hh"
#include "global.hh"
#include "const.hh"
#include "simenv.hh"

#include <fstream>
#include <iostream>
#include <functional>
#include <ctime>

bool Tstorageunit::simstep(){
	//TODO(dave): implement
}

bool Tstorageunit::simsteps(const int steps){
	return simsteps(steps,0);
}

//writes output every outputnstep-step. If outputnstep is 0, nothing will be written
bool Tstorageunit::simsteps(const int steps, const int outputnstep){
	//TODO(dave): output(writetocsv)
	int opns(outputnstep);

	if(!outputnstep)
		opns = steps;

	//avoid if with double loop
	for(int stepi(0); stepi < steps; stepi += opns){
		for(int stepj(0); stepj < opns; ++stepj){
		simstep();
		}
		//output every outputnstep steps
		writetocsv(_fluid_temperature, _simenv._numcells, true);
		writetocsv(_solid_temperature, _simenv._numcells, false);
	}
}

const int Tstorageunit::getstate(){
	return _state;
}

const bool Tstorageunit::isidle(){
	return _state > 1;
}

const bool Tstorageunit::writetocsv(precision_t* data, int size, bool fluid){
	
	//TODO(dave): Optimize

	std::ofstream* fsp = _simenv._fs_solid;
	if(fluid){
		fsp = _simenv._fs_fluid;
	}

	for (int i = 0; i < size-1; ++i)
	{
		*fsp << data[i] << ";";
	}
	*fsp << data[size-1] << "\n";

};