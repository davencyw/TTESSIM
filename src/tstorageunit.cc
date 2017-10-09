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
	//TODO(dave): implement
}

const int Tstorageunit::getstate(){
	return _state;
}

const bool Tstorageunit::writetocsv(){
	std::string filename("r_"+ std::to_string(_simenv._runhash) + "_" + std::to_string(_total_time) + ".csv");
	std::string fullpath(_simenv._outfolder + filename);

	std::ofstream filestream(fullpath, std::ofstream::out);

	//TODO(dave): write data (!!! don't use std::endl, use "\n")
	filestream << "test00,test01;test10,test11;";

	filestream.close();
};