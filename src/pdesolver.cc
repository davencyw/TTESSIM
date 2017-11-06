#include "pdesolver/pdesolver.hh"
#include "global.hh"
#include "const.hh"
#include "simenv.hh"

#include <stdio.h>

#include <fstream>
#include <algorithm>
#include <cmath>


void Pdesolver::solvefluid(precision_t* __restrict__ fluid_temperature, precision_t* __restrict__ fluid_temperature_o){
	
	//Loop over inner N-2 cells

	#ifdef __INTEL_COMPILER
	#pragma ivdep
	#elif __GNUC__
	#pragma GCC ivdep
	#endif

	for (int i = 1; i < _simenv->_numcells-1; ++i)
	{		const precision_t tfi = fluid_temperature[i];
			const precision_t tfim1 = fluid_temperature[i-1];
			const precision_t tfip1 = fluid_temperature[i+1];

			//TODO(dave): Optimize, verify
			fluid_temperature_o[i] = tfi - _dt * (_uf*_idx * (tfi - tfim1)) + _alphafidx2dt * (tfim1 - 2 * tfi + tfip1);
			#ifdef TESTING
				fluid_temperature_o[i] +=  _uf * std::sin(_k*i*_dx) - _alphaf * _k * _k * std::cos(_k*i*_dx);
			#endif
	}

	//TODO(dave):
	//Boundary cells

	//swap pointers
	std::swap(fluid_temperature,fluid_temperature_o);
}


void Pdesolver::solvesolid(precision_t* __restrict__ solid_temperature, precision_t* __restrict__ solid_temperature_o){

	//Loop over inner N-2 cells
	
	#ifdef __INTEL_COMPILER
	#pragma ivdep
	#elif __GNUC__
	#pragma GCC ivdep
	#endif

	for (int i = 0; i < _simenv->_numcells-1; ++i)
	{
		const precision_t tsi = solid_temperature[i];
		const precision_t tsim1 = solid_temperature[i-1];
		const precision_t tsip1 = solid_temperature[i+1];
		//TODO(dave): Optimize, verify
		solid_temperature_o[i] = tsi + _alphafidx2dt * (tsip1 - 2 * tsi + tsim1);
	}

	//TODO(dave):
	//Boundary cells

	//swap pointers
	std::swap(solid_temperature, solid_temperature_o);
}

	
#ifdef TESTING

void Pdesolver::testing(){

	//TODO(Dave): Order verification study fluid.
	//TODO(Dave): Order verification study solid.

	//write fluid
	std::string filename("testing_OVS_r_"+ std::to_string(_simenv->_runhash) + "_f.csv");
	std::string fullpath(_simenv->_outfolder + filename);
	std::ofstream fs;
	fs.open(fullpath, std::ofstream::out | std::ofstream::app);


	//write solid
	filename = "testing_OVS_r_" + std::to_string(_simenv->_runhash) + "_s.csv";
	fullpath = _simenv->_outfolder + filename;
	fs.open(fullpath, std::ofstream::out | std::ofstream::app);

}	

bool Pdesolver::verifyfluid(const int n, precision_t* error){

		_k = 2 * __SC_PI * n * _simenv->_storage_height;
		_n = n;
		const auto solution = [this] (precision_t x) {return std::cos(_k * x);};
		const auto slack = [this] (precision_t x, precision_t uf, precision_t alphaf) {return uf * std::sin(_k*x) - alphaf * _k * _k * std::cos(_k*x);};

  		precision_t* fluid_temperature 	= (precision_t *) _mm_malloc(sizeof(precision_t)*_numcells, 32);
  		precision_t* fluid_temperature_o 	= (precision_t *) _mm_malloc(sizeof(precision_t)*_numcells, 32);
  		precision_t* fluid_solution 	= (precision_t *) _mm_malloc(sizeof(precision_t)*_numcells, 32);


  		//Fill solution and initial values to arrays
  		precision_t x(0.0);
  		for (int i = 0; i < _numcells; ++i)
  		{	
  			fluid_temperature[i] = solution(x);
  			fluid_solution[i] = solution(x);
  			x += _dx;
  		}

  		precision_t diff(1.0);

  		//Loop while solution not stable
  		for (int i = 0; i < _maxiterations && diff > _tol; ++i)
  		{
  			solvefluid(fluid_temperature,fluid_temperature_o);
  			diff = 0.0;
  			for (int j = 0; j < _numcells; ++j)
  			{
  				diff += std::abs(fluid_temperature[j] - fluid_temperature_o[j]);
  			}
  			std::swap(fluid_temperature,fluid_temperature_o);
  		}

  		//compute difference to solution
  		precision_t diff_solution(0.0);
  		for (int i = 0; i < _numcells; ++i)
  		{
  			diff_solution += std::pow(fluid_temperature_o[i] - fluid_solution[i],2);
  		}

  		*error = diff_solution;

  		_mm_free(fluid_temperature);
  		_mm_free(fluid_temperature_o);
  		_mm_free(fluid_solution);

	};

	bool Pdesolver::verfiysolid(const int n, precision_t* error){
		//TODO(dave): implement
	};
#endif