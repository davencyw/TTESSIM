#include "pdesolver/pdesolver.hh"
#include "global.hh"
#include "const.hh"
#include "simenv.hh"

#include <stdio.h>

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

			fluid_temperature_o[i] = tfi - _dt * (_uf*_idx * (tfi - tfim1)) + _alphafidx2dt * (tfim1 - 2 * tfi + tfip1);
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

		solid_temperature_o[i] = tsi + _alphafidx2dt * (tsip1 - 2 * tsi + tsim1);
	}

	//TODO(dave):
	//Boundary cells

	//swap pointers
	std::swap(solid_temperature, solid_temperature_o);
}


	bool Pdesolver::verifyfluid(const int n){

		const int k(2 * __SC_PI * n * _simenv->_storage_height);		
		const auto solution = [&k] (precision_t x) {return std::cos(k * x);};
		const auto slack = [&k] (precision_t x, precision_t uf, precision_t alphaf) {return uf * std::sin(k*x) - alphaf * k * k * std::cos(k*x);};

		//create data and solution
		precision_t* fluid_temperature 		= (precision_t *) _mm_malloc(sizeof(precision_t)*_simenv->_numcells, 32);
  		precision_t* fluid_temperature_o 	= (precision_t *) _mm_malloc(sizeof(precision_t)*_simenv->_numcells, 32);
  		precision_t* solution_temperature	= (precision_t *) _mm_malloc(sizeof(precision_t)*_simenv->_numcells, 32);

  		for (int i = 0; i < _simenv->_numcells; ++i)
  		{
  			solution_temperature[i] = solution(i * _dx);
  			//initial condition
  			fluid_temperature[i] = solution(i * _dx);
  		}

	};
	bool Pdesolver::verfiysolid(){
		//TODO(dave): implement
	};
