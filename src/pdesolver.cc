#include "pdesolver/pdesolver.hh"
#include "global.hh"
#include "const.hh"
#include "simenv.hh"

#include <stdio.h>

#include <algorithm>

bool Pdesolver::solvefluid(precision_t* __restrict__ fluid_temperature, precision_t* __restrict__ fluid_temperature_o){
	
	//Loop over inner N-2 cells
	#pragma GCC ivdep
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


bool Pdesolver::solvesolid(precision_t* __restrict__ solid_temperature, precision_t* __restrict__ solid_temperature_o){

	//Loop over inner N-2 cells
	#pragma GCC ivdep
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