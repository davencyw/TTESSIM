#include "pdesolver/pdesolver.hh"
#include "global.hh"
#include "const.hh"
#include "simenv.hh"

#include <stdio.h>

#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>


//TODO(dave): Optimize arithmetic!! (reordering)
//TODO(dave): Solve swap problem in solvefluid/solid!!!

void Pdesolver::solvefluid(precision_t* __restrict__ fluid_temperature, precision_t* __restrict__ fluid_temperature_o, precision_t boundary){
	
	
	//*DEBUG*/std::cout<<"FDDIFF: "<<fluid_temperature[3] - fluid_temperature[4]*2 + fluid_temperature[5] <<"\n";
	//*DEBUG*/std::cout<<"FD03: "<<fluid_temperature[3] <<"\n";
	//*DEBUG*/std::cout<<"FD04: "<<fluid_temperature[4] <<"\n";
	//*DEBUG*/std::cout<<"FD05: "<<fluid_temperature[5] <<"\n\n";

	#ifdef __INTEL_COMPILER
	#pragma ivdep
	#elif __GNUC__
	#pragma GCC ivdep
	#endif
	for (int i = 1; i < _simenv->_numcells-1; ++i)
	{		const precision_t tfi = fluid_temperature[i];
			const precision_t tfim1 = fluid_temperature[i-1];
			const precision_t tfip1 = fluid_temperature[i+1];

			fluid_temperature_o[i] = tfi - _dt * _uf*_idx * (tfi - tfim1) + _alphafidx2dt * (tfim1 - 2 * tfi + tfip1);

			#ifdef TESTING
				fluid_temperature_o[i] -=  _dt* (_uf * std::cos(_k*i*_dx) - _alphaf  * _k * std::sin(_k*i*_dx));
				//fluid_temperature_o[i] -=  _dt* (_uf * _k * std::sin(_k*i*_dx) + _alphaf  * _k * _k * std::cos(_k*i*_dx));
			#endif
	}

	//Boundary cells
	fluid_temperature_o[0] = fluid_temperature[0] - _dt *_uf*_idx*(fluid_temperature[0]-boundary) + _alphafidx2dt*(fluid_temperature[1]-fluid_temperature[0]);
	 
	#ifdef TESTING
		//TODO(dave): simplify sin/cos
		fluid_temperature_o[0]-= _dt * (_uf * std::cos(0) - _alphaf *  _k * std::sin(0));
		//fluid_temperature_o[0]-= _dt * (_k * _uf * std::sin(0) + _alphaf * _k *  _k * std::cos(0));
	#endif

	const int N(_simenv->_numcells-1);
	const int Nm1(_simenv->_numcells-2);
	fluid_temperature_o[N] =  fluid_temperature[N] - _dt*_uf*_idx*(fluid_temperature[N]-fluid_temperature[Nm1]) + _alphafidx2dt*(fluid_temperature[N]-fluid_temperature[N-1]);

	#ifdef TESTING
		//TODO(dave): simplify sin/cos
		fluid_temperature_o[N]-= _dt * (_uf * std::cos(_k*N*_dx) - _alphaf * _k * std::sin(_k*N*_dx));
		//fluid_temperature_o[N]-= _dt * (_uf * _k * std::sin(_k*N*_dx) + _alphaf * _k * _k * std::cos(_k*N*_dx));
	#endif
	
	//TODO(dave): swap pointers
}


void Pdesolver::solvesolid(precision_t* __restrict__ solid_temperature, precision_t* __restrict__ solid_temperature_o, precision_t boundary){

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
		//TODO(dave): alphas
		solid_temperature_o[i] = tsi + _alphafidx2dt * (tsip1 - 2 * tsi + tsim1);
		//TODO(dave): Testing
	}

	//Boundary cells
	solid_temperature_o[0] = solid_temperature[0] + _alphafidx2dt*(solid_temperature[1]-solid_temperature[0]);
	//TODO(dave): Testing
	//TODO(dave): alphas
	const int N(_simenv->_numcells-1);
	const int Nm1(_simenv->_numcells-2);
	solid_temperature_o[N] = solid_temperature[N] + _alphafidx2dt*(solid_temperature[N]-solid_temperature[Nm1]);
	//TODO(dave): Testing
	//TODO(dave): alphas

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
	precision_t* error = new precision_t;

	verifyfluid(1,error);
	std::cout<<"ERR: "<<*error<<std::endl;
	fs.close();

	//write solid
	filename = "testing_OVS_r_" + std::to_string(_simenv->_runhash) + "_s.csv";
	fullpath = _simenv->_outfolder + filename;
	fs.open(fullpath, std::ofstream::out | std::ofstream::app);
	verfiysolid(1,error);
	fs.close();

}	

bool Pdesolver::verifyfluid(const int n, precision_t* error){

		_k = 2 * __SC_PI * n / _simenv->_storage_height;
		_n = n;
		const auto solution = [this] (precision_t x) {return std::cos(_k * x);};

  		precision_t* fluid_temperature 	= (precision_t *) _mm_malloc(sizeof(precision_t)*_numcells, 32);
  		precision_t* fluid_temperature_o 	= (precision_t *) _mm_malloc(sizeof(precision_t)*_numcells, 32);
  		precision_t* fluid_solution 	= (precision_t *) _mm_malloc(sizeof(precision_t)*_numcells, 32);


  		//Fill solution and initial values to arrays
  		precision_t leftboundary(solution(-_dx));
  		precision_t x(0.0);
  		for (int i = 0; i < _numcells; ++i)
  		{	
  			fluid_solution[i] = fluid_temperature[i] = solution(x);
  			x += _dx;
  		}

  		precision_t diff(1.0);

  		//Loop while solution not converged
  		int i(0);
  		for (; i < _maxiterations && diff > _tol; ++i)
  		{
  			solvefluid(fluid_temperature,fluid_temperature_o,leftboundary);
  			diff = 0.0;
  			for (int j = 0; j < _numcells; ++j)
  			{
  				diff += std::abs(fluid_temperature[j] - fluid_temperature_o[j]);
  			}
  			diff = diff/static_cast<precision_t>(_numcells);

  			//TODO(dave): Check pointerswap with solvefluid
  			std::swap(fluid_temperature,fluid_temperature_o);
  			/*DEBUG*/std::cout<<"DIFF: "<<diff<<"\n";
  		}

  		/*DEBUG*/std::cout<<"ITER: "<< i << std::endl;

  		//compute difference to solution
  		precision_t diff_solution(0.0);
  		for (int i = 0; i < _numcells; ++i)
  		{
  			diff_solution += std::abs((fluid_temperature_o[i] - fluid_solution[i])/fluid_solution[i]);
  		}

  		*error = diff_solution/static_cast<precision_t>(_numcells);

  		_mm_free(fluid_temperature);
  		_mm_free(fluid_temperature_o);
  		_mm_free(fluid_solution);

	};

	bool Pdesolver::verfiysolid(const int n, precision_t* error){
		//TODO(dave): implement
	};
#endif