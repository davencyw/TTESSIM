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
#ifndef __PDESOLVER_HH__
#define __PDESOLVER_HH__

#include "global.hh"
#include <functional>
#include "simenv.hh"

class Pdesolver
{
public:

	Pdesolver() = default;
	Pdesolver(const SimEnv* simenv): _simenv(simenv){
		_dt = _simenv->_deltat;
		_dx = _simenv->_storage_height / _simenv->_numcells;
		_dx2 = _dx*_dx;
		_idx = 1.0/_dx;
		_idx2 = 1.0/_dx2;
		_kf = _simenv->_kf;
		_epsilon = _simenv->_epsilon;
		_rhof = _simenv->_rhof;
		_cf = _simenv->_cf;
		_alphaf = _kf / _epsilon * _rhof * _cf;
		_alphafidx2 = _alphaf * _idx2;
		_alphafidx2dt = _alphafidx2 * _dt;
		_uf = _simenv->_uf;
		_numcells = _simenv->_numcells;

	};

	void solvefluid(precision_t* fluid_temperature, precision_t* fluid_temperature_o);
	void solvesolid(precision_t* solid_temperature, precision_t* solid_temperature_o);


	#ifdef DEBUG
	//These functions solve the governing fluid and
	//solid equations for a MMS given the slack term as
	//default or passed lambda expression.
	bool verifyfluid(const int n);
	bool verfiysolid();
	#endif


private:
	const SimEnv* _simenv;

	precision_t _dt;
	precision_t _dx;
	precision_t _dx2;
	precision_t _idx;
	precision_t _idx2;
	precision_t _kf;
	precision_t _epsilon;
	precision_t _rhof;
	precision_t _cf;
	precision_t _alphaf;
	precision_t _alphafidx2;
	precision_t _alphafidx2dt;
	precision_t _uf;
	int 		_numcells;

	//DEBUG
	#ifdef DEBUG
	int _n;
	int _k;
	static constexpr precision_t _tol = 10e-8;
	#endif
};


#endif //__PDESOLVER_HH__