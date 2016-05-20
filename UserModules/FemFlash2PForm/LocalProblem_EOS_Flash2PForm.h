/**
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LocalProblem_EOS.h
 *
 * Created on 2014-05-14 by Yonghui HUANG & Haibing Shao
 */

#ifndef LOCAL_PROBLEM_EOS_FLASH2PFORM_H
#define LOCAL_PROBLEM_EOS_FLASH2PFORM_H

#include "ChemLib/chemconst.h"
#include "AbstractEOS_Flash2PForm.h"

/**
  * TODO: describe this class
  */
class LocalProblem_EOS_Flash2PForm
{
public:
	/**
      * constructor of the class
      */
	LocalProblem_EOS_Flash2PForm();
	
	/**
      * destructor of the class
      */
	~LocalProblem_EOS_Flash2PForm(void);

	/**
	  * solve the EOS
	  * In our case, we have 2 inputs, P and X
	  * and 8 outputs, which are the secondary variables. 
	  * the initial guess will be embedded into the Output vector.
	  */
	void solve(ogsChem::LocalVector & Input, ogsChem::LocalVector & Output );
	//void deriv_solve(ogsChem::LocalVector & Input, ogsChem::LocalVector & Output)
	/**
	  * return the number of secondary variables
	  */
	std::size_t get_n() { return N; };

	/**
	  * TODO: describe the parameter
	  */
	const double eps = 1.0e-7;

	/**
	  * TODO: describe the parameter
	  */
	const double Iter_tol=200;

	/**
	  * TODO: describe the parameter
	  */
	const double theta = 0.7;

	/**
	  * tolorence of newton iteration
	  */
	const double tot = 1.0e-12; //a small tolorence is necessary for the iteration
	
	/**
	  * parameter to calc the N_G
	  */
	const double s= 1.0e-4;

	const double T = 403.15;
	const double R = 8.314;
	MathLib::LocalVector P_crit;
	P_crit << 7.375e+6, 3.39e+6, 4.599e+6, 4.654e+6, 3.609e+6, 2.504e+6, 1.502e+6, 0.76e+6;
	MathLib::LocalVector T_crit = MathLib::LocalVector::Zero(8);
	T_crit << 304.14, 126.21, 190.56, 327.81, 435.62, 574.41, 708.95, 891.47;
	MathLib::LocalVector acentric = MathLib::LocalVector::Zero(8);
	acentric << 0.239, 0.039, 0.011, 0.11783, 0.21032, 0.41752, 0.66317, 1.7276;
private:
	/**
      * parameter to define the number of secondary variables
	  */
	std::size_t N;

	/**
	  * mean pressure ---primary variable
	  */
	double P_L;

	/**
	  * Total hydrogen mass density ---primary variable
	  */
	double X_L;

	/**
	  *  initial guess vector of the secondary variables 
	  */
	ogsChem::LocalVector U_ini;

	/**
	  * internal output vector
	  */
	ogsChem::LocalVector U_cur;

	/**
	  * Internal Residual vector
	  */
	ogsChem::LocalVector Res;

	/**
	  * Internal matrix of Jacobian 
	  */
	ogsChem::LocalMatrix Matrix_Jacobian;

	/**
	  * pointer the abstract EOS class
	  */
	AbstractEOS_Flash2PForm * _EOS;

	/**
	  * minimization solve function. 
	  * it is capable of handling rank(J) < n. 
	  */
	void solve_minimization(ogsChem::LocalMatrix & J,
		ogsChem::LocalVector & b,
		ogsChem::LocalVector & dx);

	/**
	  * Newton iteration with line search
	  */
	void solve_LocalProblem_Newton_LineSearch(std::size_t & flag);
};

#endif
