/**
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file EOS_H2_H2O_ISOT.h
 *
 * Created on 2014-05-12 by Yonghui HUANG
 */

#pragma once


#include "math.h"
#include <algorithm>
#include "AbstractEOS_Flash2PForm.h"
#include <cmath>


/**
 *  EOS formulation for calculating each secondary variable
 *  PC PG PL N_L N_G Kr_g Kr_l X_m X_M Sg
 */
class EOS_Flash2PForm : public AbstractEOS_Flash2PForm
{
public:
    
	/**
	  * Constructor
	  */
	EOS_Flash2PForm() : AbstractEOS_Flash2PForm(3)
    {
    };

	/**
	  * DESTRUCTOR
	  */
	virtual ~EOS_Flash2PForm()
    {

    };

	/**
	  * realization of the eval function. 
	  * this function will evaluate the resiudal of the governing equation, 
	  * based on the unknown values given. 
	  */
	virtual void eval(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalVector & res)
	{
		// getting unknowns
		double Sg = vec_unknowns(0);
		double rho_L_h = vec_unknowns(1);
		double rho_G_h = vec_unknowns(2);

		// calculating residual
		res(0) = Calc_Res_Sg(Sg, rho_L_h, rho_G_h);
		res(1) = Calc_Res_rho_L_h(Sg, rho_L_h);
		res(2) = Calc_Res_rho_G_h(Sg, rho_G_h);
	};

	/**
	* realization of the calc_Jacobian function.
	* this function will evaluate the Jacobian matrix,
	* based on the unknown values given.
	*/
	
	virtual void calc_Jacobian(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalMatrix & J)
	{
		assert(J.cols() == this->get_n() && J.rows() == this->get_n());

		// getting unknowns
		double Sg = vec_unknowns(0);
		double rho_L_h = vec_unknowns(1);
		double rho_G_h = vec_unknowns(2);
		double PG(0.0);
		PG = getPG(Sg);
		double PG_h = Function_PG_h(PG);
		double F1 = Sg;
		double G1 = std::min(C_h*PG_h, X_L) - rho_L_h;
		double F2 = 1 - Sg;
		double G2 = rho_G_h - std::max(C_v*PG_h, X_L);
		// evaluate J
		J.setZero();
		J(0, 0) = rho_L_h - rho_G_h;//-----dF(1)/dSg
		J(0, 1) = -(1 - Sg);
		J(0, 2) = -Sg;

		if (F1 <= G1) {
			J(1, 0) = 1.0;
			J(1, 1) = 0.0;
			J(1, 2) = 0.0;
		}
		else{
			if (C_h*PG_h <= X_L)
			{
				J(1, 0) = C_h*Deriv_dPGH_dPG(Sg)*Deriv_dPGdSg(Sg);
				J(1, 1) = -1.0;
				J(1, 2) = 0.0;
			}
			else{
				J(1, 0) = 0.0;
				J(1, 1) = -1.0;
				J(1, 2) = 0.0;
			}

		}
		if (F2 <= G2) {
			J(2, 0) = -1.0;
			J(2, 1) = 0.0;
			J(2, 2) = 0.0;
		}
		else{
			if (C_v*PG_h >= X_L)
			{
				J(2, 0) = -C_v*Deriv_dPGH_dPG(Sg)*Deriv_dPGdSg(Sg);
				J(2, 1) = 0.0;
				J(2, 2) = 1.0;
			}
			else{
				J(2, 0) = 0.0;
				J(2, 1) = 0.0;
				J(2, 2) = 1.0;
			}

		}
	};
	
	
	virtual void calc_Rest_SecVar(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalVector & vec_rest_var)
	{
		
	};
	/**
	  * realization of the calc_Jacobian function.
	  * this function will evaluate the Jacobian matrix,
	  * based on the unknown values given.
	  */


	virtual void set_env_condition(ogsChem::LocalVector & env_conditon)
	{
		// mean pressure
		P_L = env_conditon(0);
		// molar fraction of the light component
		X_L = env_conditon(1);
	}; 
	/**
	* realization of the calc the residual of saturation
	* EOS--Function 1
	* X=Sl*rho_L_h+Sg*rho_G_h
	*/
	virtual double Calc_Res_Sg(double Sg, double rho_L_h, double rho_G_h)
	{
		double Res_Sg(0.0);
		Res_Sg = X_L - ((1 - Sg)*rho_L_h + Sg*rho_G_h);
		return Res_Sg;
	}
	/**
	* realization of the calc the residual of rho_L_h
	* EOS--Function 2
	* rho_L_h=min(C_h*p_g^h,X_L)
	*/
	virtual double Calc_Res_rho_L_h(double Sg, double rho_L_h)
	{
		
	}
	/**
	* realization of the calc the residual of rho_L_h
	* EOS--Function 3
	* rho_G_h=MAX(C_V*p_g^h,X_L)
	*/
	virtual double Calc_Res_rho_G_h(double Sg, double rho_G_h)
	{
		
	}

	

	virtual double computeFugacityCoefficient(double P, double Vm, double Tc,double pc,double w)
	{
		double Z = P*Vm / (R*T);
		double m = 0.37464 + 1.54226*w - 0.26992*pow(w , 2);
		double alfa = pow((1 + m*(1 - pow((T / Tc),0.5))), 2);
		double ai = 0.45724*pow(R,2)*pow(Tc , 2) / pc*alfa;
		double bi = 0.07780*R*Tc / pc;
		double bbar = bi;
		return 0.0;
			//exp((Z - 1)*bbar / b - log((V - b)*Z / V) + (a / (b*R*T)) / (e - s)*log((V + s*b) / (V + e*b))*...
			//(1 + abar / a - bbar / b));
	}
	/**
	* water vapor pressure of pure water
	* Function of Temperature
	* comes from Amaziane's paper
	*/
	virtual double get_Vapor_Pressure(double T)
	{
		// Here unit of T is Celsius;
		double P_vapor(0.0);
		/*double tmp(0.0);
		tmp = 2.786 + 0.031514*T - 1.2373e-4*pow(T, 2) + 4.2267e-7*pow(T, 3) - 8.1308e-10*pow(T, 4);
		P_vapor = pow(10, tmp);*/
		return P_vapor;
	}

	virtual double func_wilson(double P, double T, double pc, double Tc, double w)
	{
		return (pc / P)*exp(5.37*(1 + w)*(1 - (Tc / T)));
	}

	virtual double func_Rachford_Rice(double vf, MathLib::LocalVector & Z, MathLib::LocalVector & K)
	{
		double object(0.0);
		for (int nn = 0; nn < Z.size(); nn++){
			object += Z(nn) / (1 + (K(nn) - 1)*vf);
		}
	}
	virtual double solve_Rachford_Rice(double vf, MathLib::LocalVector & Z, MathLib::LocalVector & K){
		double eps(1e-6);
		double vf_pre, vf_cur;
		vf_pre = vf;
		
		for (int i = 0; i < 10; i++){
			if (i > 10){
				WARN("Solving Rachford_Rice equation does not converge! \n Using old values as seoncdary varibales. \n");
			}
			double f = func_Rachford_Rice(vf_pre, K, Z);
			double df=(func_Rachford_Rice(vf_pre + eps, Z,K) - func_Rachford_Rice(vf_pre - eps,Z,K)) / 2 / eps;
			double delta_vf = f / df;
			vf_cur = vf_pre - delta_vf;
			if (abs(vf_cur - vf_pre) / vf_pre < 1e-4){
				return vf_cur;
				break;
			}
			else
				vf_pre = vf_cur;
		}
	}
	virtual double  computefugacitycoeff(double P, double T,
		MathLib::LocalVector & X, MathLib::LocalVector & K, 
		MathLib::LocalVector & pc, MathLib::LocalVector & Tc, 
		MathLib::LocalVector & Acentric,MathLib::LocalMatrix & Theta)//calculate the cubic equation
	{
		

	}
	virtual void updateaCache()
	{
		for (int compIIdx = 0; compIIdx < numComponents; ++compIIdx) {
	         for (int compJIdx = 0; compJIdx < numComponents; ++compJIdx) {
				
				double Psi = FluidSystem::interactionCoefficient(compIIdx, compJIdx);
				aCache_[compIIdx][compJIdx] =
				std::sqrt(this->pureParams_[compIIdx].a()
					* this->pureParams_[compJIdx].a())
					* (1 - Psi);
			}
		}
	}
	virtual void updatepureParam(double P, double T,
		MathLib::LocalVector & X, MathLib::LocalVector & K,
		MathLib::LocalVector & pc, MathLib::LocalVector & Tc,
		MathLib::LocalVector & Acentric, MathLib::LocalMatrix & Theta, MathLib::LocalMatrix & matpureParam)
	{
		MathLib::LocalVector m, Alfa, ai, bi;
		MathLib::LocalMatrix Q;
		for (int nn = 0; nn < X.size(); nn++){
			m(nn) = 0.37464 + 1.54226*Acentric(nn) - 0.26992*pow(Acentric(nn), 2);
			Alfa(nn) = pow((1 + m(nn)*(1 - pow((T / Tc(nn)), 0.5))), 2);
			ai(nn) = 0.45724*(pow(R, 2))*(pow(Tc(nn), 2)) / pc(nn)*Alfa(nn);
			bi(nn) = 0.07780*R*Tc(nn) / pc(nn);
		}
		matpureParam.col(0) = ai;
		matpureParam.col(1) = bi;
	
	}
	virtual void updatemixParam()
	{

	}
	virtual void NsPol3(double p, double q, double r, MathLib::LocalVector& roots)
	{
		double eps = 7E-15;
		double a, b, h, phi, D, z[3];
		double pi = 3.1415926535897;
		double nz;
		int i;

		b = (p / 3) * (p / 3);
		a = q / 3 - b;
		b = b * p / 3 + 0.5 * (r - p / 3 * q);
		h = sqrt(fabs(a));

		if (b < 0)
			h = -h;

		D = MathLib::fastpow(a, 3) + b * b;

		if (D <= (-eps))
		{
			nz = 3;
			phi = acos(b / MathLib::fastpow(h, 3)) / 3;
			z[0] = 2 * h * cos(pi / 3 - phi) - p / 3;
			z[1] = 2 * h * cos(pi / 3 + phi) - p / 3;
			z[2] = -2 * h * cos(phi) - p / 3;
		}
		else if (D < eps)
		{
			nz = 3;
			z[0] = -2 * h - p / 3;
			z[1] = h - p / 3;
			z[2] = z[1];
		}
		else
		{
			nz = 1;
			if (a >= eps)
			{
				b = b / MathLib::fastpow(h, 3);
				phi = log(b + sqrt(b * b + 1)) / 3;
				z[0] = -2 * h * sinh(phi) - p / 3;
			}
			else if (a > (-eps))
			{
				z[0] = pow((2 * abs(b)), 1. / 3.);
				if (b > 0)
					z[0] = -z[0];
				z[0] = z[0] - p / 3;
			}
			else
			{
				b = b / MathLib::fastpow(h, 3);
				phi = log(b + sqrt(b * b - 1)) / 3;
				z[0] = -2 * h * cosh(phi) - p / 3;
			}
		}

		for (i = 0; i < nz; i++)
			roots->push_back(z[i]);
	}
	template <typename DynamicEigenMatrix>
	void push_back(DynamicEigenMatrix& m, Vector3d&& values, std::size_t row)
	{
		if (row >= m.rows()) {
			m.conservativeResize(row + 1, Eigen::NoChange);
		}
		m.row(row) = values;
	}
private:

	/**
	  * parameter for calculate the molar density of gas phase
      */
	const double R=8.314;
	const double T = 403.15;// [K]
	const double Hen = 7.65e-6; //Henry constant
	//const double P_vapor = 0.0;
	const double C_h = 1.530e-8;// Here I define a constant value C_h which is equal to Hen*M_G(Molar mass of gas) Hen* M_G
	const double C_v = 7.938638137741564e-07;// 7.939211048841233e-07;//M_G/RT
	const double C_w = 3.969319068870782e-06;//M_L/RT
	const double eps = 1e-12;
	const double rho_L_std = 1000;
	const double M_G = 0.002;
	const double M_L = 0.01;
	
	const double numComponents = 8;//8 components
	double P_L;
	double X_L;

	double _alpha;
	double _beta;

	double xi = 1e-5;
};


