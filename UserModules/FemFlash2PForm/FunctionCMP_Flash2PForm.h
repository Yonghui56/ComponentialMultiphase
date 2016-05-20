/**
* Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file FunctionCMP_2P2C.h
*
* Created on    2014-05-12 by Yonghui HUANG
*/

#ifndef FUNCTIONCMP_FLASH2PFORM_H
#define FUNCTIONCMP_FLASH2PFORM_H

#include "BaseLib/CodingTools.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "SolutionLib/Fem/FemIVBVProblem.h"
#include "SolutionLib/Fem/TemplateFemEquation.h"
#include "SolutionLib/Fem/SingleStepFEM.h"
#include "ProcessLib/AbstractTransientProcess.h"
#include "MathLib/DataType.h"
#include "NumLib/Nonlinear/DiscreteNRSolverWithStepInitFactory.h"
#include "UserModules/FemFlash2PForm/NonLinearCMP_Flash2PForm_TimeODELocalAssembler.h"
#include "UserModules/FemFlash2PForm/NonLinearCMP_Flash2PForm_JacobianLocalAssembler.h"
#include "UserModules/FemFlash2PForm/Flash2PForm_Nested_EOS_NRIterationStepInitializer.h"
#include "Fem_CMP_Flash2P_Solution.h"
#include "SingleStepCMP_Flash2PForm.h"

#include "LocalProblem_EOS_Flash2PForm.h"
#include "EOS_Flash2PForm.h"

template <class T_DISCRETE_SYSTEM, class T_LINEAR_SOLVER >
class FunctionCMP_Flash2PForm
	: public ProcessLib::Process
{
public:
	// input variable is velocity
	enum In {};//Velocity = 0 
	// no output variable
	enum Out { Pressure = 0, Total_Molar_Fraction_1 = 0, Total_Molar_Fraction_2 = 0, Total_Molar_Fraction_3 = 0, Total_Molar_Fraction_4=0 };

	enum EOS_EVL_METHOD { 
		EOS_EVL_FIN_DIFF, 
		EOS_EVL_ANALYTICAL
	};

	// local matrix and vector
	typedef MathLib::LocalMatrix LocalMatrix;
	typedef MathLib::LocalVector LocalVector;

	typedef FunctionCMP_Flash2PForm       MyFunctionData;    // FEM problem class
	typedef T_DISCRETE_SYSTEM      MyDiscreteSystem;  // Discretization
	typedef T_LINEAR_SOLVER        MyLinearSolver;    // linear solver

	// memory for discretized concentration vector
	typedef typename FemLib::FemNodalFunctionVector<MyDiscreteSystem>::type MyNodalFunctionVector;
	typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;
	typedef typename FemLib::FemNodalFunctionMatrix<MyDiscreteSystem>::type MyNodalFunctionMatrix; 
	typedef typename FemLib::FEMIntegrationPointFunctionVector<MyDiscreteSystem>::type MyIntegrationPointFunctionVector;

	// local assembler
	typedef NonLinearCMP_Flash2PForm_TimeODELocalAssembler<NumLib::ElementWiseTimeEulerEQSLocalAssembler, MyNodalFunctionScalar, FunctionCMP_Flash2PForm>      MyNonLinearAssemblerType;
	typedef NonLinearCMP_Flash2PForm_TimeODELocalAssembler<NumLib::ElementWiseTimeEulerResidualLocalAssembler, MyNodalFunctionScalar, FunctionCMP_Flash2PForm> MyNonLinearResidualAssemblerType;
	typedef NonLinearCMP_Flash2PForm_JacobianLocalAssembler<MyNodalFunctionScalar, MyFunctionData> MyNonLinearJacobianAssemblerType;
	typedef Flash2PForm_Nested_EOS_NRIterationStepInitializer<MyNodalFunctionScalar, MyFunctionData> MyNRIterationStepInitializer;
	typedef NumLib::DiscreteNRSolverWithStepInitFactory<MyNRIterationStepInitializer> MyDiscreteNonlinearSolverFactory;

	/**
	* nonlinear coupled part solving xi
	* Equation definition
	*/
	typedef SolutionLib::TemplateFemEquation<
		MyDiscreteSystem,
		MyLinearSolver,
		MyNonLinearAssemblerType,
		MyNonLinearResidualAssemblerType,
		MyNonLinearJacobianAssemblerType
	> MyNonLinearEquationType;
	/**
	* FEM IVBV problem definition
	*/
	typedef SolutionLib::FemIVBVProblem<
		MyDiscreteSystem,
		MyNonLinearEquationType
	> MyNonLinearCMPFlash2PFormProblemType;

	/**
	* Solution algorithm definition
	*/
	typedef SolutionLib::SingleStepFEM<
		MyNonLinearCMPFlash2PFormProblemType,
		MyLinearSolver,
		MyDiscreteNonlinearSolverFactory
	> MyNonLinearSolutionType;

	/**
	  * the general CompMultiPhase problem part
	  */
	typedef SolutionLib::Fem_CMP_Flash2P_Solution<MyDiscreteSystem> MyCMPFlash2PFormProblemType;
	typedef SolutionLib::SingleStepCMP_Flash2PForm<MyFunctionData,
		MyCMPFlash2PFormProblemType,
		MyNonLinearCMPFlash2PFormProblemType,
		MyNonLinearSolutionType> MyCMPFlash2PFormSolution;
	typedef typename MyCMPFlash2PFormProblemType::MyVariable MyVariableCMPFlash2PForm;

	/**
	  * constructor
	  */
	FunctionCMP_Flash2PForm()
		: Process("CMP_Flash2PForm", 0, 8),
		_feObjects(0), _n_Comp(8), _n_Phases(2)
	{
		_EOS = new EOS_Flash2PForm();
		m_EOS_EVL_METHOD = EOS_EVL_ANALYTICAL;
	};

	/**
	  * destructor, reclaim the memory
	  */
	virtual ~FunctionCMP_Flash2PForm()
	{
		BaseLib::releaseObject(_P);
		BaseLib::releaseObject(_X1);
		BaseLib::releaseObject(_X2);
		BaseLib::releaseObject(_X3);
		BaseLib::releaseObject(_X4);
		BaseLib::releaseObject(_X5);
		BaseLib::releaseObject(_X6);
		BaseLib::releaseObject(_X7);
		
		
		BaseLib::releaseObject(_S);
		BaseLib::releaseObject(_mat_secDer); 
		BaseLib::releaseObject(_EOS);

		/*
		BaseLib::releaseObject(_problem);
		BaseLib::releaseObjectsInStdVector(_concentrations);
		*/
	};

	/**
	  * initialization of the problem class
	  */
	virtual bool initialize(const BaseLib::Options &option);

	/**
	  * finalize but nothing to do here
	  */
	virtual void finalize() {};

	/**
	  * returns the convergence checker
	  */
	virtual NumLib::IConvergenceCheck* getConvergenceChecker()
	{
		return &_checker;
	}

	/**
	  * function to solve the current time step
	  */
	virtual int solveTimeStep(const NumLib::TimeStep &time)
	{
		INFO("Solving %s...", getProcessName().c_str());
		initializeTimeStep(time);
		getSolution()->solveTimeStep(time);
		updateOutputParameter(time);
		return 0;
	}

	/**
	  * function to suggest the next time step
	  */
	virtual double suggestNext(const NumLib::TimeStep &time_current)
	{
		return getSolution()->suggestNext(time_current);
	}

	/**
	  * called when this problem is awake
	  */
	virtual bool isAwake(const NumLib::TimeStep &time)
	{
		return getSolution()->isAwake(time);
	}

	virtual bool accept(const NumLib::TimeStep &time)
	{
		return getSolution()->accept(time);
	}

	/**
	  * called when this time step is accepted
	  */
	virtual void finalizeTimeStep(const NumLib::TimeStep &time)
	{
		output(time);
		getSolution()->finalizeTimeStep(time);
	};

	/**
	  * update the value of primary variable P
	  */
    void set_P(MyNodalFunctionScalar* new_nodal_values)
	{
		std::size_t node_idx;
        double nodal_P;

        for (node_idx = _P->getDiscreteData()->getRangeBegin();
             node_idx < _P->getDiscreteData()->getRangeEnd();
			 node_idx++)
		{
            nodal_P = new_nodal_values->getValue(node_idx);
            this->_P->setValue(node_idx, nodal_P);
		}
	};

    /**
      * update the value of primary variable X
      */
    void set_X1(MyNodalFunctionScalar* new_nodal_values)
    {
        std::size_t node_idx;
        double nodal_X1;

        for (node_idx = _X1->getDiscreteData()->getRangeBegin();
             node_idx < _X1->getDiscreteData()->getRangeEnd();
             node_idx++)
        {
            nodal_X1 = new_nodal_values->getValue(node_idx);
            this->_X1->setValue(node_idx, nodal_X1);
        }
    };
	/**
	* update the value of primary variable X
	*/
	void set_X2(MyNodalFunctionScalar* new_nodal_values)
	{
		std::size_t node_idx;
		double nodal_X2;

		for (node_idx = _X2->getDiscreteData()->getRangeBegin();
			node_idx < _X2->getDiscreteData()->getRangeEnd();
			node_idx++)
		{
			nodal_X2 = new_nodal_values->getValue(node_idx);
			this->_X2->setValue(node_idx, nodal_X2);
		}
	};
	/**
	* update the value of primary variable X
	*/
	void set_X3(MyNodalFunctionScalar* new_nodal_values)
	{
		std::size_t node_idx;
		double nodal_X3;

		for (node_idx = _X3->getDiscreteData()->getRangeBegin();
			node_idx < _X3->getDiscreteData()->getRangeEnd();
			node_idx++)
		{
			nodal_X3 = new_nodal_values->getValue(node_idx);
			this->_X3->setValue(node_idx, nodal_X3);
		}
	};
	/**
	* update the value of primary variable X
	*/
	void set_X4(MyNodalFunctionScalar* new_nodal_values)
	{
		std::size_t node_idx;
		double nodal_X4;

		for (node_idx = _X4->getDiscreteData()->getRangeBegin();
			node_idx < _X4->getDiscreteData()->getRangeEnd();
			node_idx++)
		{
			nodal_X4 = new_nodal_values->getValue(node_idx);
			this->_X4->setValue(node_idx, nodal_X4);
		}
	};

	/**
	* update the value of primary variable X
	*/
	void set_X5(MyNodalFunctionScalar* new_nodal_values)
	{
		std::size_t node_idx;
		double nodal_X5;

		for (node_idx = _X5->getDiscreteData()->getRangeBegin();
			node_idx < _X5->getDiscreteData()->getRangeEnd();
			node_idx++)
		{
			nodal_X5 = new_nodal_values->getValue(node_idx);
			this->_X5->setValue(node_idx, nodal_X5);
		}
	};
	/**
	* update the value of primary variable X
	*/
	void set_X6(MyNodalFunctionScalar* new_nodal_values)
	{
		std::size_t node_idx;
		double nodal_X6;

		for (node_idx = _X6->getDiscreteData()->getRangeBegin();
			node_idx < _X6->getDiscreteData()->getRangeEnd();
			node_idx++)
		{
			nodal_X6 = new_nodal_values->getValue(node_idx);
			this->_X6->setValue(node_idx, nodal_X6);
		}
	};
	/**
	* update the value of primary variable X
	*/
	void set_X7(MyNodalFunctionScalar* new_nodal_values)
	{
		std::size_t node_idx;
		double nodal_X7;

		for (node_idx = _X7->getDiscreteData()->getRangeBegin();
			node_idx < _X7->getDiscreteData()->getRangeEnd();
			node_idx++)
		{
			nodal_X7 = new_nodal_values->getValue(node_idx);
			this->_X7->setValue(node_idx, nodal_X7);
		}
	};
	/*
	MyNodalFunctionScalar* get_concentrations(size_t idx_conc)
	{
		return _concentrations[idx_conc];
	};
	*/

	/**
	  * set mean pressure nodal value
	  */
	void set_P_node_values(std::size_t node_idx, double node_value){ _P->setValue(node_idx, node_value);  };

	/**
	  * set molar fraction nodal value
	  */
	void set_X1_node_values(std::size_t node_idx, double node_value){ _X1->setValue(node_idx, node_value); };
	void set_X2_node_values(std::size_t node_idx, double node_value){ _X2->setValue(node_idx, node_value); };
	void set_X3_node_values(std::size_t node_idx, double node_value){ _X3->setValue(node_idx, node_value); };
	void set_X4_node_values(std::size_t node_idx, double node_value){ _X4->setValue(node_idx, node_value); };
	void set_X5_node_values(std::size_t node_idx, double node_value){ _X5->setValue(node_idx, node_value); };
	void set_X6_node_values(std::size_t node_idx, double node_value){ _X6->setValue(node_idx, node_value); };
	void set_X7_node_values(std::size_t node_idx, double node_value){ _X7->setValue(node_idx, node_value); };
	/**
	  * calculate nodal Equation-of-State (EOS) system
	  */
	void calc_nodal_eos_sys(double dt);
	/**
	* Here define two functions 
	* return the matrix of M and K
	* would be helpful for constructing the Jacobian matrix
	*/
	std::vector<LocalMatrix>& get_elem_M_matrix(void) { return _elem_M_matrix; };
	std::vector<LocalMatrix>& get_elem_K_matrix(void) { return _elem_K_matrix; };

	/**
	  * return the number of chemical components
	  */
	std::size_t get_n_Comp(void) { return _n_Comp; };


	/**
	* return the value of gas phase saturation
	*/
	MyIntegrationPointFunctionVector* getS(void) { return _S; };
	/**
	* return the value of capillary pressure
	*/
	MyIntegrationPointFunctionVector* get_PC(void) { return _PC; };
	/**
	* return the value of Liquid  pressure
	*/
	MyIntegrationPointFunctionVector* get_PL(void) { return _PL; };
	/**
	* return the value of Gas phase pressure
	*/
	MyIntegrationPointFunctionVector* get_PG(void) { return _PG; };
	/**
	* return the dissolved gas mass density
	*/
	MyIntegrationPointFunctionVector* get_rhoLh(void) { return _rho_L_h; };
	/**
	* return the  gas mass density IN GAS PHASE
	*/
	MyIntegrationPointFunctionVector* get_rhoGh(void) { return _rho_G_h; };
	/**
	* return the  water vapor mass density IN GAS PHASE
	*/
	MyIntegrationPointFunctionVector* get_rhoGw(void) { return _rho_G_w; };

	MyIntegrationPointFunctionVector* get_dPGh_dPg(void) { return _dPGh_dPg; };
	MyIntegrationPointFunctionVector* get_dPc_dSg(void) { return _dPcdSg; };

	MyIntegrationPointFunctionVector* get_Func_C(void)  { return _Func_C; };
	/**
	* return the matrix of derivative of secondary variable based on P and X
	*/
	MyNodalFunctionMatrix* get_mat_secDer(void) { return _mat_secDer; };
	MyNodalFunctionVector* get_vec_tempVar(void) { return _vec_tempVar; };
protected:
	virtual void initializeTimeStep(const NumLib::TimeStep &time);

	/**
	  * this function is called to exchange output parameters
	  */
	virtual void updateOutputParameter(const NumLib::TimeStep &time);

	/**
	  * get the pointer of solution class for current problem
	  */
	virtual MyCMPFlash2PFormSolution* getSolution() { return _solution; };

	/**
	  * output the result of current solution
	  */
	virtual void output(const NumLib::TimeStep &time);

private:
	DISALLOW_COPY_AND_ASSIGN(FunctionCMP_Flash2PForm);

private:
	/**
	* Nonlinear iterator
	*/
	MyNRIterationStepInitializer*              myNRIterator;

	/**
	* Nonlinear solver factory
	*/
	MyDiscreteNonlinearSolverFactory*          myNSolverFactory;

	/**
	* nonlinear equation
	*/
	MyNonLinearEquationType*                  _non_linear_eqs;

	/**
	* non-linear problem
	*/
	MyNonLinearCMPFlash2PFormProblemType*     _non_linear_problem;

	/**
	* non-linear solution
	*/
	MyNonLinearSolutionType*                  _non_linear_solution; 

	/**
	* degree of freedom equation ID talbe for the nonlinear problem
	*/
	DiscreteLib::DofEquationIdTable * _nl_sol_dofManager;

	/**
	  * Component based multiphase flow problem
	  */
	MyCMPFlash2PFormProblemType* _problem;

	/**
	  * solution class for the component based multiphase flow problem
	  */
	MyCMPFlash2PFormSolution*    _solution;

	/**
	  * FEM object
	  */
	FemLib::LagrangeFeObjectContainer* _feObjects;

	/**
	  * convergence checker
	  */
	NumLib::DiscreteDataConvergenceCheck _checker;

	/**
	  * Primary variable 1): vector of mean pressure values on each node
	  */
	MyNodalFunctionScalar* _P;

	/**
	  * Primary variable 2): vector of molar fraction values of the light component on each node
	  */
	MyNodalFunctionScalar* _X1;
	MyNodalFunctionScalar* _X2;
	MyNodalFunctionScalar* _X3;
	MyNodalFunctionScalar* _X4;
	MyNodalFunctionScalar* _X5;
	MyNodalFunctionScalar* _X6;
	MyNodalFunctionScalar* _X7;
	/**
	*		Vector _output
	*		store the eight values of the secondary variables
	*/
	
	/** 
	  * secondary variable --saturation
	  */

	MyIntegrationPointFunctionVector* _S;
	MyIntegrationPointFunctionVector* _PG;
	MyIntegrationPointFunctionVector* _PL;
	MyIntegrationPointFunctionVector* _PC;
	MyIntegrationPointFunctionVector* _rho_L;
	MyIntegrationPointFunctionVector* _rho_G;
	
	MyIntegrationPointFunctionVector* _dPGh_dPg;
	MyIntegrationPointFunctionVector* _dPcdSg;//
	MyIntegrationPointFunctionVector* _Func_C;//

	/**
	* store the local matrix M and K 
	*/
	std::vector<LocalMatrix> _elem_M_matrix; 
	std::vector<LocalMatrix> _elem_K_matrix;
		
	/**
	  * derivative 
	  * for each node 8 by 2 matrix. 
	  */
	MyNodalFunctionMatrix* _mat_secDer;
	/**
	* define a vectior to store the value of Omega M and Characteristic function
	* therefore the size should be 3.
	*/
	MyNodalFunctionVector* _vec_tempVar;

	LocalProblem_EOS_Flash2PForm* _LP_EOS;//LOCAL PROBLEM
	/**
	  * the id of the msh applied in this process
	  */
	size_t _msh_id;

	/**
	  * number of chemical components, 
	  * in the current UserModule, it is fixed to two. 
	  */
	const size_t _n_Comp;

	/**
	  * number of fluid phases,
	  * in the current UserModule, it is also fixed to two.
	  */
	const size_t _n_Phases;

	EOS_Flash2PForm* _EOS;
	EOS_EVL_METHOD m_EOS_EVL_METHOD;
	double real_x[3], gp_x[3];
};

#include "FunctionCMP_Flash2PForm.hpp"

#endif  // end of ifndef
