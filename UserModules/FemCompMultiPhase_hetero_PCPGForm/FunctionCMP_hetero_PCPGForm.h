/**
* Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file FunctionCMP_PressureForm.h
*
* Created on    2015-03-12 by Yonghui HUANG
*/

#ifndef FUNCTIONCMP_HETERO_PCPGFORM_H
#define FUNCTIONCMP_HETERO_PCPGFORM_H

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
#include "UserModules/FemCompMultiPhase_hetero_PCPGForm/NonLinearCMP_hetero_PCPGForm_TimeODELocalAssembler.h"
#include "UserModules/FemCompMultiPhase_hetero_PCPGForm/NonLinearCMP_hetero_PCPGForm_JacobianLocalAssembler.h"
#include "UserModules/FemCompMultiPhase_hetero_PCPGForm/Hetero_PCPGForm_Nested_EOS_NRIterationStepInitializer.h"
#include "Fem_CMP_hetero_PCPGForm_Solution.h"
#include "SingleStepCMP_hetero_PCPGForm.h"
#include "EOS_hetero_PCPGForm.h"

template <class T_DISCRETE_SYSTEM, class T_LINEAR_SOLVER >
class FunctionCMP_hetero_PCPGForm
	: public ProcessLib::Process
{
public:
	// input variable is velocity
	enum In {};//Velocity = 0 
	// no output variable
	enum Out { Gas_Pressure = 0 , Capillary_Pressure = 0};

	enum EOS_EVL_METHOD { 
		EOS_EVL_FIN_DIFF, 
		EOS_EVL_ANALYTICAL
	};

	// local matrix and vector
	typedef MathLib::LocalMatrix LocalMatrix;
	typedef MathLib::LocalVector LocalVector;

	typedef FunctionCMP_hetero_PCPGForm       MyFunctionData;    // FEM problem class
	typedef T_DISCRETE_SYSTEM      MyDiscreteSystem;  // Discretization
	typedef T_LINEAR_SOLVER        MyLinearSolver;    // linear solver

	// memory for discretized concentration vector
	typedef typename FemLib::FemNodalFunctionVector<MyDiscreteSystem>::type MyNodalFunctionVector;
	typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;
	typedef typename FemLib::FemNodalFunctionMatrix<MyDiscreteSystem>::type MyNodalFunctionMatrix; 
	typedef typename FemLib::FEMIntegrationPointFunctionVector<MyDiscreteSystem>::type MyIntegrationPointFunctionVector;
	// local assembler
	typedef NonLinearCMP_hetero_PCPGForm_TimeODELocalAssembler<NumLib::ElementWiseTimeEulerEQSLocalAssembler, MyNodalFunctionScalar, FunctionCMP_hetero_PCPGForm>      MyNonLinearAssemblerType;
	typedef NonLinearCMP_hetero_PCPGForm_TimeODELocalAssembler<NumLib::ElementWiseTimeEulerResidualLocalAssembler, MyNodalFunctionScalar, FunctionCMP_hetero_PCPGForm> MyNonLinearResidualAssemblerType;
	typedef NonLinearCMP_hetero_PCPGForm_JacobianLocalAssembler<MyNodalFunctionScalar, MyFunctionData> MyNonLinearJacobianAssemblerType; 
	typedef Hetero_PCPGForm_Nested_EOS_NRIterationStepInitializer<MyNodalFunctionScalar, MyFunctionData> MyNRIterationStepInitializer;
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
	> MyNonLinearCMPPCPGFORMProblemType;

	/**
	* Solution algorithm definition
	*/
	typedef SolutionLib::SingleStepFEM<
		MyNonLinearCMPPCPGFORMProblemType,
		MyLinearSolver,
		MyDiscreteNonlinearSolverFactory
	> MyNonLinearSolutionType;

	/**
	  * the general CompMultiPhase problem part
	  */

	typedef SolutionLib::Fem_CMP_hetero_PCPGForm_Solution<MyDiscreteSystem> MyCMPPCPGFormProblemType;
	typedef SolutionLib::SingleStepCMP_hetero_PCPGForm<MyFunctionData,
		MyCMPPCPGFormProblemType,
		MyNonLinearCMPPCPGFORMProblemType,
		MyNonLinearSolutionType> MyCMPPCPGFormSolution;
	typedef typename MyCMPPCPGFormProblemType::MyVariable MyVariableCMPPCPGForm;

	/**
	  * constructor
	  */
	FunctionCMP_hetero_PCPGForm()
		: Process("CMP_hetero_PCPGForm", 0, 2),
		_feObjects(0), _n_Comp(2), _n_Phases(2)
	{
		_EOS = new EOS_hetero_PCPGForm();
	};

	/**
	  * destructor, reclaim the memory
	  */
	virtual ~FunctionCMP_hetero_PCPGForm()
	{
		BaseLib::releaseObject(_PG);
		BaseLib::releaseObject(_PC);

		
		BaseLib::releaseObject(_S);
		BaseLib::releaseObject(_dPc);
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
      * update the value of primary variable PG
      */
    void set_PG(MyNodalFunctionScalar* new_nodal_values)
    {
        std::size_t node_idx;
        double nodal_PG;

        for (node_idx = _PG->getDiscreteData()->getRangeBegin();
             node_idx < _PG->getDiscreteData()->getRangeEnd();
             node_idx++)
        {
            nodal_PG = new_nodal_values->getValue(node_idx);
            this->_PG->setValue(node_idx, nodal_PG);
        }
    };
	/**
	  * update the value of primary variable PC
	  */
    void set_PC(MyNodalFunctionScalar* new_nodal_values)
	{
		std::size_t node_idx;
        double nodal_PC;

        for (node_idx = _PC->getDiscreteData()->getRangeBegin();
             node_idx < _PC->getDiscreteData()->getRangeEnd();
			 node_idx++)
		{
            nodal_PC = new_nodal_values->getValue(node_idx);
            this->_PC->setValue(node_idx, nodal_PC);
		}
	};

   

	/*
	MyNodalFunctionScalar* get_concentrations(size_t idx_conc)
	{
		return _concentrations[idx_conc];
	};
	*/


	/**
	  * set molar fraction nodal value
	  */
	void set_PG_node_values(std::size_t node_idx, double node_value){ _PG->setValue(node_idx, node_value); };
	/**
	  * set mean pressure nodal value
	  */
	void set_PC_node_values(std::size_t node_idx, double node_value){ _PC->setValue(node_idx, node_value);  };

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

	MyIntegrationPointFunctionVector* getdPC(void) { return _dPc; };
	/**
	* return the matrix of derivative of secondary variable based on P and X
	*/
	MyNodalFunctionMatrix* get_mat_secDer(void) { return _mat_secDer; };

protected:
	virtual void initializeTimeStep(const NumLib::TimeStep &time);

	/**
	  * this function is called to exchange output parameters
	  */
	virtual void updateOutputParameter(const NumLib::TimeStep &time);

	/**
	  * get the pointer of solution class for current problem
	  */
	virtual MyCMPPCPGFormSolution* getSolution() { return _solution; };

	/**
	  * output the result of current solution
	  */
	virtual void output(const NumLib::TimeStep &time);

private:
	DISALLOW_COPY_AND_ASSIGN(FunctionCMP_hetero_PCPGForm);

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
	MyNonLinearCMPPCPGFORMProblemType*     _non_linear_problem;

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
	MyCMPPCPGFormProblemType* _problem;

	/**
	  * solution class for the component based multiphase flow problem
	  */
	MyCMPPCPGFormSolution*    _solution;

	/**
	  * FEM object
	  */
	FemLib::LagrangeFeObjectContainer* _feObjects;

	/**
	  * convergence checker
	  */
	NumLib::DiscreteDataConvergenceCheck _checker;

	/**
	  * Primary variable 1): vector of CAPILLARY pressure values on each node
	  */
	MyNodalFunctionScalar* _PG;

	/**
	  * Primary variable 2): vector of gas pressure values of the light component on each node
	  */
	MyNodalFunctionScalar* _PC;
	/**
	*		Vector _output
	*		store the eight values of the secondary variables
	*/
	
	/** 
	  * secondary variable --saturation
	  */
	MyIntegrationPointFunctionVector* _S;
	MyIntegrationPointFunctionVector* _dPc;
	//MyNodalFunctionScalar* _S;
	//MyNodalFunctionScalar* _dPc;
	/**
	* store the local matrix M and K 
	*/
	std::vector<LocalMatrix> _elem_M_matrix; 
	std::vector<LocalMatrix> _elem_K_matrix;
		
	/**
	  * derivative of PC
	  * for each node 8 by 2 matrix. 
	  */
	MyNodalFunctionMatrix* _mat_secDer;

	
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

	EOS_hetero_PCPGForm* _EOS;


};

#include "FunctionCMP_hetero_PCPGForm.hpp"

#endif  // end of ifndef
