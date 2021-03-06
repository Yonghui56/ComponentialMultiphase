/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FEMEnums.cpp
 *
 * Created on 2010-09-02 by Thomas Fischer
 */

#include "FEMEnums.h"
#include <cstdlib>
#include <iostream>

namespace ogs5
{

namespace FiniteElement
{
ProcessType convertProcessType ( const std::string& pcs_type_string )
{
    if (pcs_type_string.compare ("LIQUID_FLOW") == 0)
        return LIQUID_FLOW;
    if (pcs_type_string.compare ("FLUID_FLOW") == 0)
        return FLUID_FLOW;
    if (pcs_type_string.compare ("TWO_PHASE_FLOW") == 0)
        return TWO_PHASE_FLOW;
    if (pcs_type_string.compare ("RICHARDS_FLOW") == 0)
        return RICHARDS_FLOW;
    if (pcs_type_string.compare ("OVERLAND_FLOW") == 0)
        return OVERLAND_FLOW;
    if (pcs_type_string.compare ("GROUNDWATER_FLOW") == 0)
        return GROUNDWATER_FLOW;
    if (pcs_type_string.compare ("HEAT_TRANSPORT") == 0)
        return HEAT_TRANSPORT;
    if (pcs_type_string.compare ("DEFORMATION") == 0)
        return DEFORMATION;
    if (pcs_type_string.compare ("DEFORMATION_FLOW") == 0)
        return DEFORMATION_FLOW;
    if (pcs_type_string.compare ("DEFORMATION_DYNAMIC") == 0)
        return DEFORMATION_DYNAMIC;
    if (pcs_type_string.compare ("MASS_TRANSPORT") == 0)
        return MASS_TRANSPORT;
	if (pcs_type_string.compare ("KIN_REACT_GIA") == 0)  // HS 09.2012
        return KIN_REACT_GIA;
	if (pcs_type_string.compare ("REACT_GIA") == 0)  // RZ 05.2013
        return REACT_GIA;
    if (pcs_type_string.compare ("REACT_TRANS_OPS") == 0)  // RZ 05.2013
        return REACT_TRANS_OPS;
	if (pcs_type_string.compare ("MULTI_PHASE_FLOW") == 0)
        return MULTI_PHASE_FLOW;
    if (pcs_type_string.compare ("DEFORMATION_H2") == 0)
        return DEFORMATION_H2;
    if (pcs_type_string.compare ("AIR_FLOW") == 0)
        return AIR_FLOW;
    if (pcs_type_string.compare ("FLUID_MOMENTUM") == 0)
        return FLUID_MOMENTUM;
    if (pcs_type_string.compare ("RANDOM_WALK") == 0)
        return RANDOM_WALK;
    if (pcs_type_string.compare ("FLUX") == 0)
        return FLUX;
    if (pcs_type_string.compare ("PS_GLOBAL") == 0)
        return PS_GLOBAL;
    if (pcs_type_string.compare ("NO_PCS") == 0)
        return NO_PCS;
    if (pcs_type_string.compare ("PTC_FLOW") == 0)
        return PTC_FLOW;
	if (pcs_type_string.compare("CMP_2P2C") == 0)
		return CMP_2P2C;
	if (pcs_type_string.compare("CMP_PressureForm") == 0)
		return CMP_PressureForm;
	if (pcs_type_string.compare("CMP_TotalDensityForm") == 0)
		return CMP_TotalDensityForm;
	if (pcs_type_string.compare("CMP_GlobalComplementaryForm") == 0)
		return CMP_GlobalComplementaryForm;
	if (pcs_type_string.compare("CMP_CapPresFormular") == 0)
		return CMP_CapPresFormular;
    //else
        //std::cout << "WARNING in convertProcessType: process type #" << pcs_type_string <<
        //"# unknown" << std::endl;
    return INVALID_PROCESS;
}

std::string convertProcessTypeToString ( ProcessType pcs_type )
{
    if (pcs_type == LIQUID_FLOW)
        return "LIQUID_FLOW";
    if (pcs_type == FLUID_FLOW)
        return "FLUID_FLOW";
    if (pcs_type == TWO_PHASE_FLOW)
        return "TWO_PHASE_FLOW";
    if (pcs_type == RICHARDS_FLOW)
        return "RICHARDS_FLOW";
    if (pcs_type == OVERLAND_FLOW)
        return "OVERLAND_FLOW";
    if (pcs_type == GROUNDWATER_FLOW)
        return "GROUNDWATER_FLOW";
    if (pcs_type == HEAT_TRANSPORT)
        return "HEAT_TRANSPORT";
    if (pcs_type == DEFORMATION)
        return "DEFORMATION";
    if (pcs_type == DEFORMATION_FLOW)
        return "DEFORMATION_FLOW";
    if (pcs_type == DEFORMATION_DYNAMIC)
        return "DEFORMATION_DYNAMIC";
    if (pcs_type == MASS_TRANSPORT)
        return "MASS_TRANSPORT";
	if (pcs_type == KIN_REACT_GIA)  // HS 09.2012
		return "KIN_REACT_GIA";
    if (pcs_type == REACT_TRANS_OPS)  // HS 09.2012
		return "REACT_TRANS_OPS";
	if (pcs_type == REACT_GIA)  // RZ 05.2013
		return "REACT_GIA";
    if (pcs_type == MULTI_PHASE_FLOW)
        return "MULTI_PHASE_FLOW";
    if (pcs_type == DEFORMATION_H2)
        return "DEFORMATION_H2";
    if (pcs_type == AIR_FLOW)
        return "AIR_FLOW";
    if (pcs_type == FLUID_MOMENTUM)
        return "FLUID_MOMENTUM";
    if (pcs_type == RANDOM_WALK)
        return "RANDOM_WALK";
    if (pcs_type == FLUX)
        return "FLUX";
    if (pcs_type ==   PS_GLOBAL)
        return "PS_GLOBAL";
    if (pcs_type == PTC_FLOW)
        return "PTC_FLOW";
	if (pcs_type == CMP_2P2C)
		return "CMP_2P2C";
	if (pcs_type == CMP_PressureForm)
		return "CMP_PressureForm";
	if (pcs_type == CMP_TotalDensityForm)
		return "CMP_TotalDensityForm";
	if (pcs_type == CMP_GlobalComplementaryForm)
		return "CMP_GlobalComplementaryForm";
	if (pcs_type == CMP_CapPresFormular)
		return "CMP_CapPresFormular";
	if (pcs_type ==   NO_PCS)
        return "NO_PCS";

    return "INVALID_PROCESS";
}


const std::list<std::string> getAllProcessNames()
{
    size_t count(1);
    std::list<std::string> enum_names;

    while (count != PROCESS_END)
    {
        enum_names.push_back( convertProcessTypeToString(static_cast<ProcessType>(count++)) );
    }
    return enum_names;
}


PrimaryVariable convertPrimaryVariable ( const std::string& pcs_pv_string )
{
    if (pcs_pv_string.compare ("PRESSURE1") == 0)
        return PRESSURE;
    if (pcs_pv_string.compare ("PRESSURE2") == 0)
        return PRESSURE2;
    if (pcs_pv_string.compare ("PRESSURE_RATE1") == 0)
        return PRESSURE_RATE1;
    if (pcs_pv_string.compare ("SATURATION1") == 0)
        return SATURATION;
    if (pcs_pv_string.compare ("SATURATION2") == 0)
        return SATURATION2;
    if (pcs_pv_string.compare ("TEMPERATURE1") == 0)
        return TEMPERATURE;
    if (pcs_pv_string.compare ("DISPLACEMENT_X1") == 0)
        return DISPLACEMENT_X;
    if (pcs_pv_string.compare ("DISPLACEMENT_Y1") == 0)
        return DISPLACEMENT_Y;
    if (pcs_pv_string.compare ("DISPLACEMENT_Z1") == 0)
        return DISPLACEMENT_Z;
    if (pcs_pv_string.compare ("CONCENTRATION1") == 0)
        return CONCENTRATION;
    if (pcs_pv_string.compare ("HEAD") == 0)
        return HEAD;
    if (pcs_pv_string.compare ("VELOCITY_DM_X") == 0)
        return VELOCITY_DM_X;
    if (pcs_pv_string.compare ("VELOCITY_DM_Y") == 0)
        return VELOCITY_DM_Y;
    if (pcs_pv_string.compare ("VELOCITY_DM_Z") == 0)
        return VELOCITY_DM_Z;
    if (pcs_pv_string.compare ("VELOCITY1_X") == 0)
        return VELOCITY1_X;
    if (pcs_pv_string.compare ("VELOCITY1_Y") == 0)
        return VELOCITY1_Y;
    if (pcs_pv_string.compare ("VELOCITY1_Z") == 0)
        return VELOCITY1_Z;
    if (pcs_pv_string.compare ("STRESS_XX") == 0)
        return STRESS_XX;
    if (pcs_pv_string.compare ("STRESS_XY") == 0)
        return STRESS_XY;
    if (pcs_pv_string.compare ("STRESS_XZ") == 0)
        return STRESS_XZ;
    if (pcs_pv_string.compare ("STRESS_YY") == 0)
        return STRESS_YY;
    if (pcs_pv_string.compare ("STRESS_YZ") == 0)
        return STRESS_YZ;
    if (pcs_pv_string.compare ("STRESS_ZZ") == 0)
        return STRESS_ZZ;
    if (pcs_pv_string.compare ("ACCELERATION_X1") == 0)
        return ACCELERATION_X1;
    if (pcs_pv_string.compare ("ACCELERATION_Y1") == 0)
        return ACCELERATION_Y1;
    if (pcs_pv_string.compare ("ACCELERATION_Z1") == 0)
        return ACCELERATION_Z1;
    if (pcs_pv_string.compare ("EXCAVATION") == 0)
        return EXCAVATION;
    if (pcs_pv_string.compare ("STRAIN_XX") == 0)
        return STRAIN_XX;
    if (pcs_pv_string.compare ("STRAIN_XY") == 0)
        return STRAIN_XY;
    if (pcs_pv_string.compare ("STRAIN_XZ") == 0)
        return STRAIN_XZ;
    if (pcs_pv_string.compare ("STRAIN_YY") == 0)
        return STRAIN_YY;
    if (pcs_pv_string.compare ("STRAIN_YZ") == 0)
        return STRAIN_YZ;
    if (pcs_pv_string.compare ("STRAIN_ZZ") == 0)
        return STRAIN_ZZ;
    if (pcs_pv_string.compare ("STRAIN_PLS") == 0)
        return STRAIN_PLS;
    //else
    //{
        //std::cout << "convertPrimaryVariable #" << pcs_pv_string << "# not found" << std::endl;
        //exit (1);
    //}
    return INVALID_PV;
}

std::string convertPrimaryVariableToString ( PrimaryVariable pcs_pv )
{
    if (pcs_pv == PRESSURE)
        return "PRESSURE1";
    if (pcs_pv == PRESSURE2)
        return "PRESSURE2";
    if (pcs_pv == PRESSURE_RATE1)
        return "PRESSURE_RATE1";
    if (pcs_pv == SATURATION)
        return "SATURATION1";
    if (pcs_pv == SATURATION2)
        return "SATURATION2";
    if (pcs_pv == TEMPERATURE)
        return "TEMPERATURE1";
    if (pcs_pv == DISPLACEMENT_X)
        return "DISPLACEMENT_X1";
    if (pcs_pv == DISPLACEMENT_Y)
        return "DISPLACEMENT_Y1";
    if (pcs_pv == DISPLACEMENT_Z)
        return "DISPLACEMENT_Z1";
    if (pcs_pv == CONCENTRATION)
        return "CONCENTRATION1";
    if (pcs_pv == HEAD)
        return "HEAD";
    if (pcs_pv == VELOCITY_DM_X)
        return "VELOCITY_DM_X";
    if (pcs_pv == VELOCITY_DM_Y)
        return "VELOCITY_DM_Y";
    if (pcs_pv == VELOCITY_DM_Z)
        return "VELOCITY_DM_Z";
    if (pcs_pv == VELOCITY1_X)
        return "VELOCITY1_X";
    if (pcs_pv == VELOCITY1_Y)
        return "VELOCITY1_Y";
    if (pcs_pv == VELOCITY1_Z)
        return "VELOCITY1_Z";
    if (pcs_pv == STRESS_XX)
        return "STRESS_XX";
    if (pcs_pv == STRESS_XY)
        return "STRESS_XY";
    if (pcs_pv == STRESS_XZ)
        return "STRESS_XZ";
    if (pcs_pv == STRESS_YY)
        return "STRESS_YY";
    if (pcs_pv == STRESS_YZ)
        return "STRESS_YZ";
    if (pcs_pv == STRESS_ZZ)
        return "STRESS_ZZ";
    if (pcs_pv == STRAIN_XX) return "STRAIN_XX";
    if (pcs_pv == STRAIN_XY) return "STRAIN_XY";
    if (pcs_pv == STRAIN_XZ) return "STRAIN_XZ";
    if (pcs_pv == STRAIN_YY) return "STRAIN_YY";
    if (pcs_pv == STRAIN_YZ) return "STRAIN_YZ";
    if (pcs_pv == STRAIN_ZZ) return "STRAIN_ZZ";
    if (pcs_pv == STRAIN_PLS) return "STRAIN_PLS";
    if (pcs_pv == ACCELERATION_X1)
        return "ACCELERATION_X1";
    if (pcs_pv == ACCELERATION_Y1)
        return "ACCELERATION_Y1";
    if (pcs_pv == ACCELERATION_Z1)
        return "ACCELERATION_Z1";
    if (pcs_pv == EXCAVATION)
        return "EXCAVATION";
    return "INVALID_PRIMARY_VARIABLE";
}

const std::list<std::string> getAllPrimaryVariableNames()
{
    size_t count(1);
    std::list<std::string> enum_names;

    while (count != PV_END)
    {
        enum_names.push_back( convertPrimaryVariableToString(static_cast<PrimaryVariable>(count++)) );
    }
    return enum_names;
}

DistributionType convertDisType(const std::string& dis_type_string)
{
    if (dis_type_string.compare("CONSTANT") == 0)
        return CONSTANT;
    if (dis_type_string.compare("ANALYTICAL") == 0)
        return ANALYTICAL;
    if (dis_type_string.compare("AVERAGE") == 0)
        return AVERAGE;
    if (dis_type_string.compare("CONSTANT_GEO") == 0)
        return CONSTANT_GEO;
    if (dis_type_string.compare("GRADIENT") == 0)
        return GRADIENT;
    if (dis_type_string.compare("RESTART") == 0)
        return RESTART;
    if (dis_type_string.compare("LINEAR") == 0)
        return LINEAR;
    if (dis_type_string.compare("POINT") == 0)
        return POINT;
    if (dis_type_string.compare("CONSTANT_NEUMANN") == 0)
        return CONSTANT_NEUMANN;
    if (dis_type_string.compare("LINEAR_NEUMANN") == 0)
        return LINEAR_NEUMANN;
    if (dis_type_string.compare("NORMALDEPTH") == 0)
        return NORMALDEPTH;
    if (dis_type_string.compare("CRITICALDEPTH") == 0)
        return CRITICALDEPTH;
    if (dis_type_string.compare("GREEN_AMPT") == 0)
        return GREEN_AMPT;
    if (dis_type_string.compare("SYSTEM_DEPENDENT") == 0)
        return SYSTEM_DEPENDENT;
    if (dis_type_string.compare("PRECIPITATION") == 0)
        return PRECIPITATION;
    if (dis_type_string.compare("DIRECT") == 0)
        return DIRECT;
    if (dis_type_string.compare("FUNCTION") == 0)
        return FUNCTION;                              //24.08.2011. WW
    else
    {
        std::cout << "convertDisType #" << dis_type_string << "# not found"
                  << std::endl;
        exit(1);
    }
    return INVALID_DIS_TYPE;
}

std::string convertDisTypeToString(DistributionType dis_type)
{
    if (dis_type == ANALYTICAL)
        return "ANALYTICAL";
    if (dis_type == AVERAGE)
        return "AVERAGE";
    if (dis_type == CONSTANT)
        return "CONSTANT";
    if (dis_type == CONSTANT_GEO)
        return "CONSTANT_GEO";
    if (dis_type == GRADIENT)
        return "GRADIENT";
    if (dis_type == RESTART)
        return "RESTART";
    if (dis_type == LINEAR)
        return "LINEAR";
    if (dis_type == POINT)
        return "POINT";
    if (dis_type == CONSTANT_NEUMANN)
        return "CONSTANT_NEUMANN";
    if (dis_type == LINEAR_NEUMANN)
        return "LINEAR_NEUMANN";
    if (dis_type == NORMALDEPTH)
        return "NORMALDEPTH";
    if (dis_type == CRITICALDEPTH)
        return "CRITICALDEPTH";
    if (dis_type == GREEN_AMPT)
        return "GREEN_AMPT";
    if (dis_type == SYSTEM_DEPENDENT)
        return "SYSTEM_DEPENDENT";
    if (dis_type == PRECIPITATION)
        return "PRECIPITATION";
    if (dis_type == DIRECT)
        return "DIRECT";
    if (dis_type == FUNCTION)
        return "FUNCTION";         //24.08.2011. WW

    return "INVALID_DIS_TYPE";
}

const std::list<std::string> getAllDistributionNames()
{
    size_t count(1);
    std::list<std::string> enum_names;

    while (count != DIS_END)
    {
        enum_names.push_back( convertDisTypeToString(static_cast<DistributionType>(count++)) );
    }
    return enum_names;
}

ErrorMethod convertErrorMethod(const std::string& error_method_string)
{
    if (error_method_string.compare("LMAX") == 0)
        return LMAX;
    if (error_method_string.compare("ENORM") == 0)
        return ENORM;
    if (error_method_string.compare("EVNORM") == 0)
        return EVNORM;
    if (error_method_string.compare("ERNORM") == 0)
        return ERNORM;
    if (error_method_string.compare("BNORM") == 0)
        return BNORM;
    else
    {
        std::cout << "convertErrorMethod #" << error_method_string << "# not found"<< std::endl;
        exit(1);
    }
    return INVALID_ERROR_METHOD;
}

} // end namespace FiniteElement
}
