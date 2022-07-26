/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2022 University of Utah, The Trustees of Columbia University in
 the City of New York, and others.
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.*/



#include "stdafx.h"
#include "FEBioPolarFluid.h"
#include <FECore/FECoreKernel.h>
#include "FEPolarFluidSolver.h"
#include "FEPolarFluid.h"
#include "FEPolarFluidDomain3D.h"
#include "FEPolarFluidDomainFactory.h"
#include "FEPolarFluidAnalysis.h"
#include "FETangentialFlowPFStabilization.h"
#include "FEFluidModule.h"

//-----------------------------------------------------------------------------
const char* FEBioPolarFluid::GetVariableName(FEBioPolarFluid::POLAR_FLUID_VARIABLE var)
{
    switch (var)
    {
        case DISPLACEMENT                : return "displacement"               ; break;
        case RELATIVE_FLUID_VELOCITY     : return "relative fluid velocity"    ; break;
        case RELATIVE_FLUID_ACCELERATION : return "relative fluid acceleration"; break;
        case ANGULAR_FLUID_VELOCITY      : return "angular fluid velocity"     ; break;
        case ANGULAR_FLUID_ACCELERATION  : return "angular fluid acceleration" ; break;
        case FLUID_DILATATION            : return "fluid dilatation"           ; break;
        case FLUID_DILATATION_TDERIV     : return "fluid dilatation tderiv"    ; break;
    }
    assert(false);
    return nullptr;
}

void FEBioPolarFluid::InitModule()
{
    FECoreKernel& febio = FECoreKernel::GetInstance();
    
    // register domain
    febio.RegisterDomain(new FEPolarFluidDomainFactory);
    
    // define the polar fluid module
    febio.CreateModule(new FEPolarFluidModule, "polar fluid",
                       "{"
                       "   \"title\" : \"Polar Fluid\","
                       "   \"info\"  : \"Polar fluid analysis.\""
                       "}");
    
    febio.SetModuleDependency("fluid");
    
    //-----------------------------------------------------------------------------
    // analyis classes (default type must match module name!)
    REGISTER_FECORE_CLASS(FEPolarFluidAnalysis, "polar fluid");
    
    REGISTER_FECORE_CLASS(FEPolarFluidSolver, "polar fluid");
    
    REGISTER_FECORE_CLASS(FEPolarFluid, "polar fluid");
    
    REGISTER_FECORE_CLASS(FEPolarFluidDomain3D, "polar-fluid-3D");
    
    REGISTER_FECORE_CLASS(FETangentialFlowPFStabilization, "fluid tangential stabilization");
    
    febio.SetActiveModule(0);
}