/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
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
#pragma once
#include "FEDiscreteElasticMaterial.h"
#include <FECore/FEFunction1D.h>
#include "febiomech_api.h"

//-----------------------------------------------------------------------------
//! This class implements a non-linear spring where users can choose what 
//! deformation measure to use in the force curve. 
class FEBIOMECH_API FEGeneralSpringMaterial : public FEDiscreteElasticMaterial
{
public:
	FEGeneralSpringMaterial(FEModel* pfem);

	vec3d Force(FEDiscreteMaterialPoint& mp) override;
	mat3d Stiffness(FEDiscreteMaterialPoint& mp) override;
	double StrainEnergy(FEDiscreteMaterialPoint& mp) override;

private:
	FEFunction1D*	m_F;		//!< force-displacement function
	double			m_scale;	//!< scale factor for force
	int				m_measure;	//!< deformation measure

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
