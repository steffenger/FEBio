/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
#include "FEPrescribedDisplacement.h"
#include <FECore/FEModel.h>
#include <FECore/DOFS.h>

//=======================================================================================
// NOTE: I'm setting FEBoundaryCondition is the base class since I don't want to pull
//       in the parameters of FEPrescribedDOF. 
BEGIN_FECORE_CLASS(FEPrescribedDisplacement, FEBoundaryCondition)
	ADD_PARAMETER(m_comp, "dof", 0, "x-displacement\0y-displacement\0z-displacement\0");
	ADD_PARAMETER(m_scale, "value");
	ADD_PARAMETER(m_brelative, "relative");

	ADD_PROPERTY(m_nodeSet, "node_set", FEProperty::Reference);
END_FECORE_CLASS();

FEPrescribedDisplacement::FEPrescribedDisplacement(FEModel* fem) : FEPrescribedDOF(fem)
{
	m_comp = 0;
}

bool FEPrescribedDisplacement::Init()
{
	FEModel* fem = GetFEModel();
	int ndof = -1;
	switch (m_comp)
	{
	case 0: ndof = fem->GetDOFIndex("x"); break;
	case 1: ndof = fem->GetDOFIndex("y"); break;
	case 2: ndof = fem->GetDOFIndex("z"); break;
	}
	assert(ndof >= 0);
	SetDOF(ndof);
	return FEPrescribedDOF::Init();
}
