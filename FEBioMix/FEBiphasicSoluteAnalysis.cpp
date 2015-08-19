#include "stdafx.h"
#include "FEBiphasicSoluteAnalysis.h"
#include "FEBioMech/FERigidMaterial.h"
#include "FEBioMech/FEContactInterface.h"
#include "FEBioMech/FEUncoupledMaterial.h"
#include "FEBiphasic.h"
#include "FEBiphasicSolute.h"
#include "FETriphasic.h"
#include "FECore/FERigidBody.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
void FEBiphasicSoluteAnalysis::InitNodes()
{
    // get nodal DOFS
    DOFS& fedofs = *DOFS::GetInstance();
    int MAX_CDOFS = fedofs.GetCDOFS();
    
	// open all dofs we need
	FEMesh& mesh = m_fem.GetMesh();
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (int i=0; i<(int)node.m_ID.size(); ++i) node.m_ID[i] = DOF_FIXED;

		if (node.m_bexclude == false)
		{
			if (node.m_rid < 0)
			{
				node.m_ID[DOF_X] = DOF_OPEN;
				node.m_ID[DOF_Y] = DOF_OPEN;
				node.m_ID[DOF_Z] = DOF_OPEN;
			}

			if (node.m_bshell)
			{
				node.m_ID[DOF_U] = DOF_OPEN;
				node.m_ID[DOF_V] = DOF_OPEN;
				node.m_ID[DOF_W] = DOF_OPEN;
			}

			node.m_ID[DOF_P] = DOF_OPEN;

			for (int k=0; k<MAX_CDOFS; ++k) {
				int dofc = DOF_C + k;
				node.m_ID[dofc] = DOF_OPEN;
			}
		}
	}

	// apply fixed dofs
	for (int i=0; i<m_fem.FixedBCs(); ++i)
	{
		FEFixedBC& bc = *m_fem.FixedBC(i);
		bc.Activate();
	}

	// apply prescribed dofs
	int ndis = m_fem.PrescribedBCs();
	for (int i=0; i<ndis; ++i)
	{
		FEPrescribedBC& DC = *m_fem.PrescribedBC(i);
		if (DC.IsActive())
		{
			int dof = DC.GetDOF();
			for (size_t j = 0; j<DC.Items(); ++j)
			{
				FENode& node = m_fem.GetMesh().Node(DC.NodeID(j));
				node.m_ID[dof] = DOF_PRESCRIBED;
			}
		}
	}

	// fix all mixture dofs that are not used that is, that are not part of a biphasic material.
	const int NN = mesh.Nodes();
	vector<int> tag;

	// do the pressure dofs first
	tag.assign(NN, 0);
	for (int nd = 0; nd<mesh.Domains(); ++nd)
	{
		FEDomain& dom = mesh.Domain(nd);
		FEBiphasic*		  bm  = dynamic_cast<FEBiphasic*      >(dom.GetMaterial());
		FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(dom.GetMaterial());
		FETriphasic*      btm = dynamic_cast<FETriphasic*     >(dom.GetMaterial());
		if (bm || bsm || btm)
		{
			for (int i=0; i<dom.Elements(); ++i)
			{
				FEElement& el = dom.ElementRef(i);
				int N = el.Nodes();
				int* n = &el.m_node[0];
				for (int j=0; j<N; ++j) 
					if (mesh.Node(n[j]).m_ID[DOF_P] != DOF_FIXED) tag[ n[j] ] = 1;
			}
		}
	}

	// fix pressure dofs of all unmarked nodes
	for (int nd = 0; nd<mesh.Domains(); ++nd)
	{
		FEDomain& dom = mesh.Domain(nd);
		for (int i=0; i<dom.Elements(); ++i)
		{
			FEElement& el = dom.ElementRef(i);
			int N = el.Nodes();
			int* n = &el.m_node[0];
			for (int j=0; j<N; ++j) {
				if (tag[ n[j] ] != 1) mesh.Node(n[j]).m_ID[DOF_P] = DOF_FIXED;
			}
		}
	}

	// next, do the concentration dofs
	for (int k=0; k<MAX_CDOFS; ++k)
	{
		int dofc = DOF_C + k;
		tag.assign(NN, 0);
		for (int nd = 0; nd<mesh.Domains(); ++nd)
		{
			FEDomain& dom = mesh.Domain(nd);

			// get the material
			FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(dom.GetMaterial());
			FETriphasic*      btm = dynamic_cast<FETriphasic*     >(dom.GetMaterial());

			// see if this material has this concentration dof
			bool has_c = false;
			if (bsm && (dofc == DOF_C + bsm->GetSolute()->GetSoluteID())) has_c = true;
			if (btm && (dofc == DOF_C + btm->m_pSolute[0]->GetSoluteID())) has_c = true;
			if (btm && (dofc == DOF_C + btm->m_pSolute[1]->GetSoluteID())) has_c = true;

			// if so, mark all non-fixed dofs
			if (has_c)
			{
				for (int i=0; i<dom.Elements(); ++i)
				{
					FEElement& el = dom.ElementRef(i);
					int N = el.Nodes();
					int* n = &el.m_node[0];
					for (int j=0; j<N; ++j) {
						if (mesh.Node(n[j]).m_ID[dofc] != DOF_FIXED) tag[ n[j] ] = 1;
					}
				}
			}
		}
	
		// step 2. fix concentration dofs of all unmarked nodes
		for (int nd = 0; nd<mesh.Domains(); ++nd)
		{
			FEDomain& dom = mesh.Domain(nd);
			for (int i=0; i<dom.Elements(); ++i)
			{
				FEElement& el = dom.ElementRef(i);
				int N = el.Nodes();
				int* n = &el.m_node[0];
				for (int j=0; j<N; ++j) {
					if (tag[ n[j]] != 1) mesh.Node(n[j]).m_ID[dofc] = DOF_FIXED;
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! This function is called before the analysis is solved and initializes all
//! analysis data, such as determine active boundary conditions, initializes
//! equation numbers (the latter is actually done by the FESolver class).
bool FEBiphasicSoluteAnalysis::Activate()
{
	// initialize base class data
	FEAnalysis::Activate();

	// reset nodal ID's
	InitNodes();

	// initialize equations
	// ----->
	if (m_psolver->InitEquations() == false) return false;

	// initialize linear constraints
	// Must be done after equations are initialized
	if (InitLinearConstraints() == false) return false;
	// ----->

	// Now we adjust the equation numbers of prescribed dofs according to the above rule
	// Make sure that a prescribed dof has not been fixed
	// TODO: maybe this can be moved to the FESolver::InitEquations function
	int ndis = m_fem.PrescribedBCs();
	for (int i=0; i<ndis; ++i)
	{
		FEPrescribedBC& DC = *m_fem.PrescribedBC(i);
		if (DC.IsActive()) DC.Update();
	}

	// modify the linear constraints
	if (m_fem.m_LinC.size())
	{
		list<FELinearConstraint>::iterator il = m_fem.m_LinC.begin();
		for (int l=0; l<(int) m_fem.m_LinC.size(); ++l, ++il) il->Activate();
	}

	// modify the (aug lag) nonlinear constraints
	// TODO: I think this is already done in FEM::Init. Why do I need to do this again?
	int M = m_fem.NonlinearConstraints();
	for (int m=0; m<M; ++m) 
	{
		FENLConstraint* plc = m_fem.NonlinearConstraint(m);
		if (plc->IsActive()) plc->Init();
	}

	// do one time initialization of solver data
	if (m_psolver->Init() == false)
	{
		felog.printbox("FATAL ERROR","Failed to initialize solver.\nAborting run.\n");
		return false;
	}

	return true;
}


