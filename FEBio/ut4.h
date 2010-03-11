#pragma once
#include "FEDomain.h"
#include "FENodeElemList.h"

//-----------------------------------------------------------------------------
//! Domain for nodally integrated tet elements 

//! This class implements the uniform nodal strain tetrahedron with
//! isochoric stabilization as described by Gee, Dohrmann, Key and Wall (2009)
//!
class FEUT4Domain : public FEElasticSolidDomain
{
	struct UT4NODE
	{
		int		inode;	// index of FE node
		double	Vi;		// reference nodal volume
		double	vi;		// current nodal volume
		mat3d	Fi;		// average deformation gradient
		mat3ds	si;		// nodal stress
	};

public:
	//! constructor
	FEUT4Domain(FEMesh* pm);

	//! clone function
	FEDomain* Clone()
	{
		FEUT4Domain* pd = new FEUT4Domain(m_pMesh);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh;
		pd->m_tag = m_tag;
		pd->m_alpha = m_alpha;
		pd->m_Node = m_Node;
		return pd;
	}

	//! initialization function
	bool Initialize(FEM& fem);

	//! Update stresses
	void UpdateStresses(FEM& fem);

	//! calculates the total residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolidSolver* psolver);

protected:
	//! calculates the nodal residual
	void NodalResidual(FESolidSolver* psolver, vector<double>& R);

	//! calculates the element residual
	void ElementResidual(FESolidSolver* psolver, vector<double>& R);

	//! Calculates the internal stress vector for solid elements
	void ElementInternalForces(FESolidElement& el, vector<double>& fe);

protected:
	//! Calculates the element stiffness matrix
	void ElementalStiffnessMatrix(FESolidSolver* psolver);

	//! calculates the solid element stiffness matrix
	void ElementStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! Calculates the nodal stiffness matrix
	void NodalStiffnessMatrix(FESolidSolver* psolver);

	//! geometrical stiffness (i.e. initial stress)
	void GeometricalStiffness(FESolidElement& el, matrix& ke);

	//! material stiffness component
	void MaterialStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! nodal geometry stiffness contribution
	void NodalGeometryStiffness(UT4NODE& node);

	//! nodal material stiffness contribution
	void NodalMaterialStiffness(UT4NODE& node);

protected:
	//! calculate the volume of a tet element
	double TetVolume(vec3d* r);

private:
	double	m_alpha;	//!< stabilization factor alpha

	vector<int>		m_tag;	//!< nodal tags
	vector<UT4NODE>	m_Node;	//!< Nodal data

	FENodeElemList	m_NEL;
};
