#pragma once
#include "FEBioMech/FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Neo Hookean material

//! Implementation of a neo-Hookean hyperelastic material.
class FENeoHookeanPI : public FEElasticMaterial
{
public:
	FENeoHookeanPI(FEModel* pfem) : FEElasticMaterial(pfem) {}

public:
	double	m_E;	//!< Young's modulus
	double	m_v;	//!< Poisson's ratio

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! data initialization and checking
	void Init();

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
