// FEMaterial.h: interface for the FEMaterial class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_)
#define AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "tens4d.h"
#include "FECoreBase.h"
#include "FEMaterialPoint.h"
#include "FECoordSysMap.h"
#include "DumpStream.h"
#include "FECoreKernel.h"
#include "FEModelParam.h"
#include "FEDomainList.h"
#include <string.h>
#include <stddef.h>

#define INRANGE(x, a, b) ((x)>=(a) && (x)<=(b))
#define IN_RIGHT_OPEN_RANGE(x, a, b) ((x)>=(a) && (x)<(b))

//-----------------------------------------------------------------------------
// forward declaration of some classes
class FEModel;
class FEElement;
class FEDomain;

//-----------------------------------------------------------------------------
// Forward declaration of the FEElasticMaterial class. 
// TODO: The only reason I had to do this is to define the FEMaterial::GetElasticMaterial.
// However, this is only a temporary construct so make sure to delete this forward declaration
// when no longer needed.
class FEElasticMaterial;

//-----------------------------------------------------------------------------
//! Abstract base class for material types

//! From this class all other material classes are derived.

class FECORE_API FEMaterial : public FECoreBase
{
	DECLARE_SUPER_CLASS(FEMATERIAL_ID);

public:
	FEMaterial(FEModel* fem);
	virtual ~FEMaterial();

	//! returns a pointer to a new material point object
	virtual FEMaterialPoint* CreateMaterialPointData() { return 0; };

	//! performs initialization
	bool Init();

	//! Serialize material data to archive
	void Serialize(DumpStream& ar);

	//! Return elastic material \todo I need to move this function up the hierarchy once I redesign the material library
	virtual FEElasticMaterial* GetElasticMaterial() { return 0; }
    
    //! Update specialized material points at each iteration
    virtual void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp) {}

public:
	//! Set the local coordinate system map
	void SetCoordinateSystemMap(FECoordSysMap* pmap);

	//! Get the local coordinate system
	FECoordSysMap* GetCoordinateSystemMap();

	//! Set the local coordinate for a material point
	virtual void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp);

public:
	template <class T> T* ExtractProperty();

public:
	//! Assign a domain to this material
	void AddDomain(FEDomain* dom);

	//! get the domaint list
	FEDomainList& GetDomainList() { return m_domList; }

private:
	FECoordSysMap*	m_pmap;			//!< local material coordinate system
	FEDomainList	m_domList;		//!< list of domains that use this material

	DECLARE_FECORE_CLASS();
};

template <class T> T* FEMaterial::ExtractProperty()
{
	if (dynamic_cast<T*>(this)) return dynamic_cast<T*>(this);

	int NC = Properties();
	for (int i = 0; i < NC; i++)
	{
		FEMaterial* pmi = static_cast<FEMaterial*>(GetProperty(i));
		T* pm = pmi->ExtractProperty<T>();
		if (pm) return pm;
	}

	return nullptr;
}

#endif // !defined(AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_)
