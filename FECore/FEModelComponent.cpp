#include "stdafx.h"
#include "FEModelComponent.h"
#include <string.h>
#include "DumpStream.h"

//-----------------------------------------------------------------------------
//! The constructor takes two arguments: the SUPER_CLASS_ID which defines the 
//! type of the model component and a pointer to the FEModel object this component
//! belongs to.
FEModelComponent::FEModelComponent(SUPER_CLASS_ID sid, FEModel* fem) : FECoreBase(fem, sid)
{
	// assign a class ID
	static int nid = 1;
	m_nClassID = nid++;

	// the ID can be used by derived class to define a identifier for derived classes
	// This value needs to be set in the constructor
	m_nID = 0;

	// initialize parameters
	m_bactive = true;
}

//-----------------------------------------------------------------------------
FEModelComponent::~FEModelComponent()
{
	
}

//-----------------------------------------------------------------------------
int FEModelComponent::GetID() const
{
	return m_nID;
}

//-----------------------------------------------------------------------------
void FEModelComponent::SetID(int n)
{
	m_nID = n;
}

//-----------------------------------------------------------------------------
int FEModelComponent::GetClassID() const
{
	return m_nClassID;
}

//-----------------------------------------------------------------------------
void FEModelComponent::SetClassID(int n)
{
	m_nClassID = n;
}

//-----------------------------------------------------------------------------
bool FEModelComponent::Init()
{
	return true;
}

//-----------------------------------------------------------------------------
bool FEModelComponent::IsActive() const
{ 
	return m_bactive; 
}

//-----------------------------------------------------------------------------
void FEModelComponent::Activate()
{ 
	m_bactive = true; 
}

//-----------------------------------------------------------------------------
void FEModelComponent::Deactivate()
{ 
	m_bactive = false; 
}

//-----------------------------------------------------------------------------
void FEModelComponent::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	FECoreBase::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_nID;
		ar << m_nClassID;
		ar << m_bactive;
	}
	else
	{
		ar >> m_nID;
		ar >> m_nClassID;
		ar >> m_bactive;
	}
}
