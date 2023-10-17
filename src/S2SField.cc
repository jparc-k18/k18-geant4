// -*- C++ -*-

#include "S2SField.hh"

#include <fstream>
#include <iomanip>

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4ThreeVector.hh>
#include <G4TwoVector.hh>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetSizeMan.hh"
#include "S2SFieldMap.hh"
#include "FuncName.hh"
#include "PrintHelper.hh"

namespace
{
const auto& gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gSize = DetSizeMan::GetInstance();
}

//_____________________________________________________________________________
S2SField::S2SField(const G4String& file_name)
  : m_is_ready(false),
    m_field_map(new S2SFieldMap(file_name))
{
  Initialize();
}

//_____________________________________________________________________________
S2SField::~S2SField()
{
}

//_____________________________________________________________________________
G4bool
S2SField::Initialize()
{
  G4cout << FUNC_NAME << G4endl;
  m_field_map->Initialize();
  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
void
S2SField::GetFieldValue(const G4double Point[4], G4double* Bfield) const
{
  Bfield[0] = 0.*CLHEP::tesla;
  Bfield[1] = 0.*CLHEP::tesla;
  Bfield[2] = 0.*CLHEP::tesla;

  if (!m_is_ready)
    return;

  G4double pos[3] =
    { Point[0]/CLHEP::cm, Point[1]/CLHEP::cm, Point[2]/CLHEP::cm };
  m_field_map->GetFieldValue(pos, Bfield);

#if 0
  G4ThreeVector b(Bfield[0], Bfield[1], Bfield[2]);
  if(b.mag() > 0.001*CLHEP::tesla
     || true
     ){
    PrintHelper helper(4, std::ios::fixed, G4cout);
    G4cout << FUNC_NAME << " X=("
	   << std::setw(10) << Point[0] << " "
	   << std::setw(10) << Point[1] << " "
	   << std::setw(10) << Point[2] << " "
	   << std::setw(10) << Point[3] << "), B=("
	   << std::setw(10) << Bfield[0]/CLHEP::tesla << " "
	   << std::setw(10) << Bfield[1]/CLHEP::tesla << " "
	   << std::setw(10) << Bfield[2]/CLHEP::tesla << ")" << G4endl;
  }
#endif

  return;
}

//_____________________________________________________________________________
const G4ThreeVector&
S2SField::GetSizeField() const
{
  if(m_field_map){
    return m_field_map->GetFieldSize();
  }else{
    static G4ThreeVector nullvector;
    return nullvector;
  }
}
