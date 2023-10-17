// -*- C++ -*-

#ifndef MAGNETIC_FIELD_HH
#define MAGNETIC_FIELD_HH

#include <vector>
#include <map>
#include <globals.hh>
#include <G4MagneticField.hh>
#include <G4String.hh>
#include <G4ThreeVector.hh>

class S2SFieldMap;

//_____________________________________________________________________________
class MagnetInfo
{
public:
  static G4String ClassName();
  MagnetInfo(G4String n="")
    : name(n)
  {}
  enum EMagnetType { kDipole, kQuadrupole, kShs, kKurama, NMagnetType };

public:
  G4int    type;
  G4String name;
  G4double b0;        // magnetic field
  G4double l;         // length
  G4double ra1;       // rotation angle around z-axis
  G4ThreeVector pos;  // center position
  G4ThreeVector size; // half width, height, and thickness
  // for Dipole
  G4double rho;   // bending radius
  G4double bend;  // bending angle
  G4double alpha; // angle between conetral track and entrance pole
  G4double beta;  // angle between conetral track and exit pole
  // for Quadrupole
  G4double a0; // aperture

public:
  G4bool CalcK18Field(const G4ThreeVector& point, G4double* bfield) const;
};

//_____________________________________________________________________________
inline G4String
MagnetInfo::ClassName()
{
  static G4String s_name("MagnetInfo");
  return s_name;
}

//_____________________________________________________________________________
class S2SField : public G4MagneticField
{
public:
  static G4String ClassName();
  S2SField(const G4String& file_name);
  ~S2SField();

private:
  S2SField(const S2SField&);
  S2SField& operator =(const S2SField&);

private:
  // using MagnetMap = std::map<G4String, MagnetInfo>;
  using Field = std::vector<std::vector<std::vector<G4ThreeVector>>>;
  G4bool    m_is_ready;
  S2SFieldMap* m_field_map;

public:
  virtual void GetFieldValue(const G4double Point[4], G4double* Bfield) const;

public:
  const G4ThreeVector& GetSizeField() const;
  G4bool Initialize();
  G4bool IsReady() const { return m_is_ready; }
};

//_____________________________________________________________________________
inline G4String
S2SField::ClassName()
{
  static G4String s_name("S2SField");
  return s_name;
}

#endif
