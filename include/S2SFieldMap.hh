// -*- C++ -*-

#ifndef S2S_FIELD_MAP_HH
#define S2S_FIELD_MAP_HH

#include <vector>

#include <G4String.hh>
#include <G4ThreeVector.hh>

//______________________________________________________________________________
class S2SFieldMap
{
public:
  static G4String ClassName();
  S2SFieldMap(const G4String& file_name);
  ~S2SFieldMap();

private:
  S2SFieldMap(const S2SFieldMap&);
  S2SFieldMap& operator =(const S2SFieldMap&);

private:
  typedef std::vector< std::vector< std::vector<G4ThreeVector> > > Field;
  G4bool      m_is_ready;
  G4String    m_file_name;
  Field       m_b;
  G4int       m_nx, m_ny, m_nz;
  G4double    m_xmin, m_ymin, m_zmin;
  G4double    m_xmax, m_ymax, m_zmax;
  G4double    m_dx, m_dy, m_dz;
  G4double    m_value_calc;
  G4double    m_value_nmr;
  G4ThreeVector m_field_size;

public:
  G4bool Initialize();
  G4bool IsInsideField(G4double* pos) const;
  G4bool IsInsideField(const G4ThreeVector& pos) const;
  G4bool IsReady() const { return m_is_ready; }
  const G4ThreeVector& GetFieldSize() const { return m_field_size; }
  G4bool GetFieldValue(const G4double pointCM[3],
			G4double *BfieldTesla) const;
  void   SetValueCalc(G4double v){ m_value_calc = v; }
  void   SetValueNMR(G4double v){ m_value_nmr = v; }

private:
  void ClearField();
};

//______________________________________________________________________________
inline G4String
S2SFieldMap::ClassName()
{
  static G4String s_name("S2SFieldMap");
  return s_name;
}

#endif
