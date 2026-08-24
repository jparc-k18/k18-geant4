// -*- C++ -*-

#ifndef K18_FIELD_MAP_HH
#define K18_FIELD_MAP_HH 1

#include <cstddef>
#include <vector>

#include <globals.hh>
#include <G4String.hh>
#include <G4ThreeVector.hh>

// Regular K1.8 native-coordinate field map.  File coordinates are cm and
// field components are Tesla; query positions use Geant4 internal length.
class K18FieldMap
{
public:
  static G4String ClassName();
  K18FieldMap(const G4String& file_name,
              G4double value_measure=1., G4double value_calc=1.);
  ~K18FieldMap() = default;

  G4bool Initialize();
  G4bool GetFieldValue(const G4double point[4], G4double* bfield) const;

  G4bool IsReady() const { return m_is_ready; }
  G4double Scale() const { return m_scale; }
  const G4String& FileName() const { return m_file_name; }

private:
  std::size_t Index(G4int ix, G4int iy, G4int iz) const;
  void Clear();

  G4bool m_is_ready;
  G4String m_file_name;
  std::vector<G4ThreeVector> m_field;
  G4int m_nx;
  G4int m_ny;
  G4int m_nz;
  G4double m_xmin;
  G4double m_ymin;
  G4double m_zmin;
  G4double m_dx;
  G4double m_dy;
  G4double m_dz;
  G4double m_scale;
};

inline G4String
K18FieldMap::ClassName()
{
  static const G4String name("K18FieldMap");
  return name;
}

#endif
