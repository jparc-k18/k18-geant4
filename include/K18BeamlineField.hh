// -*- C++ -*-

#ifndef K18_BEAMLINE_FIELD_HH
#define K18_BEAMLINE_FIELD_HH 1

#include <memory>

#include <G4MagneticField.hh>
#include <G4String.hh>
#include <G4ThreeVector.hh>

class S2SField;
class K18FieldMap;

// Converts points and vectors between the native QQDQQ frame and S-2S.
class K18BeamlineFrame
{
public:
  K18BeamlineFrame(const G4ThreeVector& target_global, G4double target_l);

  G4ThreeVector InternalToS2SPoint(const G4ThreeVector& point) const;
  G4ThreeVector S2SToInternalPoint(const G4ThreeVector& point) const;
  G4ThreeVector InternalToS2SVector(const G4ThreeVector& vec) const;
  G4ThreeVector S2SToInternalVector(const G4ThreeVector& vec) const;

  static G4double BendAngle();
  static G4ThreeVector Axis(G4double angle_deg);
  static G4ThreeVector XAxis(G4double angle_deg);
  static G4ThreeVector YAxis();
  static G4double VinL();
  static G4ThreeVector Vin();
  static G4double VoutL();
  static G4ThreeVector Vout();
  static G4double Q10UpstreamOffset();

private:
  G4ThreeVector m_target_global;
  G4ThreeVector m_target_internal;
};

// Provides the QQDQQ magnetic field from either an analytic model or a map.
class K18BeamlineField : public G4MagneticField
{
public:
  K18BeamlineField();
  explicit K18BeamlineField(const G4String& scale_prefix);
  ~K18BeamlineField() override;

  void GetFieldValue(const G4double point[4], G4double* bfield) const override;
  G4bool Contains(const G4ThreeVector& pos, G4int i) const;
  G4bool ContainsNominal(const G4ThreeVector& pos, G4int i) const;
  G4bool UsesFieldMap() const { return static_cast<G4bool>(m_field_map); }

private:
  struct Quad {
    G4double b0;
    G4double a0;
    G4double l;
    G4double z;
  };
  struct Sect {
    G4double b0;
    G4double rho;
    G4double width;
    G4double half_gap;
    G4ThreeVector center;
    G4double bend_angle;
    G4double alpha;
  };
  enum MagnetType { kQuad, kSect };
  struct Magnet {
    MagnetType type;
    G4String name;
    G4double scale;
    G4double theta;
    Quad quad;
    Sect sect;
  };

  std::shared_ptr<const K18FieldMap> m_field_map;
  Magnet m_magnet[5];

  G4bool InMagnet(const G4ThreeVector& pos, G4int i) const;
  void CalcQuad(G4int i, const G4ThreeVector& pos, G4double* bfield) const;
  void CalcSect(G4int i, const G4ThreeVector& pos, G4double* bfield) const;
};

// Adapts the native QQDQQ field to Geant4 global coordinates.
class K18BeamlineFieldAdapter : public G4MagneticField
{
public:
  K18BeamlineFieldAdapter(const G4ThreeVector& target_global, G4double target_l);
  ~K18BeamlineFieldAdapter() override;

  void GetFieldValue(const G4double point[4], G4double* bfield) const override;
  G4bool UsesFieldMap() const { return m_field.UsesFieldMap(); }

private:
  K18BeamlineFrame m_frame;
  K18BeamlineField m_field;
};

// Adds the upstream QQDQQ and downstream S-2S fields in the common world.
class K18CombinedField : public G4MagneticField
{
public:
  K18CombinedField(const G4String& s2s_field_map,
                   const G4ThreeVector& target_global,
                   G4double target_l);
  ~K18CombinedField() override;

  void GetFieldValue(const G4double point[4], G4double* bfield) const override;
  G4bool UsesK18FieldMap() const { return m_k18_field.UsesFieldMap(); }

private:
  std::unique_ptr<S2SField> m_s2s_field;
  K18BeamlineFieldAdapter m_k18_field;
};

#endif
