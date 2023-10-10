// -*- C++ -*-

#ifndef DC_GEOM_RECORD_HH
#define DC_GEOM_RECORD_HH

#include <G4ThreeVector.hh>

//_____________________________________________________________________________
class DCGeomRecord
{
public:
  static G4String ClassName();
  DCGeomRecord(G4int id, const G4String& name,
               G4double x, G4double y, G4double z, G4double ta,
               G4double ra1, G4double ra2, G4double length, G4double resol,
               G4double w0, G4double dd, G4double ofs);
  DCGeomRecord(G4int id, const G4String& name,
               const G4ThreeVector pos, G4double ta,
               G4double ra1, G4double ra2, G4double length, G4double resol,
               G4double w0, G4double dd, G4double ofs);
  ~DCGeomRecord();

private:
  DCGeomRecord(const DCGeomRecord&);
  DCGeomRecord& operator =(const DCGeomRecord);

private:
  G4int       m_id;
  G4String     m_name;
  G4ThreeVector m_pos;
  G4double    m_tilt_angle;
  G4double    m_rot_angle1;
  G4double    m_rot_angle2;
  G4double    m_length;
  G4double    m_resolution;
  G4double    m_w0;
  G4double    m_dd;
  G4double    m_offset;

  G4double m_dxds, m_dxdt, m_dxdu;
  G4double m_dyds, m_dydt, m_dydu;
  G4double m_dzds, m_dzdt, m_dzdu;

  G4double m_dsdx, m_dsdy, m_dsdz;
  G4double m_dtdx, m_dtdy, m_dtdz;
  G4double m_dudx, m_dudy, m_dudz;

public:
  const G4ThreeVector& Position()       const { return m_pos; }
  G4ThreeVector        NormalVector()   const;
  G4ThreeVector        UnitVector()     const;
  G4int              Id()             const { return m_id;         }
  G4String            Name()           const { return m_name;       }
  const G4ThreeVector& Pos()            const { return m_pos;        }
  G4double           TiltAngle()      const { return m_tilt_angle; }
  G4double           RotationAngle1() const { return m_rot_angle1; }
  G4double           RotationAngle2() const { return m_rot_angle2; }
  G4double           Length()         const { return m_length;     }
  G4double           Resolution()     const { return m_resolution; }
  void               SetResolution(G4double res) { m_resolution = res; }
  G4double           WirePitch()      const { return m_dd;         }

  G4double dsdx() const { return m_dsdx; }
  G4double dsdy() const { return m_dsdy; }
  G4double dsdz() const { return m_dsdz; }
  G4double dtdx() const { return m_dtdx; }
  G4double dtdy() const { return m_dtdy; }
  G4double dtdz() const { return m_dtdz; }
  G4double dudx() const { return m_dudx; }
  G4double dudy() const { return m_dudy; }
  G4double dudz() const { return m_dudz; }

  G4double dxds() const { return m_dxds; }
  G4double dxdt() const { return m_dxdt; }
  G4double dxdu() const { return m_dxdu; }
  G4double dyds() const { return m_dyds; }
  G4double dydt() const { return m_dydt; }
  G4double dydu() const { return m_dydu; }
  G4double dzds() const { return m_dzds; }
  G4double dzdt() const { return m_dzdt; }
  G4double dzdu() const { return m_dzdu; }

  G4double WirePos(G4double wire)   const;
  G4int    WireNumber(G4double pos) const;
  void     Print() const;

private:
  void CalcVectors();
};

//_____________________________________________________________________________
inline G4String
DCGeomRecord::ClassName()
{
  static const G4String s_name("DCGeomRecord");
  return s_name;
}

//______________________________________________________________________________
struct DCGeomRecordComp
  : public std::binary_function <DCGeomRecord*, DCGeomRecord*, G4bool>
{
  G4bool operator()(const DCGeomRecord* const p1,
                    const DCGeomRecord* const p2) const
  { return p1->Id() < p2->Id(); }
};

#endif
