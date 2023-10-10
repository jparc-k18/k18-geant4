/*
  DCGeomMan.hh
*/

#ifndef DCGeomMan_h
#define DCGeomMan_h 1

#include <string>
#include <vector>
#include <map>

#include "G4ThreeVector.hh"

class DCGeomRecord;

class DCGeomMan
{
public:
  static G4String   ClassName();
  static DCGeomMan& GetInstance();
  ~DCGeomMan();

private:
  DCGeomMan();
  DCGeomMan(const DCGeomMan&);
  DCGeomMan& operator=(const DCGeomMan&);

private:
  using DCGeomContainer = std::map <G4int, DCGeomRecord*>;
  using DCGeomIterator = DCGeomContainer::const_iterator;
  using IntList = std::map<G4String, G4int>;
  using DoubleList = std::map<G4String, G4double>;
  G4bool          m_is_ready;
  G4String        m_file_name;
  DCGeomContainer m_container;
  IntList         m_detector_id_map;
  DoubleList      m_global_z_map;
  DoubleList      m_local_z_map;

public:
  G4int               CalcWireNumber(G4int lnum, G4double position) const;
  G4int               CalcWireNumber(const G4String& key, G4double position) const;
  G4double            CalcWirePosition(G4int lnum, G4double wire) const;
  G4double            CalcWirePosition(const G4String& key, G4double wire) const;
  void                Clear();
  std::vector<G4int>  GetDetectorIDList() const;
  G4int               GetDetectorId(const G4String &key) const;
  const G4ThreeVector& GetGlobalPosition(G4int lnum) const;
  const G4ThreeVector& GetGlobalPosition(const G4String& key) const;
  G4double            GetLocalZ(G4int lnum) const;
  G4double            GetLocalZ(const G4String& key) const;
  const DCGeomRecord* GetRecord(G4int lnum) const;
  const DCGeomRecord* GetRecord(const G4String& key) const;
  G4double            GetResolution(G4int lnum) const;
  G4double            GetResolution(const G4String& key) const;
  G4double            GetRotAngle1(G4int lnum) const;
  G4double            GetRotAngle1(const G4String& key) const;
  G4double            GetRotAngle2(G4int lnum) const;
  G4double            GetRotAngle2(const G4String& key) const;
  G4double            GetTiltAngle(G4int lnum) const;
  G4double            GetTiltAngle(const G4String& key) const;
  G4double            GetWirePitch(G4int lnum) const;
  G4double            GetWirePitch(const G4String& key) const;
  G4ThreeVector       Global2LocalDir(G4int lnum, const G4ThreeVector &in) const;
  G4ThreeVector       Global2LocalDir(const G4String& key, const G4ThreeVector &in) const;
  G4ThreeVector       Global2LocalPos(G4int lnum, const G4ThreeVector &in) const;
  G4ThreeVector       Global2LocalPos(const G4String& key, const G4ThreeVector &in) const;
  G4bool              Initialize();
  G4bool              Initialize(const G4String& file_name);
  G4bool              IsReady() const { return m_is_ready; }
  G4ThreeVector       Local2GlobalDir(G4int lnum, const G4ThreeVector &in) const;
  G4ThreeVector       Local2GlobalDir(const G4String& key, const G4ThreeVector &in) const;
  G4ThreeVector       Local2GlobalPos(G4int lnum, const G4ThreeVector &in) const;
  G4ThreeVector       Local2GlobalPos(const G4String& key, const G4ThreeVector &in) const;
  G4ThreeVector       NormalVector(G4int lnum) const;
  G4ThreeVector       NormalVector(const G4String& key) const;
  G4ThreeVector       UnitVector(G4int lnum) const;
  G4ThreeVector       UnitVector(const G4String& key) const;
  void                SetFileName(const G4String &file_name);
  // Do not use this method except for special cases
  void                SetResolution(G4int lnum, G4double res);
  void                SetResolution(const G4String& key, G4double res);

  // Static method
  static const G4int&    DetectorId(const G4String& key);
  static const G4double& GlobalZ(const G4String& key);
  static const G4double& LocalZ(const G4String& key);
};

//_____________________________________________________________________________
inline G4String
DCGeomMan::ClassName()
{
  static const G4String s_name("DCGeomMan");
  return s_name;
}

//_____________________________________________________________________________
inline DCGeomMan&
DCGeomMan::GetInstance()
{
  static DCGeomMan g_instance;
  return g_instance;
}

//_____________________________________________________________________________
inline const G4int&
DCGeomMan::DetectorId(const G4String& key)
{
  return GetInstance().m_detector_id_map[key];
}

//_____________________________________________________________________________
inline const G4double&
DCGeomMan::GlobalZ(const G4String& key)
{
  return GetInstance().m_global_z_map[key];
}

//_____________________________________________________________________________
inline const G4double&
DCGeomMan::LocalZ(const G4String& key)
{
  return GetInstance().m_local_z_map[key];
}

#endif
