// -*- C++ -*-

#ifndef S2S_DETECTOR_CONSTRUCTION_HH
#define S2S_DETECTOR_CONSTRUCTION_HH 1

#include <map>

#include <G4VUserDetectorConstruction.hh>
#include <G4String.hh>
#include <globals.hh>

class G4Box;
class G4Element;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;

struct MaterialList;

//_____________________________________________________________________________
class S2SDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  S2SDetectorConstruction();
  ~S2SDetectorConstruction();

public:
  G4VPhysicalVolume* Construct();

private:
  G4int                           m_experiment;
  G4bool                          m_check_overlaps;
  G4LogicalVolume*                m_world_lv;
  static std::vector<G4String>    s_detector_list;

  void AddNewDetector(G4VSensitiveDetector* sd);
  void ConstructField();
  //     void ConstructCalorimeter();
  void ConstructBAC1();
  void ConstructSAC();
  void ConstructTarget ();
  void ConstructTgtHeBag();
  void ConstructQ1();
  void ConstructQ2();
  void ConstructD1();
  void ConstructSDC1();
  void ConstructSDC2();
  void ConstructSDC3();
  void ConstructSDC4();
  void ConstructSDC5();
  void ConstructSFT();
  void ConstructSSD();
  void ConstructTOF();
  void ConstructAC1();
  void ConstructWC();
  void ConstructVP();
  void ConstructKLChamber(G4int i);
  void ConstructKLChamberMylar(G4int i);

  void ConstructHTOF();
  void ConstructHypTPC();

  void ConstructRC(G4int i);
  void ConstructPD(G4int i, G4String YZ);
};

#endif
