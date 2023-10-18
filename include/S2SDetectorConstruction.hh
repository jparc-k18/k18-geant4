// -*- C++ -*-

#ifndef S2S_DETECTOR_CONSTRUCTION_HH
#define S2S_DETECTOR_CONSTRUCTION_HH 1

#include <G4VUserDetectorConstruction.hh>
#include <globals.hh>

class G4Box;
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
  void ConstructField();
  //     void ConstructCalorimeter();
  void ConstructTarget ();
  void ConstructQ1();
  void ConstructQ2();
  void ConstructD1();
  void ConstructSDC1();
  void ConstructSDC2();
  void ConstructSDC3();
  void ConstructSDC4();
  void ConstructSDC5();
  void ConstructTOF();
  void ConstructAC1();
  void ConstructWC();
  void ConstructVP();

  void SetTarget();
  void MakePositionDetector();
  void MakeTOFCounter();
  void MakeAerogelCounter();
  void MakeWaterCounter();
  G4bool IsVolumeStopper( G4VPhysicalVolume *physVol ) const ;

};

#endif
