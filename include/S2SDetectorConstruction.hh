// -*- C++ -*-

#ifndef S2SDetectorConstruction_h
#define S2SDetectorConstruction_h 1

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
  G4Box*             solidWorld;    //pointer to the solid World
  G4LogicalVolume*   logicWorld;    //pointer to the logical World
  G4VPhysicalVolume* physiWorld;    //pointer to the physical World
  G4UniformMagField* magField;      //pointer to the magnetic field
  MaterialList*      m_material_list;
  G4Material* DMaterial;
  G4bool m_check_overlaps;

  void MakeField();
  //     void ConstructCalorimeter(G4VPhysicalVolume* pMother);
  void ConstructTarget (G4VPhysicalVolume* pMother);
  void ConstructQ1(G4VPhysicalVolume* pMother);
  void ConstructQ2(G4VPhysicalVolume* pMother);
  void ConstructD1(G4VPhysicalVolume* pMother);

  void SetTarget(G4VPhysicalVolume* pMother);
  void MakePositionDetector(G4VPhysicalVolume* pMother);
  void MakeTOFCounter(G4VPhysicalVolume* pMother);
  void MakeAerogelCounter(G4VPhysicalVolume* pMother);
  void MakeWaterCounter(G4VPhysicalVolume* pMother);
  void MakeSlits(G4VPhysicalVolume* pMother);
  G4bool IsVolumeStopper( G4VPhysicalVolume *physVol ) const ;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
