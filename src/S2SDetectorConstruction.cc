/*
  "S2SDetectorConstruction.cc"

  This was originally coded by S.Kanatsuki.
  Modified by Toshi Gogami , 11Nov2014.
*/

#include "S2SDetectorConstruction.hh"

#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4TransportationManager.hh"

#include <G4Polyhedra.hh>
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4ExtrudedSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4SimpleRunge.hh"
#include <G4VisAttributes.hh>

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4UserLimits.hh"
#include "G4SystemOfUnits.hh"

#include "G4Colour.hh"
#include "G4SDManager.hh"

#include "MagnetConstant.hh"
#include "S2SField.hh"
#include "MaterialList.hh"
#include "RadDeg.hh"
#include "DCSD.hh"
#include "SlSD.hh"
#include "TOFSD.hh"
#include "ACSD.hh"
#include "WCSD.hh"
//#include "Area.hh"
#include "DetectorID.hh"

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetSizeMan.hh"

//#include "CalorimeterSD.hh"

namespace
{
using namespace CLHEP;
const auto& confMan = ConfMan::GetInstance();
const auto& geomMan = DCGeomMan::GetInstance();
const auto& sizeMan = DetSizeMan::GetInstance();
}

//_____________________________________________________________________________
S2SDetectorConstruction::S2SDetectorConstruction()
  : solidWorld(0),
    logicWorld(0),
    physiWorld(0),
    m_material_list(new MaterialList),
    m_check_overlaps(true)
{
}

//_____________________________________________________________________________
S2SDetectorConstruction::~S2SDetectorConstruction()
{
  delete m_material_list;
}

//_____________________________________________________________________________
G4VPhysicalVolume* S2SDetectorConstruction::Construct()
{
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  ///// World
  const auto& half_size = sizeMan.GetSize("World")*mm/2.;
  solidWorld = new G4Box("World", half_size.x(), half_size.y(), half_size.z());
  logicWorld = new G4LogicalVolume(solidWorld, m_material_list->at("Air"), "World");
  physiWorld = new G4PVPlacement(nullptr, G4ThreeVector(),
				 logicWorld, "World", nullptr, false, 0, m_check_overlaps);
  logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

  MakeField();

#if 1
  ConstructTarget(physiWorld);
#endif

#if 1
  ConstructQ1(physiWorld);
  ConstructQ2(physiWorld);
  ConstructD1(physiWorld);
#endif

#if 1
  ConstructSDC1(physiWorld);
  ConstructSDC2(physiWorld);
#endif

#if 1
  ConstructSDC3(physiWorld);
  ConstructSDC4(physiWorld);
  ConstructSDC5(physiWorld);
  ConstructTOF(physiWorld);
  ConstructAC1(physiWorld);
  ConstructWC(physiWorld);
  // MakeTOFCounter(physiWorld);
  // MakeAerogelCounter(physiWorld);
  // MakeWaterCounter(physiWorld);
#endif

#if 0
  MakeSlits(physiWorld);
#endif

  return physiWorld;
}

//_____________________________________________________________________________
void S2SDetectorConstruction::MakeField()
{
  return;
  S2SField *field = new S2SField(confMan.GetFieldMap());
//    				 confMan.GetMagScaleQ1(),
// 				 confMan.GetMagScaleQ2(),
// 				 confMan.GetMagScale());
  auto fieldManager =
    G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldManager->SetDetectorField(field);
  fieldManager->CreateChordFinder(field);

  // G4Mag_UsualEqRhs* fEquation = new G4Mag_UsualEqRhs(field);
  // G4MagIntegratorStepper *pStepper = new G4SimpleRunge(fEquation);
  // G4ChordFinder *pChordFinder = new G4ChordFinder(field, 1.0E-3*mm, pStepper);
  // fieldManager->SetChordFinder(pChordFinder);

  fieldManager->GetChordFinder()->SetDeltaChord(1.e-3*mm);
  //fieldManager->SetDeltaIntersection( 1.0E-6*mm );
  //fieldManager->SetDeltaOneStep( 1.0E-4*mm );
  //fieldManager->SetMaximumEpsilonStep( 1.0E-3 );
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructTarget(G4VPhysicalVolume *pMother)
{
  const auto& half_size = sizeMan.GetSize("Target")*mm/2.;
  G4Material *TargetMater = nullptr;
  auto Target = confMan.Get<G4String>("TargetMaterial");
  if(Target == "Be"){
    TargetMater = m_material_list->at("Be9");
  }else{
    G4cout << " Sorry, Target: " << Target
           << " is not defined. So Air will be used. " << G4endl;
    TargetMater = m_material_list->at("Air");
  }

  auto TargetBox = new G4Box
    ("TargetBox", half_size.x(), half_size.y(), half_size.z());
  auto logTarget = new G4LogicalVolume(TargetBox, TargetMater, "logTarget");
  G4RotationMatrix rotTarget;

  const auto& pos = geomMan.GetGlobalPosition("Target");
  new G4PVPlacement(G4Transform3D(rotTarget, pos),
                    "physTarget",
                    logTarget,
                    pMother,
                    false,
                    0,
                    m_check_overlaps);
  logTarget->SetVisAttributes(G4Color::Gray());
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructQ1(G4VPhysicalVolume *pMother)
{
  G4double x, y, z;
  G4Material* PoleMater = m_material_list->Fe;
  // G4Material* Q1GapMater = m_material_list->at("Air");
  G4Material* Q1GapMater = m_material_list->at("HeGas");
  const G4double phiStart = 22.5*deg;
  const G4double phiTotal = 360.*deg;
  const G4int numSide   = 8;
  const G4int numZPlane = 2;
  const G4double zPlane[] = { -440.0*mm, 440.0*mm };
  const G4double rInner[] = { 0., 0. };
  const G4double rOuter[] = { 1200.0*mm, 1200.0*mm };
  const G4double a0 = 155*mm; // phi=310mm
  const G4ThreeVector pos(0, 0, -3885.35*mm);
  ///// Pole
  G4VSolid* solidQ1Gap;
  auto solidQ1Box = new G4Box("solidQ1Gap", 540*mm/2, 540*mm/2, 900*mm/2);
  const G4double rCorner = (540*mm*std::sqrt(2)-2*a0)/2.;
  // auto solidCorner = new G4Tubs("solidCorner",
  //                               0*mm, rCorner, 900*mm,
  //                               phiStart, phiTotal);
  G4double zPlaneCorner[] = { -500*mm, 500*mm };
  G4double rOuterCorner[] = { rCorner, rCorner };
  auto solidCorner = new G4Polyhedra("solidCorner", 0*deg, 360*deg, 4, 2,
                                     zPlaneCorner, rInner, rOuterCorner);
  solidQ1Gap = new G4SubtractionSolid("solidQ1Gap", solidQ1Box,
                                      solidCorner, nullptr, G4ThreeVector(270*mm, 270*mm, 0));
  solidQ1Gap = new G4SubtractionSolid("solidQ1Gap", solidQ1Gap,
                                      solidCorner, nullptr, G4ThreeVector(270*mm, -270*mm, 0));
  solidQ1Gap = new G4SubtractionSolid("solidQ1Gap", solidQ1Gap,
                                      solidCorner, nullptr, G4ThreeVector(-270*mm, 270*mm, 0));
  solidQ1Gap = new G4SubtractionSolid("solidQ1Gap", solidQ1Gap,
                                      solidCorner, nullptr, G4ThreeVector(-270*mm, -270*mm, 0));
  G4RotationMatrix rotQ1Box;
  rotQ1Box.rotateZ(45*deg);
  G4double Q1Box3V = 900+1000/sqrt(2);
  G4VSolid* solidQ1Pole = nullptr;
  solidQ1Pole = new G4Box("solidQ1Pole", 2400*mm/2, 2400*mm/2, 880*mm/2);
  auto solidQ1Corner = new G4Box("tmp", 2000*mm/2, 2000*mm/2, 2000*mm/2);
  solidQ1Pole = new G4SubtractionSolid("tmp", solidQ1Pole, solidQ1Corner,
                                       G4Transform3D(rotQ1Box, G4ThreeVector(Q1Box3V,Q1Box3V,0)));
  solidQ1Pole = new G4SubtractionSolid("tmp", solidQ1Pole, solidQ1Corner,
                                       G4Transform3D(rotQ1Box, G4ThreeVector(-Q1Box3V,Q1Box3V,0)));
  solidQ1Pole = new G4SubtractionSolid("tmp", solidQ1Pole, solidQ1Corner,
                                       G4Transform3D(rotQ1Box, G4ThreeVector(Q1Box3V,-Q1Box3V,0)));
  solidQ1Pole = new G4SubtractionSolid("tmp", solidQ1Pole, solidQ1Corner,
                                       G4Transform3D(rotQ1Box, G4ThreeVector(-Q1Box3V,-Q1Box3V,0)));
  solidQ1Pole = new G4SubtractionSolid("solidQ1Pole", solidQ1Pole,
                                       solidQ1Gap, nullptr, G4ThreeVector());
  auto lvQ1Pole = new G4LogicalVolume
    (solidQ1Pole, m_material_list->Fe, "lvQ1Pole");
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), pos),
                    "pvQ1Pole", lvQ1Pole, pMother, false, 0, m_check_overlaps);
  lvQ1Pole->SetVisAttributes(G4Color::Cyan());
  ///// Coil
  G4VSolid* solidQ1Coil = nullptr;
  solidQ1Coil = new G4Box("solidQ1Coil", 1400*mm/2, 1400*mm/2, 180*mm/2);
  solidQ1Coil = new G4SubtractionSolid("solidQ1Coil", solidQ1Coil, solidQ1Box,
                                       nullptr, G4ThreeVector());
  auto lvQ1Coil = new G4LogicalVolume
    (solidQ1Coil, m_material_list->Cu, "lvQ1Coil");
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), pos - G4ThreeVector(0, 0, 530*mm)),
                    "pvQ1CoilU", lvQ1Coil, pMother, false, 0, m_check_overlaps);
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), pos + G4ThreeVector(0, 0, 530*mm)),
                    "pvQ1CoilD", lvQ1Coil, pMother, false, 1, m_check_overlaps);
  lvQ1Coil->SetVisAttributes(G4Color::Brown());
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructQ2(G4VPhysicalVolume *pMother)
{
  G4double x, y, z;
  G4Material *PoleMater = m_material_list->Fe;
  //G4Material *Q2GapMater = m_material_list->at("Air");
  G4Material *Q2GapMater = m_material_list->at("HeGas");
  const G4double a0 = 180*mm; // phi=360mm
  const G4ThreeVector pos(0, 0, -2776.5*mm);
  ///// Pole
  G4VSolid* solidQ2Pole = nullptr;
  solidQ2Pole = new G4Box("solidQ2Pole", 2100*mm/2, 1540*mm/2, 540*mm/2);
  auto solidQ2Gap = new G4Box("solidQ2Gap", 360*mm/2, 360*mm/2, 550*mm/2);
  // 1200x600 mm2;
  G4RotationMatrix rot;
  rot.rotateZ(45*deg);
  solidQ2Pole = new G4SubtractionSolid("solidQ2Pole", solidQ2Pole, solidQ2Gap,
                                       G4Transform3D(rot, G4ThreeVector()));
  auto lvQ2Pole = new G4LogicalVolume(solidQ2Pole, PoleMater, "lvQ2Pole");
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), pos),
                    "pvQ2Pole", lvQ2Pole, pMother, false, 0, m_check_overlaps);
  lvQ2Pole->SetVisAttributes(G4Color::Cyan());
  ///// Coil
  G4VSolid* solidQ2Coil = nullptr;
  solidQ2Coil = new G4Box("solidQ2Coil", 1320*mm/2, 880*mm/2, 130*mm/2);
  auto solidQ2Box = new G4Box("solidQ2Gap", 1200*mm/2, 600*mm/2, 900*mm/2);
  solidQ2Coil = new G4SubtractionSolid("solidQ2Coil", solidQ2Coil, solidQ2Box,
                                       nullptr, G4ThreeVector());
  auto lvQ2Coil = new G4LogicalVolume
    (solidQ2Coil, m_material_list->Cu, "lvQ2Coil");
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), pos - G4ThreeVector(0, 0, 335*mm)),
                    "pvQ2CoilU", lvQ2Coil, pMother, false, 0, m_check_overlaps);
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), pos + G4ThreeVector(0, 0, 335*mm)),
                    "pvQ2CoilD", lvQ2Coil, pMother, false, 1, m_check_overlaps);
  lvQ2Coil->SetVisAttributes(G4Color::Brown());
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructD1(G4VPhysicalVolume *pMother)
{
  G4double x, y, z;
  G4Material *PoleMater = m_material_list->Fe;
  //G4Material *D1GapMater = m_material_list->at("Vacuum");
  G4Material *D1GapMater = m_material_list->at("HeGas");
  // D magnet surface
  G4Tubs *solD1Tub = new G4Tubs("solD1Tub",
                                0,  Dfr2, DfHalfGap,
                                0.*degree, bendAngleDf*degree);
  // D magnet
  G4Tubs *solD1Gap = new G4Tubs("solD1Gap",
                                Dr1,  Dr2,  DHalfGap,
                                0.*degree, bendAngleD*degree);
  G4Box *tmpD1 = new G4Box("tmp", 1000, 1000, 1000);
  G4RotationMatrix tmpD1rot;
  tmpD1rot.rotateZ(35*deg);
  G4double tmpD1r = 1380-1000;
  G4ThreeVector tmpD13Y(tmpD1r*std::cos(35*deg),tmpD1r*std::sin(35*deg),0);
  ///// Yoke
  G4SubtractionSolid *tmpsolD1Pole
    = new G4SubtractionSolid("tmp", solD1Tub, solD1Gap);
  G4SubtractionSolid *solD1Pole
    = new G4SubtractionSolid("solidDy",
                             tmpsolD1Pole,
                             tmpD1,
                             G4Transform3D(tmpD1rot, tmpD13Y));
  G4LogicalVolume *lvD1Gap = new G4LogicalVolume(solD1Gap,
						  D1GapMater,
						  "lvD1Gap");
  G4LogicalVolume *lvD1Pole = new G4LogicalVolume(solD1Pole,
						   PoleMater,
						   "lvD1Pole");
  G4RotationMatrix rotD1;
  rotD1.rotateX(90.*deg);
  rotD1.rotateY(-110.*deg);
  G4ThreeVector pos(3*m, 0, -3*m*std::tan(35*deg));
  new G4PVPlacement(G4Transform3D(rotD1, pos),
                    "physD1Pole", lvD1Pole, pMother, false, 0, m_check_overlaps);
  new G4PVPlacement(G4Transform3D(rotD1, pos),
                    "physD1Gap", lvD1Gap, pMother, false, 0, m_check_overlaps);
  ///// Endguard
  G4Box *solidEG1 = new G4Box("solidEG1", 1880*mm/2., 1600*mm/2., 76*mm/2.);
  G4Box *solidEG2 = new G4Box("solidEG2", 810*mm/2., 330*mm/2., 80*mm/2.);
  auto solidD1EG = new G4SubtractionSolid("solEnd", solidEG1, solidEG2);
  auto lvD1EG = new G4LogicalVolume(solidD1EG, m_material_list->Fe, "lvD1EG");
  G4ThreeVector posEG(3*m*std::tan(35*deg) + 278.5*mm);
  posEG.rotateY(-20.*deg);
  G4RotationMatrix rotEG;
  rotEG.rotateY(70.*deg);
  new G4PVPlacement(G4Transform3D(rotEG, posEG),
        	    "pvD1EG", lvD1EG, pMother, false, 0, m_check_overlaps);
  lvD1Pole->SetVisAttributes(G4Color::Blue());
  lvD1Gap->SetVisAttributes(G4Color::Blue());
  lvD1EG->SetVisAttributes(G4Color::Blue());
  lvD1Gap->SetUserLimits(new G4UserLimits(3.*mm));
  ///// Coil
  G4VSolid* solidD1Coil = nullptr;
  solidD1Coil = new G4Box("solidD1Coil", 1600*mm/2, 1600*mm/2, 165*mm/2);
  auto solidD1Box = new G4Box("solidD1Gap", 800*mm/2, 500*mm/2, 900*mm/2);
  solidD1Coil = new G4SubtractionSolid("solidD1Coil", solidD1Coil, solidD1Box,
                                       nullptr, G4ThreeVector());
  auto lvD1Coil = new G4LogicalVolume
    (solidD1Coil, m_material_list->Cu, "lvD1Coil");
  G4ThreeVector posCoil(0, 0, -3*m*std::tan(35*deg)-165*mm/2);
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), posCoil),
                    "pvD1CoilU", lvD1Coil, pMother, false, 0, m_check_overlaps);
  posCoil.rotateY(-110*deg);
  G4RotationMatrix rotCoil;
  rotCoil.rotateY(-110*deg);
  new G4PVPlacement(G4Transform3D(rotCoil, posCoil),
                    "pvD1CoilD", lvD1Coil, pMother, false, 1, m_check_overlaps);
  lvD1Coil->SetVisAttributes(G4Color::Brown());
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructSDC1(G4VPhysicalVolume* pMother)
{
  auto sdSDC1 = new DCSD("SDC1");
  G4SDManager::GetSDMpointer()->AddNewDetector(sdSDC1);
  const auto& sdc1_pos = (geomMan.GetGlobalPosition("SDC1-V1") +
                          geomMan.GetGlobalPosition("SDC1-U2"))/2;
  const auto& frame_size = sizeMan.GetSize("Sdc1Frame")*mm/2;
  const auto& drift_size = sizeMan.GetSize("Sdc1Drift")*mm/2;
  auto sdc1_solid = new G4Box("Sdc1Solid", frame_size.x(),
                              frame_size.y(), frame_size.z());
  auto sdc1_lv = new G4LogicalVolume(sdc1_solid, m_material_list->Ar80IsoButane20Gas,
                                     "Sdc1LV", 0, 0, 0 );
  sdc1_lv->SetVisAttributes(G4Colour::Green());
  new G4PVPlacement(0, sdc1_pos,
                    "Sdc1PV", sdc1_lv, pMother, false, 0, m_check_overlaps);
  auto sdc1pl_solid = new G4Box("Sdc1PlSolid", drift_size.x(),
                                drift_size.y(), drift_size.z());
  G4String plane_name[] = { "Sdc1V1", "Sdc1V2", "Sdc1X1",
			    "Sdc1X2", "Sdc1U1", "Sdc1U2" };
  for(G4int i=0; i<NumOfLayersSDC1; ++i){
    G4ThreeVector pos;
    switch (i) {
    case 0:
      pos.setZ( -22.5985*mm );
      break;
    case 1:
      pos.setZ( -17.4015*mm );
      break;
    case 2:
      pos.setZ( -2.5985*mm );
      break;
    case 3:
      pos.setZ( 2.5985*mm );
      break;
    case 4:
      pos.setZ( 17.4015*mm );
      break;
    case 5:
      pos.setZ( 22.5985*mm );
      break;
    }
    auto sdc1pl_lv = new G4LogicalVolume(sdc1pl_solid,
                                         m_material_list->Ar80IsoButane20Gas,
                                         plane_name[i] + "LV", 0, 0, 0);
    sdc1pl_lv->SetSensitiveDetector(sdSDC1);
    new G4PVPlacement(nullptr, pos, sdc1pl_lv, plane_name[i] + "PV",
                      sdc1_lv, false, 101+i, m_check_overlaps);
  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructSDC2(G4VPhysicalVolume* pMother)
{
  auto sdSDC2 = new DCSD("SDC2");
  G4SDManager::GetSDMpointer()->AddNewDetector(sdSDC2);
  const auto& sdc2_pos = (geomMan.GetGlobalPosition("SDC2-V1") +
                          geomMan.GetGlobalPosition("SDC2-U2"))/2;
  const auto& frame_size = sizeMan.GetSize("Sdc2Frame")*mm/2;
  const auto& drift_size = sizeMan.GetSize("Sdc2Drift")*mm/2;
  auto sdc2_solid = new G4Box("Sdc2Solid", frame_size.x(),
                              frame_size.y(), frame_size.z());
  auto sdc2_lv = new G4LogicalVolume(sdc2_solid, m_material_list->Ar80IsoButane20Gas,
                                     "Sdc2LV", 0, 0, 0 );
  sdc2_lv->SetVisAttributes(G4Colour::Green());
  new G4PVPlacement(0, sdc2_pos,
                    "Sdc2PV", sdc2_lv, pMother, false, 0, m_check_overlaps);
  auto sdc2pl_solid = new G4Box("Sdc2PlSolid", drift_size.x(),
                                drift_size.y(), drift_size.z());
  G4String plane_name[] = { "Sdc2V1", "Sdc2V2",
                            "Sdc2U1", "Sdc2U2" };
  for(G4int i=0; i<NumOfLayersSDC2; ++i){
    G4ThreeVector pos;
    switch(i){
    case 0:
      pos.setZ(-12.4*mm);
      break;
    case 1:
      pos.setZ(-7.6*mm);
      break;
    case 2:
      pos.setZ(7.6*mm);
      break;
    case 3:
      pos.setZ(12.4*mm);
      break;
    }
    auto sdc2pl_lv = new G4LogicalVolume(sdc2pl_solid,
                                         m_material_list->Ar80IsoButane20Gas,
                                         plane_name[i] + "LV", 0, 0, 0);
    sdc2pl_lv->SetSensitiveDetector(sdSDC2);
    new G4PVPlacement(nullptr, pos, sdc2pl_lv, plane_name[i] + "PV",
                      sdc2_lv, false, 101+i, m_check_overlaps);
  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructSDC3(G4VPhysicalVolume* pMother)
{
  auto sdSDC3 = new DCSD("SDC3");
  G4SDManager::GetSDMpointer()->AddNewDetector(sdSDC3);
  const auto& sdc3_pos = (geomMan.GetGlobalPosition("SDC3-X1") +
                          geomMan.GetGlobalPosition("SDC3-Y2"))/2;
  const auto& frame_size = sizeMan.GetSize("Sdc3Frame")*mm/2;
  const auto& drift_size = sizeMan.GetSize("Sdc3Drift")*mm/2;
  const auto& ra2 = geomMan.GetRotAngle2("SDC3-X1")*deg;
  auto sdc3_solid = new G4Box("Sdc3Solid", frame_size.x(),
                              frame_size.y(), frame_size.z());
  auto sdc3_lv = new G4LogicalVolume(sdc3_solid, m_material_list->Ar80IsoButane20Gas,
                                     "Sdc3LV", 0, 0, 0 );
  sdc3_lv->SetVisAttributes(G4Colour::Green());
  auto rot = new G4RotationMatrix;
  rot->rotateY(-ra2);
  new G4PVPlacement(rot, sdc3_pos,
                    "Sdc3PV", sdc3_lv, pMother, false, 0, m_check_overlaps);
  auto sdc3pl_solid = new G4Box("Sdc3PlSolid", drift_size.x(),
                                drift_size.y(), drift_size.z());
  G4String plane_name[] = { "Sdc3X1", "Sdc3X2",
                            "Sdc3Y1", "Sdc3Y2" };
  for(G4int i=0; i<NumOfLayersSDC3; ++i){
    G4ThreeVector pos;
    switch(i){
    case 0:
      pos.setZ(-16.0*mm);
      break;
    case 1:
      pos.setZ(-8.206*mm);
      break;
    case 2:
      pos.setZ(8.206*mm);
      break;
    case 3:
      pos.setZ(16.0*mm);
      break;
    }
    auto sdc3pl_lv = new G4LogicalVolume(sdc3pl_solid,
                                         m_material_list->Ar80IsoButane20Gas,
                                         plane_name[i] + "LV", 0, 0, 0);
    sdc3pl_lv->SetSensitiveDetector(sdSDC3);
    new G4PVPlacement(nullptr, pos, sdc3pl_lv, plane_name[i] + "PV",
                      sdc3_lv, false, 101+i, m_check_overlaps);
  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructSDC4(G4VPhysicalVolume* pMother)
{
  auto sdSDC4 = new DCSD("SDC4");
  G4SDManager::GetSDMpointer()->AddNewDetector(sdSDC4);
  const auto& sdc4_pos = (geomMan.GetGlobalPosition("SDC4-Y1") +
                          geomMan.GetGlobalPosition("SDC4-X2"))/2;
  const auto& frame_size = sizeMan.GetSize("Sdc4Frame")*mm/2;
  const auto& drift_size = sizeMan.GetSize("Sdc4Drift")*mm/2;
  const auto& ra2 = geomMan.GetRotAngle2("SDC4-X1")*deg;
  auto sdc4_solid = new G4Box("Sdc4Solid", frame_size.x(),
                              frame_size.y(), frame_size.z());
  auto sdc4_lv = new G4LogicalVolume(sdc4_solid, m_material_list->Ar80IsoButane20Gas,
                                     "Sdc4LV", 0, 0, 0 );
  sdc4_lv->SetVisAttributes(G4Colour::Green());
  auto rot = new G4RotationMatrix;
  rot->rotateY(-ra2);
  new G4PVPlacement(rot, sdc4_pos,
                    "Sdc4PV", sdc4_lv, pMother, false, 0, m_check_overlaps);
  auto sdc4pl_solid = new G4Box("Sdc4PlSolid", drift_size.x(),
                                drift_size.y(), drift_size.z());
  G4String plane_name[] = { "Sdc4Y1", "Sdc4Y2",
                            "Sdc4X1", "Sdc4X2" };
  for(G4int i=0; i<NumOfLayersSDC4; ++i){
    G4ThreeVector pos;
    switch(i){
    case 0:
      pos.setZ(-16.0*mm);
      break;
    case 1:
      pos.setZ(-8.206*mm);
      break;
    case 2:
      pos.setZ(8.206*mm);
      break;
    case 3:
      pos.setZ(16.0*mm);
      break;
    }
    auto sdc4pl_lv = new G4LogicalVolume(sdc4pl_solid,
                                         m_material_list->Ar80IsoButane20Gas,
                                         plane_name[i] + "LV", 0, 0, 0);
    sdc4pl_lv->SetSensitiveDetector(sdSDC4);
    new G4PVPlacement(nullptr, pos, sdc4pl_lv, plane_name[i] + "PV",
                      sdc4_lv, false, 101+i, m_check_overlaps);
  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructSDC5(G4VPhysicalVolume* pMother)
{
  auto sdSDC5 = new DCSD("SDC5");
  G4SDManager::GetSDMpointer()->AddNewDetector(sdSDC5);
  const auto& sdc5_pos = (geomMan.GetGlobalPosition("SDC5-Y1") +
                          geomMan.GetGlobalPosition("SDC5-X2"))/2;
  const auto& frame_size = sizeMan.GetSize("Sdc5Frame")*mm/2;
  const auto& drift_size = sizeMan.GetSize("Sdc5Drift")*mm/2;
  const auto& ra2 = geomMan.GetRotAngle2("SDC5-X1")*deg;
  auto sdc5_solid = new G4Box("Sdc5Solid", frame_size.x(),
                              frame_size.y(), frame_size.z());
  auto sdc5_lv = new G4LogicalVolume(sdc5_solid, m_material_list->Ar80IsoButane20Gas,
                                     "Sdc5LV", 0, 0, 0 );
  sdc5_lv->SetVisAttributes(G4Colour::Green());
  auto rot = new G4RotationMatrix;
  rot->rotateY(-ra2);
  new G4PVPlacement(rot, sdc5_pos,
                    "Sdc5PV", sdc5_lv, pMother, false, 0, m_check_overlaps);
  auto sdc5pl_solid = new G4Box("Sdc5PlSolid", drift_size.x(),
                                drift_size.y(), drift_size.z());
  G4String plane_name[] = { "Sdc5Y1", "Sdc5Y2",
                            "Sdc5X1", "Sdc5X2" };
  for(G4int i=0; i<NumOfLayersSDC5; ++i){
    G4ThreeVector pos;
    switch(i){
    case 0:
      pos.setZ(-16.0*mm);
      break;
    case 1:
      pos.setZ(-8.206*mm);
      break;
    case 2:
      pos.setZ(8.206*mm);
      break;
    case 3:
      pos.setZ(16.0*mm);
      break;
    }
    auto sdc5pl_lv = new G4LogicalVolume(sdc5pl_solid,
                                         m_material_list->Ar80IsoButane20Gas,
                                         plane_name[i] + "LV", 0, 0, 0);
    sdc5pl_lv->SetSensitiveDetector(sdSDC5);
    new G4PVPlacement(nullptr, pos, sdc5pl_lv, plane_name[i] + "PV",
                      sdc5_lv, false, 101+i, m_check_overlaps);
  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructTOF(G4VPhysicalVolume* pMother)
{
  const auto& ra2 = geomMan.GetRotAngle2("TOF") * deg;
  const auto& half_size = sizeMan.GetSize("TofSeg")*mm/2;
  const G4double pitch = geomMan.GetWirePitch("TOF")*mm;
  auto sdTOF = new TOFSD("TOF");
  G4SDManager::GetSDMpointer()->AddNewDetector(sdTOF);
  // Mother
  auto mother_solid = new G4Box("TofMotherSolid",
                                half_size.x()*NumOfSegTOF + 50.*mm,
                                half_size.y() + 50.*mm,
                                half_size.z()*2 + 50.*mm);
  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       m_material_list->at("Air"),
                                       "TofMotherLV");
  auto rot = new G4RotationMatrix;
  rot->rotateY(-ra2);
  auto pos = geomMan.GetGlobalPosition("TOF");
  new G4PVPlacement(rot, pos, "TofMotherPV", mother_lv, pMother, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // Segment
  auto segment_solid = new G4Box("TofSegmentSolid", half_size.x(),
                                 half_size.y(), half_size.z());
  auto segment_lv = new G4LogicalVolume(segment_solid,
                                        m_material_list->Scin,
                                        "TofSegmentLV");
  for(G4int i=0; i<NumOfSegTOF; ++i){
    segment_lv->SetVisAttributes(G4Colour::Cyan());
    segment_lv->SetSensitiveDetector(sdTOF);
    pos = G4ThreeVector( ( -NumOfSegTOF/2 + i )*pitch,
			 0.0,
			 2.*( - i%2 + 0.5 )*half_size.z() );
    new G4PVPlacement(nullptr, pos, segment_lv,
                      "TofSegmentPV", mother_lv, false, i, m_check_overlaps);
  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructAC1(G4VPhysicalVolume* pMother)
{
  auto ac1_sd = new ACSD("AC1");
  // ac1_sd->SetRefractiveIndex( 1.05 );
  G4SDManager::GetSDMpointer()->AddNewDetector( ac1_sd );
  const auto& ra2 = geomMan.GetRotAngle2("AC1") * deg;
  const auto& frame_size = sizeMan.GetSize("Ac1Frame") * 0.5 * mm;
  const auto& radiator_size = sizeMan.GetSize("Ac1Radiator") * 0.5 * mm;
  // Mother
  auto mother_solid = new G4Box( "Ac1MotherSolid",
				 frame_size.x() + 5.*mm,
				 frame_size.y() + 5.*mm,
				 frame_size.z() + 5.*mm );
  auto mother_lv = new G4LogicalVolume( mother_solid,
					m_material_list->at("Air"),
					"Ac1MotherLV" );
  auto rot = new G4RotationMatrix;
  rot->rotateY(-ra2);
  auto pos = geomMan.GetGlobalPosition("AC1");
  G4ThreeVector offset( 0., 0., frame_size.z() - radiator_size.z());
  offset.rotateY(ra2);
  new G4PVPlacement(rot, pos + offset,
                    "Ac1MotherPV", mother_lv, pMother, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes( G4VisAttributes::GetInvisible());
  // Frame
  auto frame_solid = new G4Box( "Ac1FrameSolid", frame_size.x(),
				frame_size.y(), frame_size.z() );
  auto frame_lv = new G4LogicalVolume( frame_solid,
				       m_material_list->at("Air"),
				       "Ac1FrameLV" );
  pos.setMag( 0. );
  new G4PVPlacement( nullptr, pos, frame_lv,
		     "Ac1FramePV", mother_lv, false, 0 );
  // Radiator
  auto radiator_solid = new G4Box( "Ac1RadiatorSolid", radiator_size.x(),
				   radiator_size.y(), radiator_size.z() );
  auto radiator_lv = new G4LogicalVolume( radiator_solid,
					  m_material_list->Aerogel,
					  "Ac1RadiatorLV" );
  radiator_lv->SetSensitiveDetector( ac1_sd );
  radiator_lv->SetVisAttributes( G4Color::Magenta() );
  pos.set( 0., 0., -frame_size.z() + radiator_size.z());
  new G4PVPlacement( nullptr, pos, radiator_lv,
		     "Ac1RadiatorPV", frame_lv, false, 0 );
  // Mirror
  const G4double mirror_thickness = 1.*mm/2.;
  const G4double mirror_space = 20.*mm;
  const G4ThreeVector triangle_size( 1100*mm/2, frame_size.y(), 340.*mm );
  const G4double mirror_angle = std::atan2(triangle_size.z(),
                                           triangle_size.x());
  const G4ThreeVector mirror1_size( ( frame_size.x() - triangle_size.x() )/2.,
				    triangle_size.y(), mirror_thickness );
  const G4ThreeVector mirror2_size( std::hypot( triangle_size.x(),
						triangle_size.z() )/2.,
				    triangle_size.y(),
				    mirror_thickness );
  auto mirror1_solid = new G4Box( "Ac1Mirror1Solid", mirror1_size.x(),
				  mirror1_size.y(), mirror1_size.z() );
  auto mirror1_lv = new G4LogicalVolume( mirror1_solid,
					 m_material_list->Al,
					"Ac1Mirror1LV" );
  auto mirror2_solid = new G4Box( "Ac1Mirror2Solid", mirror2_size.x(),
				  mirror2_size.y(), mirror2_size.z() );
  auto mirror2_lv = new G4LogicalVolume( mirror2_solid,
					 m_material_list->Al,
					"Ac1Mirror2LV" );
  for( G4int i=0; i<2; ++i ){
    pos.set( ( triangle_size.x() + mirror1_size.x() ) * ( i*2 - 1 ),
	     0., frame_size.z() - mirror_space );
    new G4PVPlacement( nullptr, pos, mirror1_lv,
		       "Ac1MirrorPV", frame_lv, false, 0 );
    pos.set( triangle_size.x()/2 * ( i*2 - 1 ),
	     0., frame_size.z() - triangle_size.z()/2 - mirror_space );
    rot = new G4RotationMatrix;
    rot->rotateY( mirror_angle * ( i*2 - 1 ) );
    new G4PVPlacement( rot, pos, mirror2_lv,
		       "Ac1MirrorPV", frame_lv, false, 0 );
  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructWC(G4VPhysicalVolume* pMother)
{
  const auto& ra2 = geomMan.GetRotAngle2("WC") * deg;
  const auto& half_size_In = sizeMan.GetSize("WcSegIn") * 0.5 * mm;
  const auto& half_size_Out = sizeMan.GetSize("WcSegOut") * 0.5 * mm;
  const G4double pitch = geomMan.GetWirePitch("WC");
  auto wcSD = new WCSD("WC");
  // wcSD->SetRefractiveIndex( 1.33 );
  G4SDManager::GetSDMpointer()->AddNewDetector( wcSD );
  // Mother
  auto mother_solid = new G4Box( "WcMotherSolid",
				 half_size_Out.x()*NumOfSegWC + 200.*mm,
				 half_size_Out.y() + 200.*mm,
				 half_size_Out.z()*2 + 200.*mm );
				 // half_size_Out.x()*NumOfSegWC + 50.*mm,
				 // half_size_Out.y() + 50.*mm,
				 // half_size_Out.z()*2 + 50.*mm );

  auto mother_lv = new G4LogicalVolume( mother_solid,
					m_material_list->at("Air"),
					"WcMotherLV" );
  auto rot = new G4RotationMatrix;
  rot->rotateY(-ra2);
  auto pos = geomMan.GetGlobalPosition("WC");
  new G4PVPlacement(rot, pos,
                    "WcMotherPV",  mother_lv, pMother, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // Segment
  auto segment_solid = new G4Box( "WcSegmentSolid", half_size_In.x(),
				  half_size_In.y(), half_size_In.z() );
  auto segment_lv = new G4LogicalVolume( segment_solid,
					 m_material_list->Water,
					 "WcSegmentLV" );
  auto WCContainer     = new G4Box("WCContainer",
				   half_size_Out.x(),
				   half_size_Out.y(),
				   half_size_Out.z());
  auto WCContainer_gap = new G4Box("WCContainer_gap",
				   half_size_In.x(),
				   half_size_In.y(),
				   half_size_In.z());
  //G4RotationMatrix* rot_wccontainer_gap;
  auto rot_wccontainer_gap = new G4RotationMatrix;
  G4ThreeVector pos_wccontainer_gap(0.0, 0.0 ,0.0);
  auto solid_WCContainer
    = new G4SubtractionSolid("solid_WCContainer",
   			     WCContainer, WCContainer_gap,
   			     rot_wccontainer_gap,
			     pos_wccontainer_gap);
  auto logWCContainer = new G4LogicalVolume(solid_WCContainer,
					    m_material_list->Acrylic,
					    "logWCContainer");
  for(G4int i=0; i<NumOfSegWC; ++i){
    pos = G4ThreeVector( ( -NumOfSegWC/2 + i )*pitch,
			 0.0,
			 2.*( i%2 - 0.5 )*half_size_Out.z() );
    //for Vessel
    //    logWCContainer->SetVisAttributes( G4Colour::White() );
    logWCContainer->SetVisAttributes( G4Colour::Cyan() );
    new G4PVPlacement( nullptr, pos, logWCContainer,
		       "WcSegmentContainerPV", mother_lv, false, i );
    //for Water
    segment_lv->SetVisAttributes( G4Colour::Cyan() );
    segment_lv->SetSensitiveDetector( wcSD );
    new G4PVPlacement( nullptr, pos, segment_lv,
		       "WcSegmentPV", mother_lv, false, i );

  }
}

//_____________________________________________________________________________
void S2SDetectorConstruction::MakePositionDetector(G4VPhysicalVolume *pMother)
{
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  DCSD *dcSD = new DCSD("BcSD");
  SDMan->AddNewDetector(dcSD);

  //  bool flagInvisible = false;
  char name[100];

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  // ~~~~~~~ Material ~~~~~~~~~
  G4bool matflag = true; // Material flag for Chambers
  //G4bool matflag = false; // Material flag for Chambers
  G4bool heflag = true; // He bag flag
  //G4bool heflag = false; // He bag flag


  G4Material *DCuAreaMater = m_material_list->Ar80IsoButane20Gas;
  G4Material *DCAreaMater = m_material_list->Ar50Ethane50Gas;
  G4Material *DCLayerMater = m_material_list->at("Vacuum");//Ar80IsoButane20Gas;
  G4Material *DCuboxMater  = m_material_list->Mylar;
  G4Material *DCboxMater  = m_material_list->Mylar;
  if(matflag==false){
    DCuAreaMater = m_material_list->at("Vacuum");//Ar80IsoButane20Gas;
    DCAreaMater = m_material_list->at("Vacuum");//Ar50Ethane50Gas;
    DCuboxMater  = m_material_list->at("Vacuum");//Mylar;
    DCboxMater  = m_material_list->at("Vacuum");//Mylar;
  }


  G4double DCuAreaX = 1000.0;//mm
  G4double DCuAreaY = 1000.0;//mm
  G4double DCuAreaZ = 101.0;//mm
  G4double DCuBoxX = 1100.0;//mm
  G4double DCuBoxY = 1100.0;//mm
  G4double DCuBoxZ = 101.2;//mm

  //_____________________________________________________________________________//_____________________________________________________________________________//_____________________________________________________________________________//_____________________________________________________________________________/ BC1(x,x,v,v,u,u)
  G4Box *solidBc1Box =
    new G4Box( "BC1Box", DCuBoxZ/2.0*mm, DCuBoxX/2.0*mm, DCuBoxY/2.0*mm );
  G4Box *solidBc1Area =
    new G4Box( "BC1Area", DCuAreaZ/2.0*mm, DCuAreaX/2.0*mm, DCuAreaY/2.0*mm );
  G4SubtractionSolid *solidBc1Frame =
    new G4SubtractionSolid( "BC1Frame", solidBc1Box, solidBc1Area );
  G4Box *solidBc1Layer =
    new G4Box( "BC1Layer", 0.001/2.0*mm, DCuAreaX/2.0*mm, DCuAreaY/2.0*mm );

  G4LogicalVolume *logBc1Box =
    new G4LogicalVolume( solidBc1Box,   DCuboxMater,  "BC1Box",   0, 0, 0 );
  G4LogicalVolume *logBc1Area =
    new G4LogicalVolume( solidBc1Area,  DCuAreaMater,  "BC1Area",  0, 0, 0 );
  G4LogicalVolume *logBc1Frame =
    new G4LogicalVolume( solidBc1Frame,  DCuboxMater,  "BC1Frame",  0, 0, 0 );

  G4int lnum_bc1;
  lnum_bc1 = geomMan.GetDetectorId("BC1-u-1");
  G4ThreeVector posBc1G1 = geomMan.GetGlobalPosition( lnum_bc1 );
  lnum_bc1 = geomMan.GetDetectorId("BC1-x-2");
  G4ThreeVector posBc1G2 = geomMan.GetGlobalPosition( lnum_bc1 );
  G4ThreeVector posBc1G = (posBc1G1+posBc1G2)/2.0;

  G4double RotateAngleBc11;
  G4double RotateAngleBc12;
  RotateAngleBc11 = geomMan.GetRotAngle1( lnum_bc1 );
  RotateAngleBc12 = geomMan.GetRotAngle2( lnum_bc1 );
  G4RotationMatrix RMBc1;
  //  RMBc1.rotateX( RotateAngleBc11*degree );
  RMBc1.rotateZ( (RotateAngleBc12-90)*degree );

  //  G4RotationMatrix rotBc1;

  G4VPhysicalVolume *physBc1Box   = new G4PVPlacement( G4Transform3D(RMBc1, posBc1G),
						       "BC1Box",
						       logBc1Box,
						       pMother,
						       false,
						       0,
						       m_check_overlaps );
  G4VPhysicalVolume *physBc1Frame = new G4PVPlacement( 0,
						       //"BC1Box", logBc1Box, pMother, false, 0 );
						       G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
						       "BC1Frame",
						       logBc1Frame,
						       physBc1Box,
						       false,
						       0,
						       m_check_overlaps );
  G4VPhysicalVolume *physBc1Area = new G4PVPlacement( 0,
						      G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
						      "BC1Area",
						      logBc1Area,
						      // physBc1Frame,//physBc1Box,
						      physBc1Box,
						      false,
						      0,
						      m_check_overlaps );


  // G4String Bc1Name[6] = {"BC1-v-1", "BC1-v-2", "BC1-u-1", "BC1-u-2", "BC1-x-1", "BC1-x-2"};
  G4String Bc1Name[6] = {"BC1-x-1", "BC1-v-1", "BC1-u-1", "BC1-x-2", "BC1-v-2", "BC1-u-2"};
  G4LogicalVolume *logBc1Layer[6];
  G4VPhysicalVolume *physBc1Layer[6];
  for (int i=0; i<6; i++) {
    logBc1Layer[i] = new G4LogicalVolume( solidBc1Layer,
					  DCLayerMater,
					  Bc1Name[i]+"-LV",
					  0, 0, 0 );
    logBc1Layer[i]->SetSensitiveDetector(dcSD);

    lnum_bc1 = geomMan.GetDetectorId(Bc1Name[i]);
    G4ThreeVector posBc1L1 = geomMan.GetGlobalPosition( lnum_bc1 );
    G4double posBc1L1D = (posBc1L1-posBc1G).mag();

    if (i<3) posBc1L1D *= -1.;
    physBc1Layer[i] = new G4PVPlacement( 0,
					 G4ThreeVector( posBc1L1D*mm, 0.0*mm, 0.0*mm ),
					 Bc1Name[i]+"-PV",
					 logBc1Layer[i],
					 physBc1Area,
					 false,
					 101+i,
					 m_check_overlaps );
  }

  //logBc1Box->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logBc1Frame->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logBc1Area->SetVisAttributes(G4VisAttributes(true,G4Colour(0.3, 0.8, 0.5)));
  //logBc1Layer->SetVisAttributes(G4VisAttributes::GetInvisible());
  for(int i=0; i<6; i++) logBc1Layer[i]->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 1.0, 1.0)));

  //_____________________________________________________________________________//_____________________________________________________________________________//_____________________________________________________________________________//_____________________________________________________________________________/ BC2(x,x,v,v,u,u)
  G4Box *solidBc2Box =
    new G4Box( "BC2Box", DCuBoxZ/2.0*mm, DCuBoxX/2.0*mm, DCuBoxY/2.0*mm );
  G4Box *solidBc2Area =
    new G4Box( "BC2Area", DCuAreaZ/2.0*mm, DCuAreaX/2.0*mm, DCuAreaY/2.0*mm );
  G4SubtractionSolid *solidBc2Frame =
    new G4SubtractionSolid( "BC2Frame", solidBc2Box, solidBc2Area );
  G4Box *solidBc2Layer =
    new G4Box( "BC2Layer", 0.001/2.0*mm, DCuAreaX/2.0*mm, DCuAreaY/2.0*mm );

  G4LogicalVolume *logBc2Box =
    new G4LogicalVolume( solidBc2Box,   DCuboxMater,  "BC2Box",   0, 0, 0 );
  G4LogicalVolume *logBc2Area =
    new G4LogicalVolume( solidBc2Area,  DCuAreaMater,  "BC2Area",  0, 0, 0 );
  G4LogicalVolume *logBc2Frame =
    new G4LogicalVolume( solidBc2Frame, DCuboxMater,  "BC2Frame",  0, 0, 0 );

  G4int lnum_bc2;
  lnum_bc2 = geomMan.GetDetectorId("BC2-x-1");
  G4ThreeVector posBc2G1 = geomMan.GetGlobalPosition( lnum_bc2 );
  lnum_bc2 = geomMan.GetDetectorId("BC2-v-2");
  G4ThreeVector posBc2G2 = geomMan.GetGlobalPosition( lnum_bc2 );
  G4ThreeVector posBc2G = (posBc2G1+posBc2G2)/2.0;

  //G4cout<<"posBc2G="<<posBc2G<<G4endl;

  G4double RotateAngleBc21;
  G4double RotateAngleBc22;
  RotateAngleBc21 = geomMan.GetRotAngle1( lnum_bc2 );
  RotateAngleBc22 = geomMan.GetRotAngle2( lnum_bc2 );
  G4RotationMatrix RMBc2;
  //  RMBc2.rotateX( RotateAngleBc21*degree );
  RMBc2.rotateZ( (RotateAngleBc22-90)*degree );

  G4VPhysicalVolume *physBc2Box   = new G4PVPlacement( G4Transform3D(RMBc2, posBc2G),
						       "BC2Box",
						       logBc2Box,
						       pMother,
						       false, 0, m_check_overlaps );
  // G4VPhysicalVolume *physBc2Frame = new G4PVPlacement( 0,
  //       					       //"BC2Box",logBc2Box,
  //       					       G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
  //       					       "BC2Frame",
  //       					       logBc2Frame,
  //       					       physBc2Box,
  //       					       false, 0, m_check_overlaps );
  G4VPhysicalVolume *physBc2Area  = new G4PVPlacement( 0,
						       G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
						       "BC2Area",
						       logBc2Area,
						       physBc2Box,
						       false, 0, m_check_overlaps );

  G4String Bc2Name[4] = {"BC2-u-1", "BC2-x-1", "BC2-v-2", "BC2-u-2"};
  G4LogicalVolume *logBc2Layer[4];
  G4VPhysicalVolume *physBc2Layer[4];
  for (int i=0; i<4; i++) {
    logBc2Layer[i] = new G4LogicalVolume( solidBc2Layer,
					  DCLayerMater,
					  Bc2Name[i]+"-LV",
					  0, 0, 0 );
    logBc2Layer[i]->SetSensitiveDetector(dcSD);

    lnum_bc2 = geomMan.GetDetectorId(Bc2Name[i]);
    G4ThreeVector posBc2L1 = geomMan.GetGlobalPosition( lnum_bc2 );
    G4double posBc2L1D = (posBc2L1-posBc2G).mag();

    if (i<2) posBc2L1D *= -1.;
    physBc2Layer[i] = new G4PVPlacement( 0,
					 G4ThreeVector( posBc2L1D*mm, 0.0*mm, 0.0*mm ),
					 Bc2Name[i]+"-PV",
					 logBc2Layer[i],
					 physBc2Area,
					 false, 107+i, m_check_overlaps );
  }

  //logBc2Box->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logBc2Frame->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logBc2Area->SetVisAttributes(G4VisAttributes(true,G4Colour(0.3, 0.8, 0.5)));
  //logBc2Layer->SetVisAttributes(G4VisAttributes::GetInvisible());
  for(int i=0; i<4; i++) logBc2Layer[i]->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 1.0, 1.0)));


  G4double DCAreaX = 2000;//mm
  G4double DCAreaY = 1000;//mm
  G4double DCAreaZ = 100;//mm
  G4double DCBoxX = 2100;//mm
  G4double DCBoxY = 1100;//mm
  G4double DCBoxZ = 100.2;//mm

  //_____________________________________________________________________________//_____________________________________________________________________________//_____________________________________________________________________________//_____________________________________________________________________________/ BC3(x,x,y,y)
  G4Box *solidBc3Box =
    new G4Box( "BC3Box", DCBoxZ/2.0*mm, DCBoxX/2.0*mm, DCBoxY/2.0*mm );
  G4Box *solidBc3Area =
    new G4Box( "BC3Area", DCAreaZ/2.0*mm, DCAreaX/2.0*mm, DCAreaY/2.0*mm );
  G4SubtractionSolid *solidBc3Frame =
    new G4SubtractionSolid( "BC3Frame", solidBc3Box, solidBc3Area );
  G4Box *solidBc3Layer =
    new G4Box( "BC3Layer", 0.001/2.0*mm, DCAreaX/2.0*mm,  DCAreaY/2.0*mm );

  G4LogicalVolume *logBc3Box =
    new G4LogicalVolume( solidBc3Box,   DCboxMater,   "BC3Box",   0, 0, 0 );
  G4LogicalVolume *logBc3Area =
    new G4LogicalVolume( solidBc3Area,  DCAreaMater,  "BC3Area",  0, 0, 0 );
  G4LogicalVolume *logBc3Frame =
    new G4LogicalVolume( solidBc3Frame,  DCAreaMater,  "BC3Frame",  0, 0, 0 );

  G4int lnum_bc3;
  lnum_bc3 = geomMan.GetDetectorId("BC3-x-2");
  G4ThreeVector posBc3G1 = geomMan.GetGlobalPosition( lnum_bc3 );
//   posBc3G1 += G4ThreeVector(rhoD*tan(bendAngleD/2.*Deg2Rad), 0, 0);
//   posBc3G1.rotateZ(bendAngleD*degree);
  lnum_bc3 = geomMan.GetDetectorId("BC3-y-1");
  G4ThreeVector posBc3G2 = geomMan.GetGlobalPosition( lnum_bc3 );
//   posBc3G2 += G4ThreeVector(rhoD*tan(bendAngleD/2.*Deg2Rad), 0, 0);
//   posBc3G2.rotateZ(bendAngleD*degree);
  G4ThreeVector posBc3G = (posBc3G1+posBc3G2)/2.0;

  G4double RotateAngleBc31;
  G4double RotateAngleBc32;
  RotateAngleBc31 = geomMan.GetRotAngle1( lnum_bc3 );
  RotateAngleBc32 = geomMan.GetRotAngle2( lnum_bc3 );
  G4RotationMatrix RMBc3;
  //  RMBc3.rotateX( RotateAngleBc31*degree );
  RMBc3.rotateZ( (RotateAngleBc32-90)*degree );

  G4VPhysicalVolume *physBc3Box = new G4PVPlacement( G4Transform3D(RMBc3, posBc3G),
						     "BC3Box",
						     logBc3Box,
						     pMother,
						     false,
						     0,
						     m_check_overlaps );
  // G4VPhysicalVolume *physBc3Frame = new G4PVPlacement( G4Transform3D(RMBc3, posBc3G),
  // G4VPhysicalVolume *physBc3Frame = new G4PVPlacement( 0,
  //       					       G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
  //       					       //"BC3Box", logBc3Box,
  //       					       "BC3Frame",
  //       					       logBc3Frame,
  //       					       physBc3Box,
  //       					       false,
  //       					       0,
  //       					       m_check_overlaps );
  G4VPhysicalVolume *physBc3Area = new G4PVPlacement( 0,
						      G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
						      "BC3Area",
						      logBc3Area,
						      // physBc3Frame,// physBc3Box,
						      physBc3Box,
						      false,
						      0,
						      m_check_overlaps );

  G4String Bc3Name[4] = {"BC3-x-1", "BC3-x-2", "BC3-y-1", "BC3-y-2"};
  G4LogicalVolume *logBc3Layer[4];
  G4VPhysicalVolume *physBc3Layer[4];
  for (int i=0; i<4; i++) {
    logBc3Layer[i] = new G4LogicalVolume( solidBc3Layer,
					  DCLayerMater,
					  Bc3Name[i]+"-LV",
					  0, 0, 0 );
    logBc3Layer[i]->SetSensitiveDetector(dcSD);

    lnum_bc3 = geomMan.GetDetectorId(Bc3Name[i]);
    G4ThreeVector posBc3L1 = geomMan.GetGlobalPosition( lnum_bc3 );
    //     posBc3L1 += G4ThreeVector(rhoD*tan(bendAngleD/2.*Deg2Rad), 0, 0);
    //     posBc3L1.rotateZ(bendAngleD*degree);
    G4double posBc3L1D = (posBc3L1-posBc3G).mag();

    //G4double posBc3L1D = 10. + (i-3)*20.;

    if (i<2)
      posBc3L1D *= -1.;

    //    std::snprintf(name, sizeof(name), "BC3Layer%d", i+1);
    //    if(i==3){
    physBc3Layer[i] = new G4PVPlacement( 0,
  					 G4ThreeVector( posBc3L1D*mm, 0.0*mm, 0.0*mm ),
  					 Bc3Name[i]+"-PV",
  					 logBc3Layer[i],
  					 physBc3Area,
  					 false,
  					 111+i,
					 m_check_overlaps );
    //    }
  }

  // logBc3Box->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logBc3Frame->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logBc3Area->SetVisAttributes(G4VisAttributes(true,G4Colour(0.3, 0.8, 0.5)));
  // logBc3Area->SetVisAttributes(G4VisAttributes(true,G4Colour(1., 0., 0.)));
  for( int i = 1; i < 4; i++ ) logBc3Layer[i]->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 1.0, 1.0)));
  //logBc3Layer->SetVisAttributes(G4VisAttributes::GetInvisible());


  //_____________________________________________________________________________//_____________________________________________________________________________//_____________________________________________________________________________//_____________________________________________________________________________/ BC4(x,x,x,x,x,x)
  G4Box *solidBc4Box =
    new G4Box( "BC4Box", DCBoxZ/2.0*mm, DCBoxX/2.0*mm, DCBoxY/2.0*mm );
  G4Box *solidBc4Area =
    new G4Box( "BC4Area", DCAreaZ/2.0*mm, DCAreaX/2.0*mm, DCAreaY/2.0*mm );
  G4SubtractionSolid *solidBc4Frame =
    new G4SubtractionSolid( "BC4Frame", solidBc4Box, solidBc4Area );
  G4Box *solidBc4Layer =
    new G4Box( "BC4Layer", 0.001/2.0*mm, DCAreaX/2.0*mm,  DCAreaY/2.0*mm );

  G4LogicalVolume *logBc4Box =
    new G4LogicalVolume( solidBc4Box,   DCboxMater,  "BC4Box",   0, 0, 0 );
  G4LogicalVolume *logBc4Area =
    new G4LogicalVolume( solidBc4Area,  DCAreaMater,  "BC4Area",  0, 0, 0 );
  G4LogicalVolume *logBc4Frame =
    new G4LogicalVolume( solidBc4Frame,  DCAreaMater, "BC4Frame",  0, 0, 0 );

  G4int lnum_bc4;
  lnum_bc4 = geomMan.GetDetectorId("BC4-x-3");
  G4ThreeVector posBc4G1 = geomMan.GetGlobalPosition( lnum_bc4 );
//   posBc4G1 += G4ThreeVector(rhoD*tan(bendAngleD/2.*Deg2Rad), 0, 0);
//   posBc4G1.rotateZ(bendAngleD*degree);
  lnum_bc4 = geomMan.GetDetectorId("BC4-x-4");
  G4ThreeVector posBc4G2 = geomMan.GetGlobalPosition( lnum_bc4 );
//   posBc4G2 += G4ThreeVector(rhoD*tan(bendAngleD/2.*Deg2Rad), 0, 0);
//   posBc4G2.rotateZ(bendAngleD*degree);
  G4ThreeVector posBc4G = (posBc4G1+posBc4G2)/2.0;

  G4double RotateAngleBc41;
  G4double RotateAngleBc42;
  RotateAngleBc41 = geomMan.GetRotAngle1( lnum_bc4 );
  RotateAngleBc42 = geomMan.GetRotAngle2( lnum_bc4 );
  G4RotationMatrix RMBc4;
  //  RMBc4.rotateX( RotateAngleBc41*degree );
  RMBc4.rotateZ( (RotateAngleBc42-90)*degree );

  //  G4ThreeVector posBc4(3.6*m,0*m,0*m);
  //  posBc4G.rotateZ(bendAngleD*degree);

  //  G4RotationMatrix RMBc4;
  //  RMBc4.rotateZ(bendAngleD*degree);

  G4VPhysicalVolume *physBc4Box = new G4PVPlacement( G4Transform3D(RMBc4, posBc4G),
						     "BC4Box",
						     logBc4Box,
						     pMother,
						     false,
						     0,
						     m_check_overlaps );
//  G4VPhysicalVolume *physBc4Box = new G4PVPlacement( G4Transform3D(RMBc4, posBc4G),
//		       "BC4Box", logBc4Box, pMother, false, 0 );
  // G4VPhysicalVolume *physBc4Frame = new G4PVPlacement( 0,
  //       					       G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
  //       					       "BC4Frame",
  //       					       logBc4Frame,
  //       					       physBc4Box,
  //       					       false,
  //       					       0,
  //       					       m_check_overlaps );
  G4VPhysicalVolume *physBc4Area = new G4PVPlacement( 0,
						      G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
						      "BC4Area",
						      logBc4Area,
						      // physBc4Frame,//physBc4Box,
						      physBc4Box,
						      false,
						      0,
						      m_check_overlaps );

  G4String Bc4Name[4] = {"BC4-x-2", "BC4-x-3", "BC4-x-4", "BC4-x-5"};
  G4LogicalVolume *logBc4Layer[4];
  G4VPhysicalVolume *physBc4Layer[4];
  for (int i=0; i<4; i++) {
    logBc4Layer[i] = new G4LogicalVolume( solidBc4Layer,
					  DCLayerMater,
					  Bc4Name[i]+"-LV",
					  0, 0, 0 );
    logBc4Layer[i]->SetSensitiveDetector(dcSD);

    lnum_bc4 = geomMan.GetDetectorId(Bc4Name[i]);
    G4ThreeVector posBc4L1 = geomMan.GetGlobalPosition( lnum_bc4 );
//     posBc4L1 += G4ThreeVector(rhoD*tan(bendAngleD/2.*Deg2Rad), 0, 0);
//     posBc4L1.rotateZ(bendAngleD*degree);
    G4double posBc4L1D = (posBc4L1-posBc4G).mag();

    if (i<2)
      posBc4L1D *= -1.;

    //    if(i==99){
    physBc4Layer[i] = new G4PVPlacement( 0,
					 G4ThreeVector( posBc4L1D*mm, 0.0*mm, 0.0*mm ),
					 Bc4Name[i]+"-PV",
					 logBc4Layer[i],
					 physBc4Area,
					 false,
					 115+i,
					 m_check_overlaps );
      //    }
  }

  //logBc4Box->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logBc4Frame->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logBc4Area->SetVisAttributes(G4VisAttributes(true,G4Colour(0.3, 0.8, 0.5)));
  //  logBc4Layer->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  //logBc4Layer->SetVisAttributes(G4VisAttributes::GetInvisible());
  for(int i=1; i<4; i++) logBc4Layer[i]->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 1.0, 1.0)));



  //_____________________________________________________________________________//_____________________________________________________________________________//_____________________________________________________________________________//_____________________________________________________________________________/ Bc5(y,y,y,y,y,y)
  G4Box *solidBc5Box =
    new G4Box( "BC5Box", DCBoxZ/2.0*mm, DCBoxX/2.0*mm, DCBoxY/2.0*mm );
  G4Box *solidBc5Area =
    new G4Box( "BC5Area", DCAreaZ/2.0*mm, DCAreaX/2.0*mm, DCAreaY/2.0*mm );
  G4SubtractionSolid *solidBc5Frame =
    new G4SubtractionSolid( "BC5Frame", solidBc5Box, solidBc5Area );
  G4Box *solidBc5Layer =
    new G4Box( "BC5Layer", 0.001/2.0*mm, DCAreaX/2.0*mm,  DCAreaY/2.0*mm );

  G4LogicalVolume *logBc5Box =
    new G4LogicalVolume( solidBc5Box,   DCboxMater,  "BC5Box",   0, 0, 0 );
  G4LogicalVolume *logBc5Area =
    new G4LogicalVolume( solidBc5Area,  DCAreaMater,  "BC5Area",  0, 0, 0 );
  G4LogicalVolume *logBc5Frame =
    new G4LogicalVolume( solidBc5Frame,  DCAreaMater, "BC5Frame",  0, 0, 0 );

  G4int lnum_bc5;
  lnum_bc5 = geomMan.GetDetectorId("BC5-y-3");
  G4ThreeVector posBc5G1 = geomMan.GetGlobalPosition( lnum_bc5 );
//   posBc5G1 += G4ThreeVector(rhoD*tan(bendAngleD/2.*Deg2Rad), 0, 0);
//   posBc5G1.rotateZ(bendAngleD*degree);
  lnum_bc5 = geomMan.GetDetectorId("BC5-y-4");
  G4ThreeVector posBc5G2 = geomMan.GetGlobalPosition( lnum_bc5 );
//   posBc5G2 += G4ThreeVector(rhoD*tan(bendAngleD/2.*Deg2Rad), 0, 0);
//   posBc5G2.rotateZ(bendAngleD*degree);
  G4ThreeVector posBc5G = (posBc5G1+posBc5G2)/2.0;

  G4double RotateAngleBc51;
  G4double RotateAngleBc52;
  RotateAngleBc51 = geomMan.GetRotAngle1( lnum_bc5 );
  RotateAngleBc52 = geomMan.GetRotAngle2( lnum_bc5 );
  G4RotationMatrix RMBc5;
  //  RMBc5.rotateX( RotateAngleBc51*degree );
  RMBc5.rotateZ( (RotateAngleBc52-90)*degree );

  //  G4ThreeVector posBc5(3.6*m,0*m,0*m);
  //  posBc5G.rotateZ(bendAngleD*degree);

  //  G4RotationMatrix RMBc5;
  //  RMBc5.rotateZ(bendAngleD*degree);

  G4VPhysicalVolume *physBc5Box = new G4PVPlacement( G4Transform3D(RMBc5, posBc5G),
						     "BC5Box",
						     logBc5Box,
						     pMother,
						     false,
						     0,
						     m_check_overlaps );
  //G4VPhysicalVolume *physBc5Box = new G4PVPlacement( G4Transform3D(RMBc5, posBc5G), "BC5Box",logBc5Box,
  G4VPhysicalVolume *physBc5Frame = new G4PVPlacement( 0,
						       G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
						       "BC5Frame",
						       logBc5Frame,
						       physBc5Box,
						       false,
						       0,
						       m_check_overlaps);
  G4VPhysicalVolume *physBc5Area = new G4PVPlacement( 0,
						      G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
						      "BC5Area",
						      logBc5Area,
						      physBc5Box,
						      false,
						      0,
						      m_check_overlaps );

  G4String Bc5Name[4] = {"BC5-y-2", "BC5-y-3", "BC5-y-4", "BC5-y-5"};
  G4LogicalVolume *logBc5Layer[4];
  G4VPhysicalVolume *physBc5Layer[4];
  for (int i=0; i<4; i++) {
    logBc5Layer[i] = new G4LogicalVolume( solidBc5Layer,
					  DCLayerMater,
					  Bc5Name[i]+"-LV",
					  0, 0, 0 );
    logBc5Layer[i]->SetSensitiveDetector(dcSD);

    lnum_bc5 = geomMan.GetDetectorId(Bc5Name[i]);
    G4ThreeVector posBc5L1 = geomMan.GetGlobalPosition( lnum_bc5 );
//     posBc5L1 += G4ThreeVector(rhoD*tan(bendAngleD/2.*Deg2Rad), 0, 0);
//     posBc5L1.rotateZ(bendAngleD*degree);
    G4double posBc5L1D = (posBc5L1-posBc5G).mag();

    if (i<2)
      posBc5L1D *= -1.;

    //    if(i==5){
    physBc5Layer[i] = new G4PVPlacement( 0,
					 G4ThreeVector( posBc5L1D*mm, 0.0*mm, 0.0*mm ),
					 Bc5Name[i]+"-PV",
					 logBc5Layer[i],
					 physBc5Area,
					 false,
					 119+i,
					 m_check_overlaps );
    //    }
  }

  //logBc5Box->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logBc5Frame->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  logBc5Area->SetVisAttributes(G4VisAttributes(true,G4Colour(0.3, 0.8, 0.5)));
  //  logBc5Layer->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  //logBc5Layer->SetVisAttributes(G4VisAttributes::GetInvisible());
  for(int i=0; i<4; i++) logBc5Layer[i]->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 1.0, 1.0)));


  //_____________________________________________________________________________//_____________________________________________________________________________//_____________________________________________________________________________//_____________________________________________________________________________//_____________________________________________________________________________//_____________________________________________________________________________//_____________________________________________________________________________//_____________________________________________________________________________//_____________________________________________________________________________//_____________________________________________________________________________/
  //_____________________________________________________________________________//_____________________________________________________________________________////Other Tracking Detector

  //HeBag between
  //target and BC1 0
  //BC1 and BC2    1
  //BC2 and Q1     2
  //Q1 and Q2      3
  //Q2 and D1      4
  //D1 and BC3     5
  //BC3 and BC4    6
  G4Box *solidHeBag[7];
  G4LogicalVolume *logHeBag[7];
  // G4VPhysicalVolume *physHeBag[7];

  G4Material *HeBagMater[7];


  /*
  HeBagMater[7]={m_material_list->at("Vacuum"),//HeGas,//0
			     m_material_list->at("Vacuum"),//HeGas,//1
			     m_material_list->at("Vacuum"),//HeGas,//2
			     m_material_list->at("Vacuum"),//HeGas,//3
			     m_material_list->at("Vacuum"),//HeGas,//4
			     m_material_list->at("Vacuum"),//HeGas,//5
			     m_material_list->at("Vacuum")};//HeGas};//6
			     HeBagMater[7]={m_material_list->HeGas,//0
			     m_material_list->HeGas,//1
			     m_material_list->HeGas,//2
			     m_material_list->HeGas,//3
		 m_material_list->HeGas,//4
		 m_material_list->HeGas,//5
		 m_material_list->HeGas};//6
  */

  if(heflag==true){
    HeBagMater[0] = m_material_list->at("HeGas");
    HeBagMater[1] = m_material_list->at("HeGas");
    HeBagMater[2] = m_material_list->at("HeGas");
    HeBagMater[3] = m_material_list->at("HeGas");
    HeBagMater[4] = m_material_list->at("HeGas");
    HeBagMater[5] = m_material_list->at("HeGas");
    HeBagMater[6] = m_material_list->at("HeGas");
  }
  else{
    HeBagMater[0] = m_material_list->at("Vacuum");
    HeBagMater[1] = m_material_list->at("Vacuum");
    HeBagMater[2] = m_material_list->at("Vacuum");
    HeBagMater[3] = m_material_list->at("Vacuum");
    HeBagMater[4] = m_material_list->at("Vacuum");
    HeBagMater[5] = m_material_list->at("Vacuum");
    HeBagMater[6] = m_material_list->at("Vacuum");
  }

  G4double tgt=(rhoD*tan(bendAngleD/2.*Deg2Rad) + driftL2 + Q2z + driftL1 + Q1z + driftL0);
  G4double bc1=posBc1G.mag();
  G4double bc2=posBc2G.mag();
  G4double Q1=rhoD*tan(bendAngleD/2.*Deg2Rad) + driftL2 + Q2z + driftL1 + Q1z/2.;
  G4double Q2=rhoD*tan(bendAngleD/2.*Deg2Rad) + driftL2 + Q2z/2.;
  G4double bc3=posBc3G.mag();
  G4double bc4=posBc4G.mag();
  G4double D1=rhoD*tan(bendAngleD/2.*Deg2Rad);

  G4ThreeVector HeBagG[7];
  HeBagG[0]=G4ThreeVector(-(tgt+(bc1+DCuBoxZ/2.))/2.,0,0);
  HeBagG[1]=G4ThreeVector(-((bc1-DCuBoxZ/2.)+(bc2+DCuBoxZ/2.))/2.,0,0);
  HeBagG[2]=G4ThreeVector(-((bc2-DCuBoxZ/2.)+(Q1+Q1z/2.))/2.,0,0);
  HeBagG[3]=G4ThreeVector(-((Q1-Q1z/2.)+(Q2+Q2z/2.))/2.,0,0);
  HeBagG[4]=G4ThreeVector(-((Q2-Q2z/2.)+D1)/2.,0,0);
  HeBagG[5]=G4ThreeVector((D1+(bc3-DCBoxZ/2.))/2.,0,0);
  HeBagG[6]=G4ThreeVector(((bc3+DCBoxZ/2.)+(bc4-DCBoxZ/2.))/2.,0,0);
  G4double thick[7];
  thick[0]=tgt-(bc1+DCuBoxZ/2.);
  thick[1]=(bc1-DCuBoxZ/2.)-(bc2+DCuBoxZ/2.);
  thick[2]=(bc2-DCuBoxZ/2.)-(Q1+Q1z/2.);
  thick[3]=(Q1-Q1z/2.)-(Q2+Q2z/2.);
  thick[4]=(Q2-Q2z/2.)-rhoD*tan(bendAngleD/2.*Deg2Rad);
  thick[5]=-rhoD*tan(bendAngleD/2.*Deg2Rad)+(bc3-DCBoxZ/2.);
  thick[6]=-(bc3+DCBoxZ/2.)+(bc4-DCBoxZ/2.);

  HeBagG[5].rotateZ(bendAngleD*degree);
  HeBagG[6].rotateZ(bendAngleD*degree);
  G4RotationMatrix RMHeBag[7];
  RMHeBag[5].rotateZ(bendAngleD*degree);
  RMHeBag[6].rotateZ(bendAngleD*degree);

  char hename[7][100]={"HeBag0","HeBag1","HeBag2","HeBag3","HeBag4","HeBag5","HeBag6"};

  for(int i=0;i<7;i++){
    if(i!=0||1){
      if(i<2){
	solidHeBag[i] =  new G4Box( hename[i], thick[i]/2.,   DCuAreaX/2.*mm, DCuAreaY/2.*mm );
	// G4cout<<i<<" "<<thick[i]<<G4endl;
      }
      else{
	solidHeBag[i] =  new G4Box( hename[i], thick[i]/2.,   DCAreaX/2.*mm, DCAreaY/2.*mm );
	// G4cout<<i<<" "<<thick[i]<<G4endl;
      }
      logHeBag[i] = new G4LogicalVolume( solidHeBag[i], HeBagMater[i], hename[i], 0, 0, 0);

      // #############################################################################################
      //   Comment out, anayway (Toshi, 28Nov2014)
      // #############################################################################################
      // physHeBag[i] = new G4PVPlacement( G4Transform3D(RMHeBag[i], HeBagG[i]),
      // 					hename[i], logHeBag[i], pMother, false, 0 );
      //logHeBag[i]->SetVisAttributes(G4VisAttributes::GetInvisible());
      //logHeBag[i]->SetVisAttributes(G4VisAttributes(true,G4Colour(0.1, 0.1, 0.1)));
      // #############################################################################################
    }
  }

}//Make PositionDetectors

//_____________________________________________________________________________
void S2SDetectorConstruction::MakeSlits(G4VPhysicalVolume *pMother)
{

  G4Material *SlitMater = m_material_list->at("Vacuum");//Scin;

  //G4Box *solidSlit = new G4Box("solidSlit",4.0/2.*m, 4.0/2.*m, 0.000001/2.*mm); //Original
  //G4Box *solidSlit = new G4Box("solidSlit",1.5/2.*m, 1.5/2.*m, 0.000001/2.*mm); // Toshi , 25Nov2014
  //G4Box *solidSlit = new G4Box("solidSlit",1.5/2.*m, 0.8/2.*m, 0.000001/2.*mm); // Toshi , 25Nov2014
  G4Box *solidSlit = new G4Box("solidSlit",4.0/2.*m, 4.0/2.*m, 0.000001/2.*mm); // Toshi , 25Nov2014
  G4LogicalVolume *logicSlit = new G4LogicalVolume(solidSlit, SlitMater, "logicSlit");
  G4VPhysicalVolume *physSlit[11];

  G4RotationMatrix rotQ;
  rotQ.rotateZ(90.*deg);
  rotQ.rotateY(90.*deg);
  G4ThreeVector gloPosSlit[11];
  double x,y,z;
  y = 0;
  z = 0;

  x = -(rhoD*tan(bendAngleD/2.*Deg2Rad) + driftL2 + Q2z + driftL1 + Q1z/2.);
  gloPosSlit[0]=G4ThreeVector(x-Q1z/2.-5, y, z);
  gloPosSlit[1]=G4ThreeVector(x+Q1z/2.+5, y, z);

  x = -(rhoD*tan(bendAngleD/2.*Deg2Rad) + driftL2 + Q2z/2.);
  gloPosSlit[2]=G4ThreeVector(x-Q2z/2.-5, y, z);
  gloPosSlit[3]=G4ThreeVector(x+Q2z/2.+5, y, z);
  //gloPosSlit[2]=G4ThreeVector(x-Q2z/2.-50.0, y, z);
  //gloPosSlit[3]=G4ThreeVector(x+Q2z/2.+50.0, y, z);

  gloPosSlit[4]=G4ThreeVector(-rhoD*tan(bendAngleD/2.*Deg2Rad)-5, 0, 0);//D entrance
  gloPosSlit[5]=G4ThreeVector( rhoD*tan(bendAngleD/2.*Deg2Rad)+5, 0, 0);//D exit
  gloPosSlit[6]=G4ThreeVector( rhoD*tan(bendAngleD/2.*Deg2Rad)+240+76, 0, 0);//End guard

  // Order -->
  // vd8 | TOF | vd9 | AC | vd10 | WC | vd11 (vd[10])
  gloPosSlit[7]=G4ThreeVector( rhoD*tan(bendAngleD/2.*Deg2Rad)+500+1200, 0, 0);        // Toshi , 27Nov2014
  gloPosSlit[8]=G4ThreeVector( rhoD*tan(bendAngleD/2.*Deg2Rad)+500+1200+100, 0, 0);    // Toshi , 27Nov2014
  gloPosSlit[9]=G4ThreeVector( rhoD*tan(bendAngleD/2.*Deg2Rad)+500+1200+100+580, 0, 0);// Toshi , 27Nov2014
  gloPosSlit[10]=G4ThreeVector( rhoD*tan(bendAngleD/2.*Deg2Rad)+500+1000+1400+140, 0, 0);   // Toshi , 25Nov2014

  gloPosSlit[5].rotateZ(bendAngleD*Deg2Rad);
  gloPosSlit[6].rotateZ(bendAngleD*Deg2Rad);
  gloPosSlit[7].rotateZ(bendAngleD*Deg2Rad);
  gloPosSlit[8].rotateZ(bendAngleD*Deg2Rad);
  gloPosSlit[9].rotateZ(bendAngleD*Deg2Rad);
  gloPosSlit[10].rotateZ(bendAngleD*Deg2Rad);

  G4RotationMatrix rotSlit5;
  rotSlit5.rotateZ(90.*degree);
  rotSlit5.rotateY(90.*degree);
  rotSlit5.rotateZ((bendAngleD)*Deg2Rad);

  /*
    G4RotationMatrix rotSlit6;
    rotSlit6.rotateZ(90.*degree);
    rotSlit6.rotateY(90.*degree);
    rotSlit6.rotateZ((bendAngleD)*Deg2Rad);
    G4RotationMatrix rotSlit7;
    rotSlit7.rotateZ(90.*degree);
    rotSlit7.rotateY(90.*degree);
    rotSlit7.rotateZ((bendAngleD)*Deg2Rad);
    G4RotationMatrix rotSlit8;
    rotSlit8.rotateZ(90.*degree);
    rotSlit8.rotateY(90.*degree);
    rotSlit8.rotateZ((bendAngleD)*Deg2Rad);
  */

  char name[11][100]={"slitNo.1","slitNo.2","slitNo.3","slitNo.4",
		      "slitNo.5","slitNo.6","slitNo.7","slitNo.8",
		      "slitNo.9","slitNo.10","slitNo.11"};

//   physSlit[0] =
//     new G4PVPlacement( G4Transform3D(rotQ, gloPosSlit[0]),
// 		       "physSlit[0]", logicSlit, pMother, false, 0 );
  for(int i=0;i<5;i++){
    physSlit[i] =
      new G4PVPlacement( G4Transform3D(rotQ, gloPosSlit[i]),
			 name[i], logicSlit, pMother, false, i, m_check_overlaps );
  }

  for(int i=5;i<11;i++){
    physSlit[i] =
      new G4PVPlacement( G4Transform3D(rotSlit5, gloPosSlit[i]),
			 name[i], logicSlit, pMother, false, i, m_check_overlaps );
  }
  //   physSlit[6] =
//     new G4PVPlacement( G4Transform3D(rotSlit6, gloPosSlit[6]),
// 		       "physSlit[6]", logicSlit, pMother, false, 6 );
//   physSlit[7] =
//     new G4PVPlacement( G4Transform3D(rotSlit7, gloPosSlit[7]),
// 		       "physSlit[7]", logicSlit, pMother, false, 7 );

  // logicSlit->SetVisAttributes(G4Color::Gray());
  logicSlit->SetVisAttributes(G4VisAttributes::GetInvisible());

  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  SlSD *slSD = new SlSD("SlSD");
  SDMan->AddNewDetector(slSD);

  logicSlit->SetSensitiveDetector(slSD);
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::MakeTOFCounter(G4VPhysicalVolume *pMother)
{
  // const auto& ra2 = geomMan.GetRotAngle2("TOF") * CLHEP::deg;
  // const auto& half_size = sizeMan.GetSize("TofSeg") * 0.5 * mm;
  // const G4double pitch = geomMan.GetWirePitch("TOF") * mm;
  // auto tofSD = new TOFSD("TOF");
  // AddNewDetector(tofSD);
  // // Mother
  // auto mother_solid = new G4Box("TofMotherSolid",
  //                               half_size.x()*NumOfSegTOF + 50.*mm,
  //                               half_size.y() + 50.*mm,
  //                               half_size.z()*2 + 50.*mm);
  // auto mother_lv = new G4LogicalVolume(mother_solid,
  //                                      m_material_map["Air"],
  //                                      "FtofMotherLV");
  // auto rot = new G4RotationMatrix;
  // rot->rotateY(- ra2 - m_rotation_angle);
  // auto pos = (geomMan.GetGlobalPosition("KURAMA") +
  //             geomMan.GetGlobalPosition("TOF"));
  // pos.rotateY(m_rotation_angle);
  // new G4PVPlacement(rot, pos, mother_lv,
  //                   "FtofMotherPV", m_world_lv, false, 0, m_check_overlaps);
  // mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // // Segment
  // auto segment_solid = new G4Box("FtofSegmentSolid", half_size.x(),
  //                                half_size.y(), half_size.z());
  // auto segment_lv = new G4LogicalVolume(segment_solid,
  //                                       m_material_map["Scintillator"],
  //                                       "FtofSegmentLV");
  // for(G4int i=0; i<NumOfSegFTOF; ++i){
  //   segment_lv->SetVisAttributes(G4Colour::Cyan());
  //   segment_lv->SetSensitiveDetector(tofSD);
  //   pos = G4ThreeVector((-NumOfSegFTOF/2 + i)*pitch,
  //                       0.0,
  //                       2.*(- i%2 + 0.5)*half_size.z());
  //   new G4PVPlacement(nullptr, pos, segment_lv,
  //                     "FtofSegmentPV", mother_lv, false, i, m_check_overlaps);
  // }

  const int SegNumTOF = 18;
  double TOFX = 70.*mm;
  double TOFY = 20.*mm;
  double TOFZ = 600.*mm;

  // ~~~~~ Material ~~~~~~~~
  G4Material *TOFMatter    = m_material_list->Scin; // Plastic scintillator (original).
  //G4Material *TOFMatter    = m_material_list->at("Vacuum"); // Vacuum for test.

  G4ThreeVector localPosTOF[SegNumTOF];
  G4ThreeVector globalPosTOF[SegNumTOF];


  //G4int aaa = 0;
  //G4double overlap = 4.0;   // [mm] (default)
  //G4double distance = 10.0; // [mm] (default)
  // ~~~~~~~NO zig-zag~~~~~~~~~~~~~~~~~~~~~~
  //G4double overlap =  0.0;          // [mm] (gap 0mm)
  //G4double overlap = -1.0;          // [mm] (gap 1mm)
  //G4double overlap = -2.0;          // [mm] (gap 2mm)
  //G4double overlap = -3.0;          // [mm] (gap 3mm)
  //G4double overlap = -4.0;          // [mm] (gap 4mm)
  //G4double distance = -1.0*TOFY/mm; // [mm] (No zig-zag)
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  G4double overlap = confMan.GetTOFOverlap();
  G4double distance= confMan.GetTOFDistance();

  distance = distance + TOFY/mm;

  //localPosTOF[12] = G4ThreeVector(0, 0, 0);
  //globalPosTOF[12] = geomMan.Local2GlobalPos(61, localPosTOF[12]);
  //localPosTOF[8]  = G4ThreeVector(0, 0, -1.0*distance/2.0*mm);
  //localPosTOF[9]  = G4ThreeVector(TOFX/2.0 -overlap/2.0,
  localPosTOF[9]  = G4ThreeVector(0.0,
				  0.0,
				  -1.0*distance/2.0*mm);
  //localPosTOF[8]  = G4ThreeVector(0, 0, -2.0*cm);
  //globalPosTOF[9] = geomMan.Local2GlobalPos(61, localPosTOF[8]);
  globalPosTOF[9] = geomMan.Local2GlobalPos(61, localPosTOF[9]);

  //  G4double TOFpos = 3800.;
  //   localPosTOF[8] = G4ThreeVector(0, 0, 0);
  //   globalPosTOF[8] = G4ThreeVector(TOFpos, 0, 0);

  // ---- segment i=8-0 -----
  for (int i=SegNumTOF/2-1; i>=0; i--) {
    double localx = localPosTOF[i+1].x() - TOFX;
    //localx = localx + (overlap*aaa)*mm;
    localx = localx + overlap*mm;
    //aaa++;
    double localz = 0;
    if(i%2 == 0){
      //localz = -2.0*cm; //
      localz = 1.0*(distance)/2.0*mm;
    }
    else{
      //localz = 2.0*cm; //
      localz =  -1.0*(distance)/2.0*mm;
    }
    localPosTOF[i] = G4ThreeVector(localx, 0, localz);
    globalPosTOF[i] = geomMan.Local2GlobalPos(61, localPosTOF[i]);
    //globalPosTOF[i] = G4ThreeVector(TOFpos+localz, localx, 0);
  }

  // ---- segment i=10-18 -----
  //aaa = 0;
  //for (int i=SegNumTOF/2+1; i<=SegNumTOF; i++) {
  for (int i=SegNumTOF/2+1; i<SegNumTOF; i++) {
    double localx = localPosTOF[i-1].x() + TOFX;
    //localx = localx - (overlap*aaa)*mm;
    localx = localx - overlap*mm;
    //aaa++;
    double localz = 0;
    if(i%2 == 0){
      //localz = -2.0*cm;
      localz = 1.0*(distance)/2.0*mm;
    }
    else{
      //localz = 2.0*cm;
      localz = -1.0*(distance)/2.0*mm;
    }

    localPosTOF[i] = G4ThreeVector(localx, 0, localz);
    globalPosTOF[i] = geomMan.Local2GlobalPos(61, localPosTOF[i]);
    //globalPosTOF[i] = G4ThreeVector(TOFpos+localz, localx, 0);
  }

  G4double RotateAngleTOF = 70;//geomMan.GetRotAngle2( lnum_bh1 );
  G4RotationMatrix RMTOF;
  RMTOF.rotateZ( (RotateAngleTOF)*degree );

  G4Box *solidTOF = new G4Box( "solidTOF", TOFY/2.*mm,  TOFX/2.*mm,  TOFZ/2.*mm );
  G4LogicalVolume *logTOF = new G4LogicalVolume( solidTOF, TOFMatter, "logicTOF", 0, 0, 0 );
  G4VPhysicalVolume *physTOF[SegNumTOF];
  char buf[100];

  for(int i=0; i<SegNumTOF; i++) {

    std::snprintf(buf, sizeof(buf), "TOF_%d", i);
    //    globalPosTOF[i].rotateZ(70*degree);
    physTOF[i] =
      new G4PVPlacement( G4Transform3D(RMTOF, globalPosTOF[i]),
			 buf, logTOF, pMother, false, i, m_check_overlaps );
  }

  logTOF->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 1.0)));
  G4UserLimits* TOFLimit = new G4UserLimits(1.0 * mm);
  logTOF->SetUserLimits(TOFLimit);

  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  TOFSD *tofSD = new TOFSD("TOF");
  SDMan->AddNewDetector(tofSD);
  logTOF->SetSensitiveDetector(tofSD);

}


//_____________________________________________________________________________
void
S2SDetectorConstruction::MakeAerogelCounter(G4VPhysicalVolume* pMother)
{
  //G4Material *ACMatter = m_material_list->at("Vacuum");//Aerogel;
  G4Material *ACMatter    = m_material_list->Aerogel;
  G4Material* ACFrameMat  = m_material_list->at("Vacuum");

  // ~~~~ beam-direction, dispersion-direction, gravity-direction ~~~~~
  G4Box *solidACBox =
    new G4Box( "ACBox",   510.0/2.0*mm, 1200.0/2.0*mm, 600.0/2.0*mm );
  G4Box *solidACArea =
    new G4Box( "ACArea",  500.0/2.0*mm, 1100.0/2.0*mm, 550.0/2.0*mm );
  G4Box *solidACRad =
    new G4Box( "ACRad",   100.0/2.0*mm, 1100.0/2.0*mm, 550.0/2.0*mm );
  G4Box *solidACLayer =
    new G4Box( "ACLayer", 0.01/2.0*mm,  1100.0/2.0*mm, 550.0/2.0*mm );
  G4SubtractionSolid *solidACFrame =
    new G4SubtractionSolid( "ACFrame", solidACBox, solidACArea );

  // G4LogicalVolume *logACBox =
  //   new G4LogicalVolume( solidACBox,   m_material_list->at("Vacuum")/*Air*/,  "ACBox",   0, 0, 0 );
  G4LogicalVolume *logACFrame =
    new G4LogicalVolume( solidACFrame, ACFrameMat,  "ACFrame", 0, 0, 0 );
  G4LogicalVolume *logACArea =
    new G4LogicalVolume( solidACArea,  m_material_list->at("Air"),  "ACArea",  0, 0, 0 );
  G4LogicalVolume *logACRad =
    new G4LogicalVolume( solidACRad,   ACMatter,  "ACRad", 0, 0, 0 );
  G4LogicalVolume *logACLayer =
    new G4LogicalVolume( solidACLayer, ACMatter,  "ACLayer", 0, 0, 0 );

  /*
  G4int lnum_ac;
  lnum_ac = geomMan.GetDetectorId("AC-v-1");
  G4ThreeVector posACG1 = geomMan.GetGlobalPosition( lnum_ac );
  lnum_ac = geomMan.GetDetectorId("AC-v-2");
  G4ThreeVector posACG2 = geomMan.GetGlobalPosition( lnum_ac );
  G4ThreeVector posACG = (posACG1+posACG2)/2.0;

  G4double RotateAngleAC;
  RotateAngleAC = geomMan.GetRotAngle2( lnum_ac );

  G4RotationMatrix RMAC; RMAC.rotateZ( RotateAngleAC*degree );
  */

  G4ThreeVector posAC(4.2*m,0*m,0*m);
  posAC.rotateZ(bendAngleD*degree);
  G4RotationMatrix RMAC;
  RMAC.rotateZ(bendAngleD*degree);

  //G4VPhysicalVolume *physACBox =
  // G4VPhysicalVolume *physACFrame =
  //   new G4PVPlacement( G4Transform3D(RMAC,posAC),
  //       	       "ACFrame",
  //       	       logACFrame,//logACBox,
  //       	       pMother,
  //       	       false,
  //       	       0,
  //       	       m_check_overlaps );

//   G4VPhysicalVolume *physACFrame =
//     new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
//                        "ACFrame", logACFrame, physACBox, false, 0 );
  // G4VPhysicalVolume *physACArea =
  //   new G4PVPlacement( 0,
  //       	       G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
  //                      "ACArea",
  //       	       logACArea,
  //       	       physACFrame,//physACBox,
  //       	       false,
  //       	       0,
  //       	       m_check_overlaps );
  // G4VPhysicalVolume *physACRad =
  //   new G4PVPlacement( 0,
  //       	       G4ThreeVector( -200.0*mm, 0.0*mm, 0.0*mm ),
  //                      "ACRad",
  //       	       logACRad,
  //       	       physACArea, //physACFrame,
  //       	       false,
  //       	       0,
  //       	       m_check_overlaps );

  //  G4VPhysicalVolume *physACLayer =
  //new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
  //                     "ACLayer", logACLayer, physACRad, false, 0 );


  //char ACName[6]/*[100] = {*/ = "AC";//"AC-x-1", "AC-x-2", "AC-v-1", "AC-v-2", "AC-u-1", "AC-u-2"};

  //  G4VPhysicalVolume *physACLayer[6];
  /*
  for (int i=0; i<6; i++) {

    lnum_ac = geomMan.GetDetectorId(ACName[i]);
    G4ThreeVector posACL1 = geomMan.GetGlobalPosition( lnum_ac );
    G4double posACL1D = (posACL1-posACG).mag();

    //    G4double posACL1D = 0.5 + (i-3)*1.0;

    if (i>=3)
      posACL1D *= -1.;

    std::snprintf(name, sizeof(buf), "ACLayer%d", i+1);
    physACLayer[i] =
      new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, posACL1D*mm, 0.0*mm ),
			 name, logACLayer, physACArea, false,  113+i );
  }
  */

  G4VisAttributes* BoxAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
  BoxAtt -> SetForceWireframe(true);
  //logACBox ->SetVisAttributes(BoxAtt);
  logACFrame -> SetVisAttributes(BoxAtt);
  //logACBox->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 1.0)));
  //logACFrame->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0)));
  G4VisAttributes* AreaAtt = new G4VisAttributes(G4Colour(0.7, 0.0, 1.0));
  AreaAtt -> SetForceWireframe(true);
  logACArea ->SetVisAttributes(AreaAtt);
  //logACArea->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 1.0)));
  logACRad->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.0, 1.0)));
  // logACLayer->SetVisAttributes(G4VisAttributes::GetInvisible());

  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  ACSD *acSD = new ACSD("AC");
  SDMan->AddNewDetector(acSD);
  //logACLayer->SetSensitiveDetector(acSD);
  logACRad->SetSensitiveDetector(acSD);
}


//MakeWaterCherenkovCounter
void S2SDetectorConstruction::MakeWaterCounter(G4VPhysicalVolume* pMother){

  //const int SegNumWC = 16;
  /*
    double WCframeX = 80.*mm;
    double WCframeY = 156.*mm;
    double WCframeZ = 506.*mm;
    double WCradX = 74.*mm;
    double WCradY = 150.*mm;
    double WCradZ = 500.*mm;
  */

  // ~~~~~ Ver. 2014/11/25 ~~~~~~
  const int SegNumWC1 = 6; // upstream
  const int SegNumWC2 = 6; // downstream
  const int nseg = SegNumWC1 + SegNumWC2;
  const int SegNumWC  = nseg;
  double WCframeX = 180.0*mm;
  double WCframeY = 230.0*mm;
  double WCframeZ = 730.0*mm;
  //double WCradX = 150.0*mm;
  //double WCradY = 200.0*mm;
  //double WCradZ = 700.0*mm;
  double WCradX = 140.0*mm; // Toshi Gogami, 23July2015
  double WCradY = 190.0*mm; // Toshi Gogami, 23July2015
  double WCradZ = 690.0*mm; // Toshi Gogami, 23July2015

  double zdiff   = WCframeX/2.0;

  // ------ Material -------------
  // G4Material* WCMatter  = m_material_list->at("Vacuum");
  G4Material* RadMatter = m_material_list->Water; // (original)
  //G4Material* RadFrame  = m_material_list->Polyethylene;
  G4Material* RadFrame  = m_material_list->Acrylic; // (original)
  //G4Material* RadMatter = m_material_list->at("Vacuum");//Water;
  //G4Material* RadFrame  = m_material_list->at("Vacuum");//Acrylic;

  G4ThreeVector localPosWC[SegNumWC];
  G4ThreeVector globalPosWC[SegNumWC];

  //  localPosWC[5] = G4ThreeVector(0, 0, -zdiff);
  //  globalPosWC[5] = geomMan.Local2GlobalPos(lnum_bh1, localPosWC[5]);

  //G4double WCpos = 4.6*m;
  G4double WCpos = 4.7*m;
  //G4double WCtpos = -WCframeY/2.0;// (before 16Dec2014)
  G4double WCtpos = 0.0;// (after 16Dec2014)
  // ~~~~ WC1 ~~~~~
  G4int WCstart1 = 2;
  localPosWC[WCstart1] = G4ThreeVector(WCtpos-WCframeY/2.0, 0, 0);
  globalPosWC[WCstart1] = G4ThreeVector(WCpos-zdiff, WCtpos-WCframeY/2.0, 0);
  double localx=0.0 , localz=0.0;
  for(int i=WCstart1-1 ; i>=0 ; i--){
    //localx = localPosWC[i+1].x() - WCframeY/2;
    localx = 0.0;
    localz = 0.0;
    localx = localPosWC[i+1].x() - WCframeY;
    localz = -zdiff;
    localPosWC[i]  = G4ThreeVector(localx, 0, localz);
    globalPosWC[i] = G4ThreeVector(WCpos+localz, localx, 0);
  }
  for(int i=WCstart1+1 ; i<SegNumWC1 ; i++){
    //localx = localPosWC[i+1].x() - WCframeY/2;
    localx = 0.0;
    localz = 0.0;
    localx = localPosWC[i-1].x() + WCframeY;
    localz = -zdiff;
    localPosWC[i]  = G4ThreeVector(localx, 0, localz);
    globalPosWC[i] = G4ThreeVector(WCpos+localz, localx, 0);
  }
  // ~~~~ WC2 ~~~~~
  G4int WCstart2 = 8;
  //localPosWC[WCstart2] = G4ThreeVector(WCtpos-WCframeY/2.0, 0, 0);
  localPosWC[WCstart2] = G4ThreeVector(WCtpos, 0, 0);
  globalPosWC[WCstart2] = G4ThreeVector(WCpos+zdiff, WCtpos, 0);
  for(int i=WCstart2-1 ; i>=SegNumWC1 ; i--){
    localx = 0.0;
    localz = 0.0;
    localx = localPosWC[i+1].x() - WCframeY;
    localz = zdiff;
    localPosWC[i]  = G4ThreeVector(localx, 0, localz);
    globalPosWC[i] = G4ThreeVector(WCpos+localz, localx, 0);
  }
  for(int i=WCstart2+1 ; i<SegNumWC ; i++){
    localx = 0.0;
    localz = 0.0;
    localx = localPosWC[i-1].x() + WCframeY;
    localz = zdiff;
    localPosWC[i]  = G4ThreeVector(localx, 0, localz);
    globalPosWC[i] = G4ThreeVector(WCpos+localz, localx, 0);
  }

  /*
    G4int WCstart = 7;
    localPosWC[WCstart] = G4ThreeVector(WCtpos, 0, 0);
    globalPosWC[WCstart] = G4ThreeVector(WCpos+zdiff, WCtpos, 0);
    //for (int i=WCstart-1; i>=0; i--) {
    for (int i=WCstart-1; i>=0; i--) {
    double localx = localPosWC[i+1].x() - WCframeY/2;
    double localz;
    if(i%2 == 0)  localz = -zdiff;
    else  localz = zdiff;

    localPosWC[i] = G4ThreeVector(localx, 0, localz);
    //    globalPosWC[i] = geomMan.Local2GlobalPos(lnum_bh1, localPosWC[i]);
    globalPosWC[i] = G4ThreeVector(WCpos+localz, localx, 0);
    }

    for (int i=WCstart+1; i<=SegNumWC; i++) {
    double localx = localPosWC[i-1].x() + WCframeY/2;
    double localz;
    if (i%2 == 0)
    localz = -zdiff;
    else
    localz = zdiff;

    localPosWC[i] = G4ThreeVector(localx, 0, localz);
    //    globalPosWC[i] = geomMan.Local2GlobalPos(lnum_bh1, localPosWC[i]);
    globalPosWC[i] = G4ThreeVector(WCpos+localz, localx, 0);
    }
  */
  char buf[100];

  G4double RotateAngleWC = 70;//geomMan.GetRotAngle2( lnum_bh1 );
  G4RotationMatrix RMWC;
  RMWC.rotateZ( RotateAngleWC*degree );

  //G4Box *solidWCBox = new G4Box( "solidWCBox", WCframeX/2.*mm,  WCframeY/2.*mm,  WCframeZ/2.*mm );
  //G4Box *solidWCArea = new G4Box( "solidWCArea", WCradX/2.*mm,  WCradY/2.*mm,  WCradZ/2.*mm );
  //G4SubtractionSolid *solidWCFrame = new G4SubtractionSolid( "solidWCFrame", solidWCBox,  solidWCArea );
  //G4Box *solidWCLayer = new G4Box( "solidWCLayer", 0.01/2.*mm,  WCradY/2.*mm,  WCradZ/2.*mm );
  G4Box *solidWCBox = new G4Box( "solidWCBox", WCframeX/2.,  WCframeY/2.,  WCframeZ/2. );
  G4Box *solidWCArea = new G4Box( "solidWCArea", WCradX/2.,  WCradY/2.,  WCradZ/2. );
  G4SubtractionSolid *solidWCFrame = new G4SubtractionSolid( "solidWCFrame", solidWCBox,  solidWCArea );
  // G4Box *solidWCLayer = new G4Box( "solidWCLayer", 0.01/2.,  WCradY/2.,  WCradZ/2.);

  //G4LogicalVolume *logWCBox = new G4LogicalVolume( solidWCBox, WCMater, "logicWCBox", 0, 0, 0 );
  G4LogicalVolume *logWCArea = new G4LogicalVolume( solidWCArea, RadMatter, "logicWCArea", 0, 0, 0 );
  G4LogicalVolume *logWCFrame = new G4LogicalVolume( solidWCFrame, RadFrame, "logicWCFrame", 0, 0, 0 );
  // G4LogicalVolume *logWCLayer = new G4LogicalVolume( solidWCLayer, WCMatter, "logicWCLayer", 0, 0, 0 );

  //G4VPhysicalVolume *physWCBox[SegNumWC];
  G4VPhysicalVolume *physWCArea[SegNumWC];
  G4VPhysicalVolume *physWCFrame[SegNumWC];
  // G4VPhysicalVolume *physWCLayer[SegNumWC];

  for (int i=0; i<SegNumWC; i++) {

    globalPosWC[i].rotateZ(70*degree);
    std::snprintf(buf, sizeof(buf), "WCBox_%d", i);
    //physWCBox[i] =
    physWCFrame[i] = new G4PVPlacement( G4Transform3D(RMWC,globalPosWC[i]),
					//buf, logWCBox, pMother, false, 0 );
					buf,
					logWCFrame,
					pMother,
					false,
					i,
					m_check_overlaps );
    //     std::snprintf(buf, "WCFrame_%d", i);
    //     physWCFrame[i] =
    //       new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
    // 			 buf, logWCFrame, physWCBox[i], false, 0 );

    std::snprintf(buf, sizeof(buf), "WCArea_%d", i);
    physWCArea[i] = new G4PVPlacement( G4Transform3D(RMWC,globalPosWC[i]),//0,
				       //G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
				       buf,
				       logWCArea,
				       pMother,//physWCFrame[i],
				       false,
				       i,
				       m_check_overlaps );
    /*
      std::snprintf(buf, sizeof(buf), "WCLayer_%d", i);
      physWCLayer[i] = new G4PVPlacement( 0,
      //G4ThreeVector( -10.0*mm, 0.0*mm, 0.0*mm ),
      G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
      buf,
      logWCLayer,
      physWCArea[i],
      false,
      i );
    */

    G4VisAttributes* FrameAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    FrameAtt -> SetForceWireframe(true);
    //logWCBox->SetVisAttributes(G4VisAttributes::GetInvisible());
    //logWCFrame->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 1.0, 1.0)));
    //logWCFrame->SetVisAttributes(G4VisAttributes::SetForceWireframe(true));
    logWCFrame->SetVisAttributes(FrameAtt);
    //logWCFrame->SetVisAttributes(G4VisAttributes::GetInvisible());
    //logWCArea->SetVisAttributes(G4VisAttributes(true,G4Colour(0.0, 1.0, 1.0)));
    logWCArea->SetVisAttributes(G4VisAttributes(true,G4Colour(1.0, 0.5, 0.3)));
    //logWCArea->SetVisAttributes(G4VisAttributes::GetInvisible());
    //logWCLayer->SetVisAttributes(G4VisAttributes::GetInvisible());
    //G4VisAttributes::SetForceWireframe (G4bool force)
  }


  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  WCSD *wcSD = new WCSD("WCSD");
  SDMan->AddNewDetector(wcSD);
  //logWCLayer->SetSensitiveDetector(wcSD);
  logWCArea->SetSensitiveDetector(wcSD);

}//MakeWaterCounter


//_____________________________________________________________________________
G4bool S2SDetectorConstruction::
IsVolumeStopper( G4VPhysicalVolume *physVol ) const
{
  G4String name = physVol->GetName();
  if( name=="physQ1" || name=="physQ2" || name=="physDy" ) return true;
  else return false;
}
