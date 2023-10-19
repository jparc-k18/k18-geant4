// -*- C++ -*-

#include "S2SDetectorConstruction.hh"

#include <G4FieldManager.hh>
#include <G4ChordFinder.hh>
#include <G4TransportationManager.hh>

#include <G4Polyhedra.hh>
#include <G4Material.hh>
#include <G4Tubs.hh>
#include <G4EllipticalTube.hh>
#include <G4Box.hh>
#include <G4Trd.hh>
#include <G4ExtrudedSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4UnionSolid.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4PVReplica.hh>
#include <G4UniformMagField.hh>
#include <G4Mag_UsualEqRhs.hh>
#include <G4SimpleRunge.hh>
#include <G4VisAttributes.hh>

#include <G4GeometryManager.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4SolidStore.hh>
#include <G4UserLimits.hh>
#include <G4SystemOfUnits.hh>

#include <G4Colour.hh>
#include <G4SDManager.hh>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "DetSizeMan.hh"
#include "MagnetConstant.hh"
#include "MaterialList.hh"
#include "S2SField.hh"

#include "DCSD.hh"
#include "VPSD.hh"
#include "TOFSD.hh"
#include "ACSD.hh"
#include "WCSD.hh"
//#include "CalorimeterSD.hh"

namespace
{
const auto& confMan = ConfMan::GetInstance();
const auto& geomMan = DCGeomMan::GetInstance();
const auto& sizeMan = DetSizeMan::GetInstance();
const G4bool check_overlaps = true;
const G4bool use_hebag = true;
const MaterialList mlist;
G4PVPlacement* physWorld;
// G4Material* matlHeBag = mlist.at("HeGas");
}

//_____________________________________________________________________________
S2SDetectorConstruction::S2SDetectorConstruction()
{
}

//_____________________________________________________________________________
S2SDetectorConstruction::~S2SDetectorConstruction()
{
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
  auto solidWorld = new G4Box
    ("World", half_size.x(), half_size.y(), half_size.z());
  auto logicWorld = new G4LogicalVolume
    (solidWorld, mlist.at("Air"), "World");
  physWorld = new G4PVPlacement
    (nullptr, G4ThreeVector(), logicWorld, "World",
     nullptr, false, 0, check_overlaps);
  logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

#if 1
  ConstructTarget();
#endif

#if 1
  ConstructQ1();
  ConstructQ2();
  ConstructD1();
  ConstructField();
#endif

#if 1
  ConstructSDC1();
  ConstructSDC2();
#endif

#if 1
  ConstructSDC3();
  ConstructSDC4();
  ConstructSDC5();
#endif

#if 1
  ConstructTOF();
  ConstructAC1();
  ConstructWC();
#endif

#if 0
  ConstructVP();
#endif

  return physWorld;
}

//_____________________________________________________________________________
void S2SDetectorConstruction::ConstructField()
{
  return;
  S2SField *field = new S2SField(confMan.Get<G4String>("FLDMAP"));
  auto fieldManager =
    G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldManager->SetDetectorField(field);
  fieldManager->CreateChordFinder(field);

  // G4Mag_UsualEqRhs* fEquation = new G4Mag_UsualEqRhs(field);
  // G4MagIntegratorStepper *pStepper = new G4SimpleRunge(fEquation);
  // G4ChordFinder *pChordFinder = new G4ChordFinder(field, 1.0E-3*mm, pStepper);
  // fieldManager->SetChordFinder(pChordFinder);

  fieldManager->GetChordFinder()->SetDeltaChord(1.e-3*mm);
  //fieldManager->SetDeltaIntersection(1.0E-6*mm);
  //fieldManager->SetDeltaOneStep(1.0E-4*mm);
  //fieldManager->SetMaximumEpsilonStep(1.0E-3);
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructTarget()
{
  const auto& half_size = sizeMan.GetSize("Target")*mm/2.;
  G4Material *TargetMater = nullptr;
  auto Target = confMan.Get<G4String>("TargetMaterial");
  if(Target == "Be"){
    TargetMater = mlist.at("Be9");
  }else{
    G4cout << " Sorry, Target: " << Target
           << " is not defined. So Air will be used. " << G4endl;
    TargetMater = mlist.at("Air");
  }

  auto TargetBox = new G4Box
    ("TargetBox", half_size.x(), half_size.y(), half_size.z());
  auto logTarget = new G4LogicalVolume(TargetBox, TargetMater, "logTarget");
  G4RotationMatrix rotTarget;

  const auto& pos = geomMan.GetGlobalPosition("Target");
  new G4PVPlacement(G4Transform3D(rotTarget, pos),
                    "physTarget",
                    logTarget,
                    physWorld,
                    false,
                    0,
                    check_overlaps);
  logTarget->SetVisAttributes(G4Color::Gray());
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructQ1()
{
  const G4double a0 = 155*mm; // phi=310mm
  const G4ThreeVector pos(0, 0, -3885.35*mm);
  ///// Pole
  G4VSolid* solidQ1Gap;
  auto solidQ1Box = new G4Box("solidQ1Gap", 540*mm/2, 540*mm/2, 1240*mm/2);
  const G4double rCorner = (540*mm*std::sqrt(2)-2*a0)/2.;
  // auto solidCorner = new G4Tubs("solidCorner",
  //                               0*mm, rCorner, 900*mm,
  //                               phiStart, phiTotal);
  const G4double rInner[] = { 0., 0. };
  G4double zPlaneCorner[] = { -630*mm, 630*mm };
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
    (solidQ1Pole, mlist.at("Fe"), "lvQ1Pole");
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), pos),
                    "pvQ1Pole", lvQ1Pole, physWorld, false, 0, check_overlaps);
  lvQ1Pole->SetVisAttributes(G4Color::Cyan());
  ///// Coil
  G4VSolid* solidQ1Coil = nullptr;
  solidQ1Coil = new G4Box("solidQ1Coil", 1400*mm/2, 1400*mm/2, 180*mm/2);
  solidQ1Coil = new G4SubtractionSolid("solidQ1Coil", solidQ1Coil, solidQ1Box,
                                       nullptr, G4ThreeVector());
  auto lvQ1Coil = new G4LogicalVolume
    (solidQ1Coil, mlist.Cu, "lvQ1Coil");
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), pos - G4ThreeVector(0, 0, 530*mm)),
                    "pvQ1CoilU", lvQ1Coil, physWorld, false, 0, check_overlaps);
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), pos + G4ThreeVector(0, 0, 530*mm)),
                    "pvQ1CoilD", lvQ1Coil, physWorld, false, 1, check_overlaps);
  lvQ1Coil->SetVisAttributes(G4Color::Brown());
  ///// HeBag
  if(!use_hebag) return;
  auto lvQ1HeBag = new G4LogicalVolume(solidQ1Gap, mlist.at("HeGas"), "lvQ1HeBag");
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), pos),
                    "pvQ1HeBag", lvQ1HeBag, physWorld, false, 0, check_overlaps);
  const auto mylar_thickness = sizeMan.Get("HeBagMylarThickness")*mm;
  auto solidQ1Mylar = new G4Box("solidQ1Mylar", a0, a0, mylar_thickness/2);
  auto logicQ1Mylar = new G4LogicalVolume(solidQ1Mylar, mlist.at("Mylar"), "logicQ1Mylar");
  G4RotationMatrix rot;
  rot.rotateZ(45*deg);
  new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(0, 0, -620*mm+mylar_thickness/2)),
                    logicQ1Mylar, "physQ1Mylar", lvQ1HeBag, false, 0, check_overlaps);
  new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(0, 0,  620*mm-mylar_thickness/2)),
                    logicQ1Mylar, "physQ1Mylar", lvQ1HeBag, false, 1, check_overlaps);
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructQ2()
{
  const G4double a0 = 180*mm; // phi=360mm
  const G4ThreeVector pos(0, 0, -2776.5*mm);
  ///// Pole
  G4VSolid* solidQ2Pole = nullptr;
  solidQ2Pole = new G4Box("solidQ2Pole", 2100*mm/2, 1540*mm/2, 540*mm/2);
  const G4double zlength = 970*mm/2;
  auto solidQ2Gap = new G4Box("solidQ2Gap", a0, a0, zlength);
  // 1200x600 mm2;
  G4RotationMatrix rot;
  rot.rotateZ(45*deg);
  solidQ2Pole = new G4SubtractionSolid("solidQ2Pole", solidQ2Pole, solidQ2Gap,
                                       G4Transform3D(rot, G4ThreeVector()));
  auto lvQ2Pole = new G4LogicalVolume(solidQ2Pole, mlist.at("Fe"), "lvQ2Pole");
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), pos),
                    "pvQ2Pole", lvQ2Pole, physWorld, false, 0, check_overlaps);
  lvQ2Pole->SetVisAttributes(G4Color::Cyan());
  ///// Coil
  G4VSolid* solidQ2Coil = nullptr;
  solidQ2Coil = new G4Box("solidQ2Coil", 1320*mm/2, 880*mm/2, 130*mm/2);
  auto solidQ2Box = new G4Box("solidQ2Gap", 1200*mm/2, 600*mm/2, 900*mm/2);
  solidQ2Coil = new G4SubtractionSolid("solidQ2Coil", solidQ2Coil, solidQ2Box,
                                       nullptr, G4ThreeVector());
  auto lvQ2Coil = new G4LogicalVolume
    (solidQ2Coil, mlist.Cu, "lvQ2Coil");
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), pos - G4ThreeVector(0, 0, 335*mm)),
                    "pvQ2CoilU", lvQ2Coil, physWorld, false, 0, check_overlaps);
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), pos + G4ThreeVector(0, 0, 335*mm)),
                    "pvQ2CoilD", lvQ2Coil, physWorld, false, 1, check_overlaps);
  lvQ2Coil->SetVisAttributes(G4Color::Brown());
  ///// HeBag
  if(!use_hebag) return;
  auto lvQ2HeBag = new G4LogicalVolume(solidQ2Gap, mlist.at("HeGas"), "lvQ2HeBag");
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(rot), pos),
                    "pvQ2HeBag", lvQ2HeBag, physWorld, false, 0, check_overlaps);
  const auto mylar_thickness = sizeMan.Get("HeBagMylarThickness")*mm;
  auto solidQ2Mylar = new G4Box("solidQ2Mylar", a0, a0, mylar_thickness/2);
  auto logicQ2Mylar = new G4LogicalVolume(solidQ2Mylar, mlist.at("Mylar"), "logicQ2Mylar");
  rot.rotateZ(-45*deg);
  new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(0, 0, -zlength+mylar_thickness/2)),
                    logicQ2Mylar, "physQ2Mylar", lvQ2HeBag, false, 0, check_overlaps);
  new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(0, 0,  zlength-mylar_thickness/2)),
                    logicQ2Mylar, "physQ2Mylar", lvQ2HeBag, false, 1, check_overlaps);
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructD1()
{
  G4Material *mD1Gap = use_hebag ? mlist.at("HeGas") : mlist.at("Air");
  // D magnet surface
  G4Tubs *solD1Tub = new G4Tubs("solD1Tub", 0, Dfr2, DfHalfGap,
                                0.*degree, bendAngleDf*degree);
  // D magnet
  G4Tubs *solD1Gap = new G4Tubs("solD1Gap", Dr1, Dr2, DHalfGap,
                                0.*degree, bendAngleD*degree);
  G4Box *tmpD1 = new G4Box("tmp", 1000, 1000, 1000);
  G4RotationMatrix tmpD1rot;
  tmpD1rot.rotateZ(35*deg);
  G4double tmpD1r = 1380-1000;
  G4ThreeVector tmpD13Y(tmpD1r*std::cos(35*deg),tmpD1r*std::sin(35*deg),0);
  ///// Yoke
  auto tmpsolD1Pole = new G4SubtractionSolid("tmp", solD1Tub, solD1Gap);
  auto solD1Pole = new G4SubtractionSolid
    ("solidDy", tmpsolD1Pole, tmpD1, G4Transform3D(tmpD1rot, tmpD13Y));
  auto lvD1Gap = new G4LogicalVolume(solD1Gap, mD1Gap, "lvD1Gap");
  auto lvD1Pole = new G4LogicalVolume(solD1Pole, mlist.at("Fe"), "lvD1Pole");
  G4RotationMatrix rotD1;
  rotD1.rotateX(90.*deg);
  rotD1.rotateY(-110.*deg);
  G4ThreeVector pos(3*m, 0, -3*m*std::tan(35*deg));
  new G4PVPlacement(G4Transform3D(rotD1, pos),
                    "physD1Pole", lvD1Pole, physWorld, false, 0, check_overlaps);
  new G4PVPlacement(G4Transform3D(rotD1, pos),
                    "physD1Gap", lvD1Gap, physWorld, false, 0, check_overlaps);
  ///// Endguard
  G4Box *solidEG1 = new G4Box("solidEG1", 1880*mm/2., 1600*mm/2., 76*mm/2.);
  G4Box *solidEG2 = new G4Box("solidEG2", 810*mm/2., 330*mm/2., 80*mm/2.);
  auto solidD1EG = new G4SubtractionSolid("solEnd", solidEG1, solidEG2);
  auto lvD1EG = new G4LogicalVolume(solidD1EG, mlist.at("Fe"), "lvD1EG");
  G4ThreeVector posEG(3*m*std::tan(35*deg) + 278.5*mm);
  posEG.rotateY(-20.*deg);
  G4RotationMatrix rotEG;
  rotEG.rotateY(70.*deg);
  new G4PVPlacement(G4Transform3D(rotEG, posEG),
        	    "pvD1EG", lvD1EG, physWorld, false, 0, check_overlaps);
  lvD1Pole->SetVisAttributes(G4Color::Blue());
  // lvD1Gap->SetVisAttributes(G4Color::Blue());
  lvD1EG->SetVisAttributes(G4Color::Blue());
  lvD1Gap->SetUserLimits(new G4UserLimits(3.*mm));
  ///// Coil
  G4VSolid* solidD1Coil = nullptr;
  solidD1Coil = new G4Box("solidD1Coil", 1600*mm/2, 1600*mm/2, 165*mm/2);
  auto solidD1Box = new G4Box("solidD1Gap", 800*mm/2, 500*mm/2, 900*mm/2);
  solidD1Coil = new G4SubtractionSolid("solidD1Coil", solidD1Coil, solidD1Box,
                                       nullptr, G4ThreeVector());
  auto lvD1Coil = new G4LogicalVolume
    (solidD1Coil, mlist.Cu, "lvD1Coil");
  G4ThreeVector posCoil(0, 0, -3*m*std::tan(35*deg)-165*mm/2);
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(), posCoil),
                    "pvD1CoilU", lvD1Coil, physWorld, false, 0, check_overlaps);
  posCoil.rotateY(-110*deg);
  G4RotationMatrix rotCoil;
  rotCoil.rotateY(-110*deg);
  new G4PVPlacement(G4Transform3D(rotCoil, posCoil),
                    "pvD1CoilD", lvD1Coil, physWorld, false, 1, check_overlaps);
  lvD1Coil->SetVisAttributes(G4Color::Brown());
  ///// HeBag UpStream
  if(!use_hebag) return;
  const G4double zlengthU = 180*mm;
  auto solidHeBagU = new G4Box("solidD1HeBag", (Dr2-Dr1)/2, DHalfGap, zlengthU/2);
  auto lvD1HeBagU = new G4LogicalVolume(solidHeBagU, mlist.at("HeGas"), "lvD1HeBag");
  pos.setX(0);
  pos += G4ThreeVector(0, 0, -zlengthU/2);
  G4RotationMatrix rot;
  new G4PVPlacement(G4Transform3D(rot, pos),
                    "pvD1HeBagU", lvD1HeBagU, physWorld, false, 0, check_overlaps);
  const auto mylar_thickness = sizeMan.Get("HeBagMylarThickness")*mm;
  auto solidD1MylarU = new G4Box("solidD1Mylar", (Dr2-Dr1)/2, DHalfGap, mylar_thickness/2);
  auto logicD1MylarU = new G4LogicalVolume(solidD1MylarU, mlist.at("Mylar"), "logicD1MylarU");
  new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(0, 0, -zlengthU/2+mylar_thickness/2)),
                    logicD1MylarU, "physD1Mylar", lvD1HeBagU, false, 0, check_overlaps);
  ///// HeBag Downstream
  const G4double zlengthD = 316.5*mm;
  auto solidHeBagD = new G4Box("solidD1HeBag", (Dr2-Dr1)/2, DHalfGap, zlengthD/2);
  auto lvD1HeBagD = new G4LogicalVolume(solidHeBagD, mlist.at("HeGas"), "lvD1HeBag");
  pos.set(3*m*std::tan(35*deg) + zlengthD/2, 0, 0);
  pos.rotateY(-20.*deg);
  rot.rotateY(70.*deg);
  new G4PVPlacement(G4Transform3D(rot, pos),
                    "pvD1HeBagD", lvD1HeBagD, physWorld, false, 0, check_overlaps);
  auto solidD1MylarD = new G4Box("solidD1MylarD", (Dr2-Dr1)/2, DHalfGap, mylar_thickness/2);
  auto logicD1MylarD = new G4LogicalVolume(solidD1MylarD, mlist.at("Mylar"), "logicD1MylarD");
  rot.rotateY(-70.*deg);
  new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(0, 0, zlengthU/2-mylar_thickness/2)),
                    logicD1MylarD, "physD1MylarD", lvD1HeBagD, false, 1, check_overlaps);
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructSDC1()
{
  auto sdSDC1 = new DCSD("SDC");
  G4SDManager::GetSDMpointer()->AddNewDetector(sdSDC1);
  const auto& sdc1_pos = (geomMan.GetGlobalPosition("SDC1-V1") +
                          geomMan.GetGlobalPosition("SDC1-U2"))/2;
  const auto& frame_size = sizeMan.GetSize("Sdc1Frame")*mm/2;
  const auto& drift_size = sizeMan.GetSize("Sdc1Drift")*mm/2;
  auto sdc1_solid = new G4Box("Sdc1Solid", frame_size.x(),
                              frame_size.y(), frame_size.z());
  auto sdc1_lv = new G4LogicalVolume(sdc1_solid, mlist.Ar80IsoButane20Gas,
                                     "Sdc1LV", 0, 0, 0);
  sdc1_lv->SetVisAttributes(G4Colour::Green());
  new G4PVPlacement(nullptr, sdc1_pos,
                    "Sdc1PV", sdc1_lv, physWorld, false, 0, check_overlaps);
  auto sdc1pl_solid = new G4Box("Sdc1PlSolid", drift_size.x(),
                                drift_size.y(), drift_size.z());
  G4String plane_name[] = { "Sdc1V1", "Sdc1V2", "Sdc1X1",
			    "Sdc1X2", "Sdc1U1", "Sdc1U2" };
  for(G4int i=0; i<NumOfLayersSDC1; ++i){
    G4ThreeVector pos;
    switch (i) {
    case 0:
      pos.setZ(-22.5985*mm);
      break;
    case 1:
      pos.setZ(-17.4015*mm);
      break;
    case 2:
      pos.setZ(-2.5985*mm);
      break;
    case 3:
      pos.setZ(2.5985*mm);
      break;
    case 4:
      pos.setZ(17.4015*mm);
      break;
    case 5:
      pos.setZ(22.5985*mm);
      break;
    }
    auto sdc1pl_lv = new G4LogicalVolume(sdc1pl_solid,
                                         mlist.Ar80IsoButane20Gas,
                                         plane_name[i] + "LV", 0, 0, 0);
    sdc1pl_lv->SetSensitiveDetector(sdSDC1);
    new G4PVPlacement(nullptr, pos, sdc1pl_lv, plane_name[i] + "PV",
                      sdc1_lv, false, 101+i, check_overlaps);
  }
  ///// Mylar
  const G4double almylar_thickness = sizeMan.Get("AlMylarThickness")*mm;
  const G4double alplate_thickness = sizeMan.Get("AlPlateThickness")*mm;
  const G4double mylar_thickness = sizeMan.Get("Sdc1MylarThickness")*mm;
  auto solidMylar = new G4Box("solidMylar", frame_size.x(),
                                frame_size.y(), mylar_thickness/2);
  auto logicMylar = new G4LogicalVolume
    (solidMylar, mlist.at("Mylar"), "logicMylar");
  G4ThreeVector pos = sdc1_pos;
  pos.setZ(sdc1_pos.z()-frame_size.z()-mylar_thickness/2);
  new G4PVPlacement(nullptr, pos, "physMylar", logicMylar,
                    physWorld, false, 0, check_overlaps);
  pos.setZ(sdc1_pos.z()+frame_size.z()+mylar_thickness/2);
  new G4PVPlacement(nullptr, pos, "physMylar", logicMylar,
                    physWorld, false, 1, check_overlaps);
  ///// Al mylar
  auto solidAlMylar = new G4Box("solidAlMylar", frame_size.x(),
                                frame_size.y(), almylar_thickness/2);
  auto logicAlMylar = new G4LogicalVolume
    (solidAlMylar, mlist.at("Mylar"), "logicAlMylar");
  pos = sdc1_pos;
  pos.setZ(sdc1_pos.z()-frame_size.z()-mylar_thickness-almylar_thickness/2);
  new G4PVPlacement(nullptr, pos, "physAlMylar", logicAlMylar,
                    physWorld, false, 0, check_overlaps);
  pos.setZ(sdc1_pos.z()+frame_size.z()+mylar_thickness+almylar_thickness/2);
  new G4PVPlacement(nullptr, pos, "physAlMylar", logicAlMylar,
                    physWorld, false, 1, check_overlaps);
  ///// Al plate
  auto solidAlPlate = new G4Box("solidAlMylar", frame_size.x(),
                                frame_size.y(), alplate_thickness/2);
  auto logicAlPlate = new G4LogicalVolume
    (solidAlPlate, mlist.at("Al"), "logicAlPlate");
  pos = sdc1_pos;
  pos.setZ(sdc1_pos.z()-frame_size.z()-mylar_thickness-almylar_thickness-alplate_thickness/2);
  new G4PVPlacement(nullptr, pos, "physAlPlate", logicAlPlate,
                    physWorld, false, 0, check_overlaps);
  pos.setZ(sdc1_pos.z()+frame_size.z()+mylar_thickness+almylar_thickness+alplate_thickness/2);
  new G4PVPlacement(nullptr, pos, "physAlPlate", logicAlPlate,
                    physWorld, false, 1, check_overlaps);
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructSDC2()
{
  auto sdSDC2 = new DCSD("SDC2");
  G4SDManager::GetSDMpointer()->AddNewDetector(sdSDC2);
  const auto& sdc2_pos = (geomMan.GetGlobalPosition("SDC2-V1") +
                          geomMan.GetGlobalPosition("SDC2-U2"))/2;
  const auto& frame_size = sizeMan.GetSize("Sdc2Frame")*mm/2;
  const auto& drift_size = sizeMan.GetSize("Sdc2Drift")*mm/2;
  auto sdc2_solid = new G4Box("Sdc2Solid", frame_size.x(),
                              frame_size.y(), frame_size.z());
  auto sdc2_lv = new G4LogicalVolume(sdc2_solid, mlist.Ar80IsoButane20Gas,
                                     "Sdc2LV", 0, 0, 0);
  sdc2_lv->SetVisAttributes(G4Colour::Green());
  new G4PVPlacement(0, sdc2_pos,
                    "Sdc2PV", sdc2_lv, physWorld, false, 0, check_overlaps);
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
                                         mlist.Ar80IsoButane20Gas,
                                         plane_name[i] + "LV", 0, 0, 0);
    sdc2pl_lv->SetSensitiveDetector(sdSDC2);
    new G4PVPlacement(nullptr, pos, sdc2pl_lv, plane_name[i] + "PV",
                      sdc2_lv, false, 101+i, check_overlaps);
  }

  // ///// Mylar
  // const G4double almylar_thickness = sizeMan.Get("AlMylarThickness")*mm;
  // const G4double alplate_thickness = sizeMan.Get("AlPlateThickness")*mm;
  // const G4double mylar_thickness = sizeMan.Get("Sdc2MylarThickness")*mm;
  // auto solidMylar = new G4Box("solidMylar", frame_size.x(),
  //                               frame_size.y(), mylar_thickness/2);
  // auto logicMylar = new G4LogicalVolume
  //   (solidMylar, mlist.at("Mylar"), "logicMylar");
  // G4ThreeVector pos = sdc2_pos;
  // pos.setZ(sdc2_pos.z()-frame_size.z()-mylar_thickness/2);
  // new G4PVPlacement(nullptr, pos, "physMylar", logicMylar,
  //                   physWorld, false, 0, check_overlaps);
  // pos.setZ(sdc2_pos.z()+frame_size.z()+mylar_thickness/2);
  // new G4PVPlacement(nullptr, pos, "physMylar", logicMylar,
  //                   physWorld, false, 1, check_overlaps);
  // ///// Al mylar
  // auto solidAlMylar = new G4Box("solidAlMylar", frame_size.x(),
  //                               frame_size.y(), almylar_thickness/2);
  // auto logicAlMylar = new G4LogicalVolume
  //   (solidAlMylar, mlist.at("Mylar"), "logicAlMylar");
  // pos = sdc2_pos;
  // pos.setZ(sdc2_pos.z()-frame_size.z()-mylar_thickness-almylar_thickness/2);
  // new G4PVPlacement(nullptr, pos, "physAlMylar", logicAlMylar,
  //                   physWorld, false, 0, check_overlaps);
  // pos.setZ(sdc2_pos.z()+frame_size.z()+mylar_thickness+almylar_thickness/2);
  // new G4PVPlacement(nullptr, pos, "physAlMylar", logicAlMylar,
  //                   physWorld, false, 1, check_overlaps);
  // ///// Al plate
  // auto solidAlPlate = new G4Box("solidAlMylar", frame_size.x(),
  //                               frame_size.y(), alplate_thickness/2);
  // auto logicAlPlate = new G4LogicalVolume
  //   (solidAlPlate, mlist.at("Al"), "logicAlPlate");
  // pos = sdc2_pos;
  // pos.setZ(sdc2_pos.z()-frame_size.z()-mylar_thickness-almylar_thickness-alplate_thickness/2);
  // new G4PVPlacement(nullptr, pos, "physAlPlate", logicAlPlate,
  //                   physWorld, false, 0, check_overlaps);
  // pos.setZ(sdc2_pos.z()+frame_size.z()+mylar_thickness+almylar_thickness+alplate_thickness/2);
  // new G4PVPlacement(nullptr, pos, "physAlPlate", logicAlPlate,
  //                   physWorld, false, 1, check_overlaps);
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructSDC3()
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
  auto sdc3_lv = new G4LogicalVolume(sdc3_solid, mlist.Ar80IsoButane20Gas,
                                     "Sdc3LV", 0, 0, 0);
  sdc3_lv->SetVisAttributes(G4Colour::Green());
  auto rot = new G4RotationMatrix;
  rot->rotateY(-ra2);
  new G4PVPlacement(rot, sdc3_pos,
                    "Sdc3PV", sdc3_lv, physWorld, false, 0, check_overlaps);
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
                                         mlist.Ar80IsoButane20Gas,
                                         plane_name[i] + "LV", 0, 0, 0);
    sdc3pl_lv->SetSensitiveDetector(sdSDC3);
    new G4PVPlacement(nullptr, pos, sdc3pl_lv, plane_name[i] + "PV",
                      sdc3_lv, false, 101+i, check_overlaps);
  }
  ///// Mylar
  const G4double almylar_thickness = sizeMan.Get("AlMylarThickness")*mm;
  const G4double alplate_thickness = sizeMan.Get("AlPlateThickness")*mm;
  const G4double mylar_thickness = sizeMan.Get("Sdc3MylarThickness")*mm;
  auto solidMylar = new G4Box("solidMylar", frame_size.x(),
                                frame_size.y(), mylar_thickness/2);
  auto logicMylar = new G4LogicalVolume
    (solidMylar, mlist.at("Mylar"), "logicMylar");
  G4ThreeVector pos = sdc3_pos;
  pos.setZ(sdc3_pos.z()-frame_size.z()-mylar_thickness/2);
  new G4PVPlacement(rot, pos, "physMylar", logicMylar,
                    physWorld, false, 0, check_overlaps);
  pos.setZ(sdc3_pos.z()+frame_size.z()+mylar_thickness/2);
  new G4PVPlacement(rot, pos, "physMylar", logicMylar,
                    physWorld, false, 1, check_overlaps);
  ///// Al mylar
  auto solidAlMylar = new G4Box("solidAlMylar", frame_size.x(),
                                frame_size.y(), almylar_thickness/2);
  auto logicAlMylar = new G4LogicalVolume
    (solidAlMylar, mlist.at("Mylar"), "logicAlMylar");
  pos = sdc3_pos;
  pos.setZ(sdc3_pos.z()-frame_size.z()-mylar_thickness-almylar_thickness/2);
  new G4PVPlacement(rot, pos, "physAlMylar", logicAlMylar,
                    physWorld, false, 0, check_overlaps);
  pos.setZ(sdc3_pos.z()+frame_size.z()+mylar_thickness+almylar_thickness/2);
  new G4PVPlacement(rot, pos, "physAlMylar", logicAlMylar,
                    physWorld, false, 1, check_overlaps);
  ///// Al plate
  auto solidAlPlate = new G4Box("solidAlMylar", frame_size.x(),
                                frame_size.y(), alplate_thickness/2);
  auto logicAlPlate = new G4LogicalVolume
    (solidAlPlate, mlist.at("Al"), "logicAlPlate");
  pos = sdc3_pos;
  pos.setZ(sdc3_pos.z()-frame_size.z()-mylar_thickness-almylar_thickness-alplate_thickness/2);
  new G4PVPlacement(rot, pos, "physAlPlate", logicAlPlate,
                    physWorld, false, 0, check_overlaps);
  pos.setZ(sdc3_pos.z()+frame_size.z()+mylar_thickness+almylar_thickness+alplate_thickness/2);
  new G4PVPlacement(rot, pos, "physAlPlate", logicAlPlate,
                    physWorld, false, 1, check_overlaps);
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructSDC4()
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
  auto sdc4_lv = new G4LogicalVolume(sdc4_solid, mlist.Ar80IsoButane20Gas,
                                     "Sdc4LV", 0, 0, 0);
  sdc4_lv->SetVisAttributes(G4Colour::Green());
  auto rot = new G4RotationMatrix;
  rot->rotateY(-ra2);
  new G4PVPlacement(rot, sdc4_pos,
                    "Sdc4PV", sdc4_lv, physWorld, false, 0, check_overlaps);
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
                                         mlist.Ar80IsoButane20Gas,
                                         plane_name[i] + "LV", 0, 0, 0);
    sdc4pl_lv->SetSensitiveDetector(sdSDC4);
    new G4PVPlacement(nullptr, pos, sdc4pl_lv, plane_name[i] + "PV",
                      sdc4_lv, false, 101+i, check_overlaps);
  }
  ///// Mylar
  const G4double almylar_thickness = sizeMan.Get("AlMylarThickness")*mm;
  const G4double alplate_thickness = sizeMan.Get("AlPlateThickness")*mm;
  const G4double mylar_thickness = sizeMan.Get("Sdc4MylarThickness")*mm;
  auto solidMylar = new G4Box("solidMylar", frame_size.x(),
                                frame_size.y(), mylar_thickness/2);
  auto logicMylar = new G4LogicalVolume
    (solidMylar, mlist.at("Mylar"), "logicMylar");
  G4ThreeVector pos = sdc4_pos;
  pos.setZ(sdc4_pos.z()-frame_size.z()-mylar_thickness/2);
  new G4PVPlacement(rot, pos, "physMylar", logicMylar,
                    physWorld, false, 0, check_overlaps);
  pos.setZ(sdc4_pos.z()+frame_size.z()+mylar_thickness/2);
  new G4PVPlacement(rot, pos, "physMylar", logicMylar,
                    physWorld, false, 1, check_overlaps);
  ///// Al mylar
  auto solidAlMylar = new G4Box("solidAlMylar", frame_size.x(),
                                frame_size.y(), almylar_thickness/2);
  auto logicAlMylar = new G4LogicalVolume
    (solidAlMylar, mlist.at("Mylar"), "logicAlMylar");
  pos = sdc4_pos;
  pos.setZ(sdc4_pos.z()-frame_size.z()-mylar_thickness-almylar_thickness/2);
  new G4PVPlacement(rot, pos, "physAlMylar", logicAlMylar,
                    physWorld, false, 0, check_overlaps);
  pos.setZ(sdc4_pos.z()+frame_size.z()+mylar_thickness+almylar_thickness/2);
  new G4PVPlacement(rot, pos, "physAlMylar", logicAlMylar,
                    physWorld, false, 1, check_overlaps);
  ///// Al plate
  auto solidAlPlate = new G4Box("solidAlMylar", frame_size.x(),
                                frame_size.y(), alplate_thickness/2);
  auto logicAlPlate = new G4LogicalVolume
    (solidAlPlate, mlist.at("Al"), "logicAlPlate");
  pos = sdc4_pos;
  pos.setZ(sdc4_pos.z()-frame_size.z()-mylar_thickness-almylar_thickness-alplate_thickness/2);
  new G4PVPlacement(rot, pos, "physAlPlate", logicAlPlate,
                    physWorld, false, 0, check_overlaps);
  pos.setZ(sdc4_pos.z()+frame_size.z()+mylar_thickness+almylar_thickness+alplate_thickness/2);
  new G4PVPlacement(rot, pos, "physAlPlate", logicAlPlate,
                    physWorld, false, 1, check_overlaps);
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructSDC5()
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
  auto sdc5_lv = new G4LogicalVolume(sdc5_solid, mlist.Ar80IsoButane20Gas,
                                     "Sdc5LV", 0, 0, 0);
  sdc5_lv->SetVisAttributes(G4Colour::Green());
  auto rot = new G4RotationMatrix;
  rot->rotateY(-ra2);
  new G4PVPlacement(rot, sdc5_pos,
                    "Sdc5PV", sdc5_lv, physWorld, false, 0, check_overlaps);
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
                                         mlist.Ar80IsoButane20Gas,
                                         plane_name[i] + "LV", 0, 0, 0);
    sdc5pl_lv->SetSensitiveDetector(sdSDC5);
    new G4PVPlacement(nullptr, pos, sdc5pl_lv, plane_name[i] + "PV",
                      sdc5_lv, false, 101+i, check_overlaps);
  }
  ///// Mylar
  const G4double almylar_thickness = sizeMan.Get("AlMylarThickness")*mm;
  const G4double alplate_thickness = sizeMan.Get("AlPlateThickness")*mm;
  const G4double mylar_thickness = sizeMan.Get("Sdc5MylarThickness")*mm;
  auto solidMylar = new G4Box("solidMylar", frame_size.x(),
                                frame_size.y(), mylar_thickness/2);
  auto logicMylar = new G4LogicalVolume
    (solidMylar, mlist.at("Mylar"), "logicMylar");
  G4ThreeVector pos = sdc5_pos;
  pos.setZ(sdc5_pos.z()-frame_size.z()-mylar_thickness/2);
  new G4PVPlacement(rot, pos, "physMylar", logicMylar,
                    physWorld, false, 0, check_overlaps);
  pos.setZ(sdc5_pos.z()+frame_size.z()+mylar_thickness/2);
  new G4PVPlacement(rot, pos, "physMylar", logicMylar,
                    physWorld, false, 1, check_overlaps);
  ///// Al mylar
  auto solidAlMylar = new G4Box("solidAlMylar", frame_size.x(),
                                frame_size.y(), almylar_thickness/2);
  auto logicAlMylar = new G4LogicalVolume
    (solidAlMylar, mlist.at("Mylar"), "logicAlMylar");
  pos = sdc5_pos;
  pos.setZ(sdc5_pos.z()-frame_size.z()-mylar_thickness-almylar_thickness/2);
  new G4PVPlacement(rot, pos, "physAlMylar", logicAlMylar,
                    physWorld, false, 0, check_overlaps);
  pos.setZ(sdc5_pos.z()+frame_size.z()+mylar_thickness+almylar_thickness/2);
  new G4PVPlacement(rot, pos, "physAlMylar", logicAlMylar,
                    physWorld, false, 1, check_overlaps);
  ///// Al plate
  auto solidAlPlate = new G4Box("solidAlMylar", frame_size.x(),
                                frame_size.y(), alplate_thickness/2);
  auto logicAlPlate = new G4LogicalVolume
    (solidAlPlate, mlist.at("Al"), "logicAlPlate");
  pos = sdc5_pos;
  pos.setZ(sdc5_pos.z()-frame_size.z()-mylar_thickness-almylar_thickness-alplate_thickness/2);
  new G4PVPlacement(rot, pos, "physAlPlate", logicAlPlate,
                    physWorld, false, 0, check_overlaps);
  pos.setZ(sdc5_pos.z()+frame_size.z()+mylar_thickness+almylar_thickness+alplate_thickness/2);
  new G4PVPlacement(rot, pos, "physAlPlate", logicAlPlate,
                    physWorld, false, 1, check_overlaps);
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructTOF()
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
                                       mlist.at("Air"),
                                       "TofMotherLV");
  auto rot = new G4RotationMatrix;
  rot->rotateY(-ra2);
  auto pos = geomMan.GetGlobalPosition("TOF");
  new G4PVPlacement(rot, pos, "TofMotherPV", mother_lv, physWorld, false, 0, check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // Segment
  auto segment_solid = new G4Box("TofSegmentSolid", half_size.x(),
                                 half_size.y(), half_size.z());
  auto segment_lv = new G4LogicalVolume(segment_solid,
                                        mlist.Scin,
                                        "TofSegmentLV");
  for(G4int i=0; i<NumOfSegTOF; ++i){
    segment_lv->SetVisAttributes(G4Colour::Cyan());
    segment_lv->SetSensitiveDetector(sdTOF);
    pos = G4ThreeVector((-NumOfSegTOF/2 + i)*pitch,
			 0.0,
			 2.*(- i%2 + 0.5)*half_size.z());
    new G4PVPlacement(nullptr, pos, segment_lv,
                      "TofSegmentPV", mother_lv, false, i, check_overlaps);
  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructAC1()
{
  auto ac1_sd = new ACSD("AC1");
  // ac1_sd->SetRefractiveIndex(1.05);
  G4SDManager::GetSDMpointer()->AddNewDetector(ac1_sd);
  const auto& ra2 = geomMan.GetRotAngle2("AC1") * deg;
  const auto& frame_size = sizeMan.GetSize("Ac1Frame") * 0.5 * mm;
  const auto& radiator_size = sizeMan.GetSize("Ac1Radiator") * 0.5 * mm;
  // Mother
  auto mother_solid = new G4Box("Ac1MotherSolid",
				 frame_size.x() + 5.*mm,
				 frame_size.y() + 5.*mm,
				 frame_size.z() + 5.*mm);
  auto mother_lv = new G4LogicalVolume(mother_solid,
					mlist.at("Air"),
					"Ac1MotherLV");
  auto rot = new G4RotationMatrix;
  rot->rotateY(-ra2);
  auto pos = geomMan.GetGlobalPosition("AC1");
  G4ThreeVector offset(0., 0., frame_size.z() - radiator_size.z());
  offset.rotateY(ra2);
  new G4PVPlacement(rot, pos + offset,
                    "Ac1MotherPV", mother_lv, physWorld, false, 0, check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // Frame
  auto frame_solid = new G4Box("Ac1FrameSolid", frame_size.x(),
				frame_size.y(), frame_size.z());
  auto frame_lv = new G4LogicalVolume(frame_solid,
				       mlist.at("Air"),
				       "Ac1FrameLV");
  pos.setMag(0.);
  new G4PVPlacement(nullptr, pos, frame_lv,
		     "Ac1FramePV", mother_lv, false, 0);
  // Radiator
  auto radiator_solid = new G4Box("Ac1RadiatorSolid", radiator_size.x(),
				   radiator_size.y(), radiator_size.z());
  auto radiator_lv = new G4LogicalVolume(radiator_solid,
					  mlist.Aerogel,
					  "Ac1RadiatorLV");
  radiator_lv->SetSensitiveDetector(ac1_sd);
  radiator_lv->SetVisAttributes(G4Color::Magenta());
  pos.set(0., 0., -frame_size.z() + radiator_size.z());
  new G4PVPlacement(nullptr, pos, radiator_lv,
		     "Ac1RadiatorPV", frame_lv, false, 0);
  // Mirror
  const G4double mirror_thickness = 1.*mm/2.;
  const G4double mirror_space = 20.*mm;
  const G4ThreeVector triangle_size(1100*mm/2, frame_size.y(), 340.*mm);
  const G4double mirror_angle = std::atan2(triangle_size.z(),
                                           triangle_size.x());
  const G4ThreeVector mirror1_size((frame_size.x() - triangle_size.x())/2.,
				    triangle_size.y(), mirror_thickness);
  const G4ThreeVector mirror2_size(std::hypot(triangle_size.x(),
						triangle_size.z())/2.,
				    triangle_size.y(),
				    mirror_thickness);
  auto mirror1_solid = new G4Box("Ac1Mirror1Solid", mirror1_size.x(),
				  mirror1_size.y(), mirror1_size.z());
  auto mirror1_lv = new G4LogicalVolume(mirror1_solid,
                                        mlist.at("Al"),
					"Ac1Mirror1LV");
  auto mirror2_solid = new G4Box("Ac1Mirror2Solid", mirror2_size.x(),
				  mirror2_size.y(), mirror2_size.z());
  auto mirror2_lv = new G4LogicalVolume(mirror2_solid,
                                        mlist.at("Al"),
					"Ac1Mirror2LV");
  for(G4int i=0; i<2; ++i){
    pos.set((triangle_size.x() + mirror1_size.x()) * (i*2 - 1),
	     0., frame_size.z() - mirror_space);
    new G4PVPlacement(nullptr, pos, mirror1_lv,
		       "Ac1MirrorPV", frame_lv, false, 0);
    pos.set(triangle_size.x()/2 * (i*2 - 1),
	     0., frame_size.z() - triangle_size.z()/2 - mirror_space);
    rot = new G4RotationMatrix;
    rot->rotateY(mirror_angle * (i*2 - 1));
    new G4PVPlacement(rot, pos, mirror2_lv,
		       "Ac1MirrorPV", frame_lv, false, 0);
  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructWC()
{
  const auto& ra2 = geomMan.GetRotAngle2("WC") * deg;
  const auto& half_size_In = sizeMan.GetSize("WcSegIn") * 0.5 * mm;
  const auto& half_size_Out = sizeMan.GetSize("WcSegOut") * 0.5 * mm;
  const G4double pitch = geomMan.GetWirePitch("WC");
  auto wcSD = new WCSD("WC");
  // wcSD->SetRefractiveIndex(1.33);
  G4SDManager::GetSDMpointer()->AddNewDetector(wcSD);
  // Mother
  auto mother_solid = new G4Box("WcMotherSolid",
				 half_size_Out.x()*NumOfSegWC + 200.*mm,
				 half_size_Out.y() + 200.*mm,
				 half_size_Out.z()*2 + 200.*mm);
				 // half_size_Out.x()*NumOfSegWC + 50.*mm,
				 // half_size_Out.y() + 50.*mm,
				 // half_size_Out.z()*2 + 50.*mm);

  auto mother_lv = new G4LogicalVolume(mother_solid,
					mlist.at("Air"),
					"WcMotherLV");
  auto rot = new G4RotationMatrix;
  rot->rotateY(-ra2);
  auto pos = geomMan.GetGlobalPosition("WC");
  new G4PVPlacement(rot, pos,
                    "WcMotherPV",  mother_lv, physWorld, false, 0, check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // Segment
  auto segment_solid = new G4Box("WcSegmentSolid", half_size_In.x(),
				  half_size_In.y(), half_size_In.z());
  auto segment_lv = new G4LogicalVolume(segment_solid,
					 mlist.Water,
					 "WcSegmentLV");
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
					    mlist.Acrylic,
					    "logWCContainer");
  for(G4int i=0; i<NumOfSegWC; ++i){
    pos = G4ThreeVector((-NumOfSegWC/2 + i)*pitch,
			 0.0,
			 2.*(i%2 - 0.5)*half_size_Out.z());
    //for Vessel
    //    logWCContainer->SetVisAttributes(G4Colour::White());
    logWCContainer->SetVisAttributes(G4Colour::Cyan());
    new G4PVPlacement(nullptr, pos, logWCContainer,
		       "WcSegmentContainerPV", mother_lv, false, i);
    //for Water
    segment_lv->SetVisAttributes(G4Colour::Cyan());
    segment_lv->SetSensitiveDetector(wcSD);
    new G4PVPlacement(nullptr, pos, segment_lv,
		       "WcSegmentPV", mother_lv, false, i);

  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructVP()
{
  G4Material *VPMater = mlist.at("Air");
  auto solidVP = new G4Box("solidVP", 4*m/2, 4*m/2, 1*CLHEP::um/2);
  auto logicVP = new G4LogicalVolume(solidVP, VPMater, "logicVP");
  auto sdVP = new VPSD("VP");
  G4SDManager::GetSDMpointer()->AddNewDetector(sdVP);
  logicVP->SetSensitiveDetector(sdVP);
  logicVP->SetVisAttributes(G4VisAttributes::GetInvisible());
  for(G4int i=0; i<NumOfLayersVP; ++i){
    auto pos = geomMan.GetGlobalPosition("VP"+std::to_string(i+1));
    auto rot = new G4RotationMatrix;
    rot->rotateY(-geomMan.GetRotAngle2("VP"+std::to_string(i+1))*deg);
    new G4PVPlacement(rot, pos, "pvVP", logicVP, physWorld, false, i,
                      check_overlaps);
  }
}

//_____________________________________________________________________________
G4bool
S2SDetectorConstruction::IsVolumeStopper(G4VPhysicalVolume *physVol) const
{
  G4String name = physVol->GetName();
  if(name=="physQ1" || name=="physQ2" || name=="physDy") return true;
  else return false;
}
