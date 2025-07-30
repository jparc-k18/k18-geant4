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
#include <G4Hype.hh>
#include <G4IntersectionSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4UnionSolid.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4PVReplica.hh>
#include <G4UniformMagField.hh>
#include <G4Mag_UsualEqRhs.hh>
#include <G4SimpleRunge.hh>
#include <G4VisAttributes.hh>
#include <G4GenericTrap.hh>

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
#include "MaterialList.hh"

#include "RCSD.hh"
#include "PDSD.hh"

namespace
{
  const auto& gConf = ConfMan::GetInstance();
  const auto& gGeom = DCGeomMan::GetInstance();
  const auto& gSize = DetSizeMan::GetInstance();
  constexpr G4double RadToDeg = 180./M_PI;
  const G4Colour ORANGE(1.0, 0.55, 0.0);

  inline G4double square(G4double x) { return x*x; }
  inline G4Material* mat(const G4String& name) {
    return G4Material::GetMaterial(name, /*warning*/ true);
  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructRC(G4int i)
{
  const G4String name = "RC"+std::to_string(i);
  const auto& ra2 = gGeom.GetRotAngle2(name) * deg;
  const auto& half_size = gSize.GetSize("RcSeg")*mm/2;
  const G4double pitch = gGeom.GetWirePitch(name)*mm;
  auto sdRC = new RCSD(name);
  AddNewDetector(sdRC);
  // Mother
  auto mother_solid = new G4Box("RcMotherSolid",
				half_size.x()*mm,
                                half_size.y()*mm,
                                half_size.z()*NumOfSegRC*mm);
  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       mat("Air"),
                                       "RcMotherLV");
  auto rot = new G4RotationMatrix;
  rot->rotateY(-ra2);
  auto pos = gGeom.GetGlobalPosition(name);
  new G4PVPlacement(rot, pos, mother_lv, "RcMotherPHYS", m_world_lv, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // Segment
  auto segment_solid = new G4Box("RcSegmentSolid", half_size.x(),
                                 half_size.y(), half_size.z());
  auto segment_lv = new G4LogicalVolume(segment_solid,
                                        mat("Scintillator"),
                                        "RcSegmentLV");
  for(G4int i=0; i<NumOfSegRC; ++i){
    segment_lv->SetVisAttributes(G4Colour::Cyan());
    segment_lv->SetSensitiveDetector(sdRC);
    pos = G4ThreeVector(0.0,
                        0.0,
			(-NumOfSegRC/2 + i)*pitch + pitch/2);
    new G4PVPlacement(nullptr, pos, segment_lv,
                      "RcSegmentPHYS", mother_lv, false, i, m_check_overlaps);
  }
}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructPD(G4int i, G4String YZ)
{
  const G4String name = "PD"+std::to_string(i)+YZ;
  const G4String RCname = "RC"+std::to_string(i);
  const auto& ra2 = gGeom.GetRotAngle2(RCname) * deg;
  const auto& half_size = gSize.GetSize("Pd"+YZ)*mm/2;
  const auto& seg_size = gSize.GetSize("Pd"+YZ+"Seg")*mm/2;
  const auto& Rcseg_size = gSize.GetSize("RcSeg")*mm/2;
  G4double pitch;
  G4int NumOfSegPD = 0;
  auto sdPD = new PDSD(name);
  AddNewDetector(sdPD);
  // Mother
  auto mother_solid = new G4Box("PdMotherSolid",
				half_size.x()*mm,
                                half_size.y()*mm,
                                half_size.z()*mm);
  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       mat("Air"),
                                       "PdMotherLV");
  auto rot = new G4RotationMatrix;
  rot->rotateY(std::abs(ra2));
  auto pos = gGeom.GetGlobalPosition(RCname);
  G4double GlobalX_pm;
  if(i == 1) GlobalX_pm =  1.;
  else       GlobalX_pm = -1.;
  if(YZ == "Y"){
    NumOfSegPD = NumOfSegPDY;
    pitch = seg_size.y()*2*mm;
    pos.setX(pos.x()
	     + GlobalX_pm*(NumOfSegRC/2*Rcseg_size.z()*2 + seg_size.z()*3));
  }
  else{
    NumOfSegPD = NumOfSegPDZ;
    pitch = seg_size.x()*2*mm;
    pos.setX(pos.x()
	     + GlobalX_pm*(NumOfSegRC/2*Rcseg_size.z()*2 + seg_size.z()));
  }
  new G4PVPlacement(rot, pos, mother_lv, "PdMotherPHYS", m_world_lv, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // Segment
  auto segment_solid = new G4Box("PdSegmentSolid", seg_size.x(),
                                 seg_size.y(), seg_size.z());
  auto segment_lv = new G4LogicalVolume(segment_solid,
                                        mat("Scintillator"),
                                        "PdSegmentLV");
  for(G4int i=0; i<NumOfSegPD; ++i){
    segment_lv->SetSensitiveDetector(sdPD);
    //pos = G4ThreeVector(0.0, 0.0, 0.0); // in case of plane
    if(YZ == "Y"){
      segment_lv->SetVisAttributes(G4Colour::Green());
      pos = G4ThreeVector(0.0, (-NumOfSegPD/2 + i)*pitch + pitch/2, 0.0);
    }
    else{
      segment_lv->SetVisAttributes(G4Colour::Magenta());
      pos = G4ThreeVector((-NumOfSegPD/2 + i)*pitch + pitch/2, 0.0, 0.0);
    }
    new G4PVPlacement(nullptr, pos, segment_lv,
                      "PdSegmentPHYS", mother_lv, false, i, m_check_overlaps);
  }
}
