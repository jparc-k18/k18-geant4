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

#include "FuncName.hh"
#include "HTOFSD.hh"
#include "TPCPadHelper.hh"
#include "TPCSD.hh"

namespace
{
namespace tpc = TPCPadHelper;
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
void S2SDetectorConstruction::ConstructHTOF()
{
  using CLHEP::mm;
  using CLHEP::deg;
  auto htof_sd = new HTOFSD("HTOF");
  AddNewDetector(htof_sd);
  const auto& htof_pos = gGeom.GetGlobalPosition("HTOF");
  const auto& half_size = gSize.GetSize("HtofSeg") * 0.5 * mm;
  const G4double L = gGeom.GetLocalZ("HTOF");
  const G4double dXdW = gGeom.GetWirePitch("HTOF");
  // Segment
  // Sintillator
  auto htof_scintilltor = new G4Box("HtofScint", half_size.x(),
                                    half_size.y(), half_size.z());
  // Light-guides
  //const G4ThreeVector lg_size(35.0,100.0,4.0);
  const G4ThreeVector lg_size(17.5,50.0,2.0);
  // Upper one
  std::vector<G4TwoVector> upper_lg_vertices;
  upper_lg_vertices.push_back(G4TwoVector(-half_size.x() , -half_size.z()));
  upper_lg_vertices.push_back(G4TwoVector(-half_size.x() , half_size.z()));
  upper_lg_vertices.push_back(G4TwoVector(half_size.x() , half_size.z()));
  upper_lg_vertices.push_back(G4TwoVector(half_size.x() , -half_size.z()));
  upper_lg_vertices.push_back(G4TwoVector(-lg_size.x() , -half_size.z()));
  upper_lg_vertices.push_back(G4TwoVector(-lg_size.x() , -half_size.z() + lg_size.z()));
  upper_lg_vertices.push_back(G4TwoVector(lg_size.x() , -half_size.z() + lg_size.z()));
  upper_lg_vertices.push_back(G4TwoVector(lg_size.x() , -half_size.z()));
  auto htof_upper_lg = new G4GenericTrap("HtofLG_upper", lg_size.y(), upper_lg_vertices);
  auto rotM_upper_lg = new G4RotationMatrix;
  rotM_upper_lg->rotateX(90.0 *deg);
  rotM_upper_lg->rotateZ(- 180.0 *deg);
  G4ThreeVector trans_upper_lg(0.*mm , half_size.y() + lg_size.y(),  0.*mm);

  // Lower one
  std::vector<G4TwoVector> lower_lg_vertices;
  lower_lg_vertices.push_back(G4TwoVector(-half_size.x() , -half_size.z()));
  lower_lg_vertices.push_back(G4TwoVector(-half_size.x() , half_size.z()));
  lower_lg_vertices.push_back(G4TwoVector(half_size.x() , half_size.z()));
  lower_lg_vertices.push_back(G4TwoVector(half_size.x() , -half_size.z()));
  lower_lg_vertices.push_back(G4TwoVector(-lg_size.x() , half_size.z() - lg_size.z()));
  lower_lg_vertices.push_back(G4TwoVector(-lg_size.x() , half_size.z()));
  lower_lg_vertices.push_back(G4TwoVector(lg_size.x() , half_size.z()));
  lower_lg_vertices.push_back(G4TwoVector(lg_size.x() , half_size.z() - lg_size.z()));
  auto htof_lower_lg = new G4GenericTrap("HtofLG_lower", lg_size.y(), lower_lg_vertices);
  auto rotM_lower_lg = new G4RotationMatrix;
  rotM_lower_lg->rotateX(- 90.0 *deg);
  rotM_lower_lg->rotateZ(- 180.0 *deg);
  G4ThreeVector trans_lower_lg(0.*mm , -half_size.y() - lg_size.y() , 0.*mm);
  auto solid_dummy = new G4UnionSolid("dummysolid", htof_scintilltor, htof_upper_lg,
                                      rotM_upper_lg, trans_upper_lg);

  auto htof_solid = new G4UnionSolid("HtofSolid", solid_dummy,
                                     htof_lower_lg, rotM_lower_lg,
                                     trans_lower_lg);
  //Common slats
  auto htof_lv = new G4LogicalVolume(htof_solid, mat("Scintillator"),
                                     "HtofLV");

  //HTOF beam-through part
  G4double HTOF_window=110.*mm;
  auto window_dummy = new G4Box("window_dummy", half_size.x(), half_size.y()/2. - HTOF_window/4., half_size.z());

  //Upper slats(Beam-through)
  G4ThreeVector trans_upper_window(0.*mm , half_size.y()/2. - HTOF_window/4. + lg_size.y(),  0.*mm);
  auto htof_solid_upper = new G4UnionSolid("HtofSolid_upper", window_dummy, htof_upper_lg,
                                           rotM_upper_lg, trans_upper_window);
  auto htof_upper_lv = new G4LogicalVolume(htof_solid_upper, mat("Scintillator"),
                                           "HtofUpperLV");
  //Lower slats(Beam-through)
  G4ThreeVector trans_lower_window(0.*mm , -half_size.y()/2. + HTOF_window/4. - lg_size.y(),  0.*mm);
  auto htof_solid_lower = new G4UnionSolid("HtofSolid_lower", window_dummy, htof_lower_lg,
                                           rotM_lower_lg, trans_lower_window);
  auto htof_lower_lv = new G4LogicalVolume(htof_solid_lower, mat("Scintillator"),
                                           "HtofLowerLV");

  for(G4int i=0; i<NumOfPlaneHTOF; ++i){
    for(G4int j=0; j<NumOfSegHTOFOnePlane; ++j){
      G4int seg = i*NumOfSegHTOFOnePlane + j;
      // (lateral, height, radial)
      G4ThreeVector seg_pos(-dXdW * (j - (NumOfSegHTOFOnePlane - 1)/2.),
                            0.*mm,
                            -L);
      auto rotMOutP = new G4RotationMatrix;
      rotMOutP->rotateY(- i * 360./NumOfPlaneHTOF*deg);
      seg_pos.rotateY(i * 360./NumOfPlaneHTOF*deg);
      seg_pos += htof_pos;
      G4int copy_no = seg+2;

      G4ThreeVector window_pos(0.*mm, half_size.y()/2. + HTOF_window/4., 0.*mm);
      //common slats
      if(i!=0)	new G4PVPlacement(rotMOutP, seg_pos, htof_lv, Form("HtofPV%d", copy_no), m_world_lv, false, copy_no, m_check_overlaps);
      else if(j==0) new G4PVPlacement(rotMOutP, seg_pos, htof_lv, Form("HtofPV%d", 0), m_world_lv, false, 0, m_check_overlaps);
      else if(j==3) new G4PVPlacement(rotMOutP, seg_pos, htof_lv, Form("HtofPV%d", 5), m_world_lv, false, 5, m_check_overlaps);
      //Beam-through slats
      else if(j==1){
	new G4PVPlacement(rotMOutP, seg_pos + window_pos, htof_upper_lv, Form("HtofPV%d", 1), m_world_lv, false, 1, m_check_overlaps);
	new G4PVPlacement(rotMOutP, seg_pos - window_pos, htof_lower_lv, Form("HtofPV%d", 2), m_world_lv, false, 2, m_check_overlaps);
      }
      else if(j==2){
	new G4PVPlacement(rotMOutP, seg_pos + window_pos, htof_upper_lv, Form("HtofPV%d", seg), m_world_lv, false, 3, m_check_overlaps);
	new G4PVPlacement(rotMOutP, seg_pos - window_pos, htof_lower_lv, Form("HtofPV%d", seg+31), m_world_lv, false, 4, m_check_overlaps);
      }
    }
  }

  htof_lv->SetSensitiveDetector(htof_sd);
  htof_upper_lv->SetSensitiveDetector(htof_sd);
  htof_lower_lv->SetSensitiveDetector(htof_sd);
  htof_lv->SetVisAttributes(G4Colour::Cyan());
  htof_upper_lv->SetVisAttributes(G4Colour::Cyan());
  htof_lower_lv->SetVisAttributes(G4Colour::Cyan());

#if 0
  // Supporting frame parts
  // Dummy for subtraction
  auto RingHoleOut = new G4Box("RingHoleOut", 127.5*mm, 6.0*mm, 5.*mm);
  auto RingHoleIn = new G4Box("RingHoleIn", 10.*mm, 6.*mm, 5.*mm);
  auto RingHole = new G4SubtractionSolid("RingHole" , RingHoleOut, RingHoleIn);

  // Top Ring
  const G4double RingIn_zPlane[2] = { -5.*mm, 5.*mm };
  const G4double RingIn_rInner[2] = { 0.*mm, 0.*mm };
  const G4double RingIn_rOuter[2] =  { 334.13*mm, 334.13*mm } ;
  auto TopRingOutSolid = new G4Tubs("TopRingOutSolid", 0.*mm, 391.97*mm,
                                    5.*mm, 0.*deg, 360.*deg);
  auto TopRingInSolid = new G4Polyhedra("TopRingInSolid", 22.5*deg, (360. + 22.5)*deg,
                                        8, 2,
                                        RingIn_zPlane, RingIn_rInner, RingIn_rOuter);

  G4SubtractionSolid* TopRingSubSolid[8];
  for(G4int i=0; i<8; i++){
    auto rotMHole = new G4RotationMatrix;
    rotMHole->rotateZ(i * 45.0 * deg);
    G4ThreeVector trans_hole(0*mm, 348.13*mm, 0.*mm);
    trans_hole.rotateZ(- i * 45.0 * deg);
    if(i==0) TopRingSubSolid[0] = new G4SubtractionSolid("TopRingSubSolid_0", TopRingOutSolid, RingHole, rotMHole, trans_hole);
    else TopRingSubSolid[i] = new G4SubtractionSolid(Form("TopRingSubSolid_%d", i), TopRingSubSolid[i-1], RingHole, rotMHole, trans_hole);
  }

  auto TopRingSolid = new G4SubtractionSolid("TopRingSoild" , TopRingSubSolid[7], TopRingInSolid);
  auto TopRing_lv = new G4LogicalVolume(TopRingSolid, mat("Iron"), "TopRingLV");
  TopRing_lv->SetVisAttributes(ORANGE);

  // Bottom Ring
  const G4double RingOut_zPlane[2] = { -5.*mm, 5.*mm };
  const G4double RingOut_rInner[2] = { 334.13*mm, 334.13*mm };
  const G4double RingOut_rOuter[2] = { 362.13*mm, 362.13*mm };
  auto BotRingOutSolid = new G4Polyhedra("BotRingOutSolid", 22.5*deg, (360. + 22.5)*deg,
                                         8, 2,
                                         RingOut_zPlane, RingOut_rInner, RingOut_rOuter);

  G4SubtractionSolid* BotRingSubSolid[8];
  for(G4int i=0; i<8; i++){
    auto rotMHole = new G4RotationMatrix;
    rotMHole->rotateZ(i * 45.0 * deg);
    G4ThreeVector trans_hole(0.*mm, 348.13*mm, 0.*mm);
    trans_hole.rotateZ(- i * 45.0 * deg);
    if(i==0) BotRingSubSolid[i] = new G4SubtractionSolid(Form("BotRingSubSolid_%d", i) , BotRingOutSolid, RingHole, rotMHole, trans_hole);
    else BotRingSubSolid[i] = new G4SubtractionSolid(Form("BotRingSubSolid_%d", i) , BotRingSubSolid[i-1], RingHole, rotMHole, trans_hole);
  }
  auto BotRing_lv = new G4LogicalVolume(BotRingSubSolid[7], mat("Iron"),"BotRingLV");
  BotRing_lv->SetVisAttributes(G4Colour::Blue());

  // Bracket
  auto BraInSolid = new G4Box("BraInSolid", 137.5*mm, 75.*mm, 2.5*mm);
  auto BraInHole = new G4Box("BraInHole", 127.5*mm, 30.*mm, 2.5*mm);
  G4ThreeVector brain_trans_hole(0.*mm, 25.*mm, 0.*mm);
  auto BraInSubSolid = new G4SubtractionSolid("BraInSubSolid", BraInSolid, BraInHole,
                                              0, brain_trans_hole);

  auto BraOutSolid = new G4Box("BraOutSolid", 137.5*mm, 2.5*mm, 16.5*mm);
  G4ThreeVector braout_trans_hole(0.*mm, 0.*mm, 3.5*mm);
  auto BraOutSubSolid = new G4SubtractionSolid("BraOutSubSolid", BraOutSolid, RingHole,
                                               0, braout_trans_hole);

  G4ThreeVector trans_braket(0.*mm, 72.5*mm, 14.*mm);
  auto BraSolid = new G4UnionSolid("BraSolid", BraInSubSolid, BraOutSubSolid,
                                   0, trans_braket);
  auto Bra_lv = new G4LogicalVolume(BraSolid, mat("Iron"),
                                    "BraLV");
  Bra_lv->SetVisAttributes(G4Colour::Green());

  // Preamp Support Frame
  const G4double PreFrameIn_zPlane[2] = { -37.*mm, 37.*mm };
  const G4double PreFrameIn_rInner[2] = { 334.13*mm, 334.13*mm };
  const G4double PreFrameIn_rOuter[2] = { 339.13*mm, 339.13*mm } ;
  auto PreFrameInSolid = new G4Polyhedra("PreFrameInSolid", 22.5*deg, (360. + 22.5)*deg,
                                         8, 2,
                                         PreFrameIn_zPlane,  PreFrameIn_rInner,  PreFrameIn_rOuter);

  const G4double PreFrameOut_zPlane[2] = { -32.5*mm, 32.5*mm };
  const G4double PreFrameOut_rInner[2] = { 354.13*mm, 354.13*mm };
  const G4double PreFrameOut_rOuter[2] = { 362.13*mm, 362.13*mm } ;
  auto PreFrameOutSolid = new G4Polyhedra("PreFrameOutSolid", 22.5*deg, (360. + 22.5)*deg,
                                          8, 2,
                                          PreFrameOut_zPlane,  PreFrameOut_rInner,  PreFrameOut_rOuter);
  auto PreFrameHole = new G4Box("PreFrameHole", 3.5*mm, 4.0*mm, 32.5*mm);

  G4SubtractionSolid* PreFrameOutSubSolid[8];
  for(G4int i=0; i<8; i++){
    auto rotMHole = new G4RotationMatrix;
    rotMHole->rotateZ(i * 45.0 * deg);
    G4ThreeVector trans_hole(0.*mm, 348.13*mm, 0.*mm);
    trans_hole.rotateZ(- i * 45.0 * deg);
    if(i==0) PreFrameOutSubSolid[0] = new G4SubtractionSolid("PreFrameOutSubSolid_0",
                                                             PreFrameOutSolid, PreFrameHole);
    else PreFrameOutSubSolid[i] = new G4SubtractionSolid(Form("PreFrameOutSubSolid_%d", i),
                                                         PreFrameOutSubSolid[i-1], PreFrameHole,
                                                         rotMHole, trans_hole);
  }
  G4ThreeVector trans_preframe(0.*mm, 0.*mm, -9.0*mm);
  auto PreFrameSolid = new G4UnionSolid("PreFrameSolid", PreFrameInSolid, PreFrameOutSubSolid[7],
                                        0, trans_preframe);
  auto PreFrame_lv = new G4LogicalVolume(PreFrameSolid, mat("Aluminum"),
					 "PreFrameLV");
  PreFrame_lv->SetVisAttributes(G4Colour::Red());

  // Bar
  auto BarMainSolid = new G4Box("BarMainSolid", 45.*mm, 618.25*mm, 5.*mm);
  auto BarSideSolid = new G4Box("BarSideSolid", 5.*mm, 618.25*mm, 5.*mm);
  G4ThreeVector trans_bar_side(-40.*mm, 0.*mm, -10.*mm);
  auto BarUniSolid = new G4UnionSolid("BarUniSolid", BarMainSolid, BarSideSolid,
                                      0, trans_bar_side);

  auto BarBotSolid = new G4Box("BarSideSolid", 45.*mm, 5.*mm, 5.*mm);
  G4ThreeVector trans_bar_bot(0.*mm, -613.25*mm, -10.*mm);

  auto BarSolid = new G4UnionSolid("BarSolid", BarUniSolid, BarBotSolid,
                                   0, trans_bar_bot);
  auto Bar_lv = new G4LogicalVolume(BarSolid, mat("Iron"),
                                    "BarLV");
  Bar_lv->SetVisAttributes(G4Colour::Blue());

  //if(true){   // Htof Frame Placement
  if(false){   // Htof Frame Placement

    auto rotMOutRing = new G4RotationMatrix;
    rotMOutRing->rotateX(- 90. *deg);

    // Top Ring
    G4ThreeVector TopRing_pos(0.*mm, 586.72*mm, 0.*mm);
    TopRing_pos += htof_pos;
    new G4PVPlacement(rotMOutRing, TopRing_pos, TopRing_lv, "TopRingPV", m_world_lv, false, 0, m_check_overlaps);

    // Bottom Ring
    G4ThreeVector BotRing_pos(0.*mm, -586.72*mm, 0.*mm);
    BotRing_pos += htof_pos;
    new G4PVPlacement(rotMOutRing, BotRing_pos, BotRing_lv, "BotRingPV",
                      m_world_lv, false, 0, m_check_overlaps);

    // Preamp Support Frame
    for(G4int i=0; i<2; ++i){
      auto rotMOutP_PreFrame = new G4RotationMatrix;
      rotMOutP_PreFrame->rotateX((1 - 2 * i) * 90.*deg);
      rotMOutP_PreFrame->rotateZ(i * 180.*deg);

      G4ThreeVector PreFrame_pos(0.*mm, 453.22*mm, 0.*mm);
      PreFrame_pos.rotateZ(- i * 180.*deg);
      PreFrame_pos += htof_pos;
      new G4PVPlacement(rotMOutP_PreFrame, PreFrame_pos, PreFrame_lv,
                        Form("PreFramePV%d", i),
                        m_world_lv, false, i, m_check_overlaps);
    }

    G4int Bar_seg=0;
    for(G4int i=0; i<NumOfPlaneHTOF; ++i){
      //Bar
      auto rotMOutP_bar = new G4RotationMatrix;
      rotMOutP_bar->rotateY(- i * 45.*deg);
      G4ThreeVector Bar_pos(0.*mm, - 36.53*mm, - 372.13*mm);
      Bar_pos.rotateY(i * 45.*deg);
      Bar_pos += htof_pos;
      if(i!=0 && i!=4){ //Beam through
	new G4PVPlacement(rotMOutP_bar, Bar_pos, Bar_lv,
                          Form("BarPV%d", Bar_seg),
                          m_world_lv, false, Bar_seg, m_check_overlaps);
	Bar_seg++;
      }
      //Bracket
      for(G4int j=0; j<2; ++j){
	auto rotMOutP_bra = new G4RotationMatrix;
	rotMOutP_bra->rotateY( - (1 - 2 * j) * i * 45.*deg);
	rotMOutP_bra->rotateZ( - j * 180.*deg);
	G4ThreeVector Bra_pos(0.*mm, 506.72*mm, - 364.63*mm);
	Bra_pos.rotateY(i * 45.*deg);
	Bra_pos.rotateZ(j * 180.*deg);
	Bra_pos += htof_pos;
	G4int Bra_seg = 2 * i + j;
	new G4PVPlacement(rotMOutP_bra, Bra_pos, Bra_lv,
                          Form("BraPV%d", Bra_seg),
                          m_world_lv, false, Bra_seg, m_check_overlaps);
      }
    }
  }
#endif
}

//_____________________________________________________________________________
void S2SDetectorConstruction::ConstructHypTPC()
{
  using CLHEP::mm;
  using CLHEP::deg;
  auto tpc_sd = new TPCSD("TPC");
  AddNewDetector(tpc_sd);
  const auto tpc_pos = gGeom.GetGlobalPosition("HypTPC")*mm;
  const auto target_pos = tpc_pos + G4ThreeVector(0, 0, -143*mm);
  // gGeom.GetGlobalPosition("Target")*mm;
  const auto& gas_vessel_size = gSize.GetSize("TpcGasVessel")*mm*0.5;
  const auto& gas_vessel_window = gSize.GetSize("TpcGasVesselWindow")*mm*0.5;
  const auto& target_holder_size = gSize.GetSize("TargetHolder")*mm*0.5;
  const auto& target_holder_yofs = gas_vessel_size[2]-target_holder_size[2];
  const auto& target_holder_window = gSize.GetSize("TargetHolderWindowPhiDz")*mm;
  const auto target_holder_thickness =
    (target_holder_size[1] - target_holder_size[0])/2./2.; // 3*mm/2.
  const auto& p10_size = gSize.GetSize("TpcP10Volume")*mm/2.;
  const auto& field_cage_size = gSize.GetSize("TpcFieldCage")*mm/2.;
  const auto& eff_volume_size = gSize.GetSize("TpcEffectiveVolume")*mm/2.;
  const G4double phiStart = 22.5*deg;
  const G4double phiTotal = 360*deg;
  const G4int numSide   = 8;
  const G4int numZPlane = 2;
  // Gas Vessel
  const G4double zPlaneGV[] = { -gas_vessel_size.z(),
                                gas_vessel_size.z() };
  const G4double rInnerGV[] = { 0., 0. };
  const G4double rOuterGV[] = { gas_vessel_size[1],
                                gas_vessel_size[1] };
  G4VSolid* gv_solid;
  gv_solid = new G4Polyhedra("TpcGasVesselSolid",
                             phiStart, phiTotal, numSide, numZPlane,
                             zPlaneGV, rInnerGV, rOuterGV);
  auto space_solid = new G4Tubs("TpcTargetSpaceSolid",
                                0.,
                                target_holder_size[0],
                                target_holder_size[2],
                                phiStart, phiTotal);
  G4ThreeVector pos = target_pos;
  pos.setY(target_holder_yofs + target_holder_thickness*2.);
  pos.rotateX(90.*deg);
  gv_solid = new G4SubtractionSolid("TpcGasVesselSolid",
                                    gv_solid, space_solid,
                                    nullptr, pos);
  auto rot = new G4RotationMatrix;
  rot->rotateX(90.*deg);
  rot->rotateZ(180.*deg); // Target is downstream
  auto gv_lv = new G4LogicalVolume(gv_solid,
                                   mat("Al"),
                                   "TpcGasVesselLV");
  new G4PVPlacement(rot, tpc_pos, gv_lv, "TpcGasVesselPV",
                    m_world_lv, false, 0, m_check_overlaps);
  gv_lv->SetVisAttributes(G4Colour::White());
  // Target holder
  auto th_side_solid = new G4Tubs("TargetHolderSideSolid",
                                  target_holder_size[0],
                                  target_holder_size[1],
                                  target_holder_size[2],
                                  phiStart, phiTotal);
  auto th_bottom_solid = new G4Tubs("TargetHolderBottomSolid",
                                    0.,
                                    target_holder_size[0],
                                    target_holder_thickness,
                                    phiStart, phiTotal);
  pos.set(0, 0, -target_holder_size[2] + target_holder_thickness);
  auto th_solid = new G4UnionSolid("TargetholderSolid",
                                   th_side_solid, th_bottom_solid,
                                   nullptr, pos);
  auto th_lv = new G4LogicalVolume(th_solid,
                                   mat("G10"),
                                   "TargetHolderLV");
  pos = (target_pos - tpc_pos);
  pos.setY(target_holder_yofs);
  pos.rotateX(90.*deg);
  new G4PVPlacement(nullptr, pos, th_lv, "TargetHolderPV",
                    gv_lv, false, 0, m_check_overlaps);
  th_lv->SetVisAttributes(G4Colour::Green());
  // Target holder window
  auto th_window_solid = new G4Tubs("TargetHolderWindowSolid",
                                    target_holder_size[0],
                                    target_holder_size[1],
                                    target_holder_window[2]/2.,
                                    target_holder_window[0]*deg,
                                    target_holder_window[1]*deg);
  auto th_window_lv = new G4LogicalVolume(th_window_solid,
                                          mat("P10"),
                                          "TargetHolderWindowLV");
  rot = new G4RotationMatrix;
  rot->rotateZ(180.*deg);
  pos.set(0, 0, -target_holder_yofs);
  new G4PVPlacement(nullptr, pos, th_window_lv,
                    "TargetHolderWindowPV",
                    th_lv, false, 0, m_check_overlaps);
  new G4PVPlacement(rot, pos, th_window_lv,
                    "TargetHolderWindowPV",
                    th_lv, false, 1, m_check_overlaps);
  th_window_lv->SetVisAttributes(G4Colour::White());

  auto gv_window_solid = new G4Box("GasVesselWindowSolid",
                                   gas_vessel_window.x(),
                                   gas_vessel_window.y(),
                                   // gas_vessel_window.z()
                                   (gas_vessel_size[1] - p10_size[1])/2.
                                   );
  auto gv_window_lv = new G4LogicalVolume(gv_window_solid,
                                          mat("P10"),
                                          "GasVesselWindowLV");
  gv_window_lv->SetVisAttributes(G4Colour::White());
  pos.set(0., (gas_vessel_size[1]+p10_size[1])/2., 0.);
  for (G4int i=0; i<numSide; ++i) {
    rot = new G4RotationMatrix;
    rot->rotateX(90.*deg);
    rot->rotateY((i+1)*360.*deg/numSide);
    pos.rotateZ(360.*deg/numSide);
    new G4PVPlacement(rot, pos, gv_window_lv,
                      "GasVesselWindowPV"+std::to_string(i),
                      gv_lv, false, i, m_check_overlaps);
  }

  // P10
  const G4double rInnerP10[] = { p10_size[0], p10_size[0] };
  const G4double rOuterP10[] = { p10_size[1], p10_size[1] };
  const G4double zPlaneP10[] = { -p10_size[2], p10_size[2] };
  G4VSolid* p10_solid;
  p10_solid = new G4Polyhedra("TpcP10Solid",
                              phiStart, phiTotal, numSide, numZPlane,
                              zPlaneP10, rInnerP10, rOuterP10);
  pos = (target_pos - tpc_pos);
  pos.setY(target_holder_yofs);
  pos.rotateX(90.*deg);
  p10_solid = new G4SubtractionSolid("TpcP10Solid",
                                     p10_solid, space_solid,
                                     nullptr, pos);
  p10_solid = new G4SubtractionSolid("TpcP10Solid",
                                     p10_solid, th_side_solid,
                                     nullptr, pos);
  auto p10_lv = new G4LogicalVolume(p10_solid,
                                    mat("P10"),
                                    "TpcP10LV");
  new G4PVPlacement(nullptr, G4ThreeVector(), p10_lv, "TpcP10PV",
                    gv_lv, false, 0, m_check_overlaps);
  p10_lv->SetVisAttributes(G4Colour::Yellow());
  // p10_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  p10_lv->SetSensitiveDetector(tpc_sd);
  // Field Cage
  const G4double rInnerFC[] = { field_cage_size[0], field_cage_size[0] };
  const G4double rOuterFC[] = { field_cage_size[1], field_cage_size[1] };
  const G4double zPlaneFC[] = { -field_cage_size[2], field_cage_size[2] };
  auto fc_solid = new G4Polyhedra("FieldCageSolid",
                                  phiStart, phiTotal, numSide, numZPlane,
                                  zPlaneFC, rInnerFC, rOuterFC);
  auto fc_lv = new G4LogicalVolume(fc_solid, mat("G10"),
                                   "FieldCageLV");
  new G4PVPlacement(nullptr, G4ThreeVector(), fc_lv,
                    "FieldCagePV", p10_lv, false, 0, m_check_overlaps);
  fc_lv->SetVisAttributes(G4Colour::Green());
  // Virtual pads
  G4LogicalVolume* pad_lv[NumOfPadTPC];
  [[maybe_unused]] G4double angle[NumOfPadTPC] = {};
  const G4double pad_center_z = gSize.Get("TpcPadCenterZ")*mm;
  G4double pad_in[NumOfPadTPC] = {};
  G4double pad_out[NumOfPadTPC] = {};
  G4double tpc_rad = 250;
  // out side less 100 mm. 10+5*x < 100 mm is pad_in_num
  const G4double pad_length_in = gSize.Get("TpcPadLengthIn");
  const G4double pad_length_out = gSize.Get("TpcPadLengthOut");
  const G4double pad_gap = gSize.Get("TpcPadGap");
  const G4int pad_configure = gSize.Get("TpcPadConfigure");
  switch (pad_configure) {
  case 1:
    for (G4int i=0; i<NumOfPadTPC; ++i) {
      if (i<NumOfPadTPCIn) {
	pad_in[i]  = 10.+(pad_length_in+pad_gap)*i;
	pad_out[i] = 10.+(pad_length_in+pad_gap)*i+pad_length_in;
	angle[i]   = 360.;
      } else {
	pad_in[i] = 10.+(pad_length_in+pad_gap)*NumOfPadTPCIn +
	  (pad_length_out+pad_gap)*(i-NumOfPadTPCIn);
	pad_out[i] = 10.+(pad_length_in+pad_gap)*NumOfPadTPCIn +
	  (pad_length_out+pad_gap)*(i-NumOfPadTPCIn) + pad_length_out;
	angle[i] =
	  -std::acos((square(pad_out[i]) +
                      square(pad_center_z) -
                      square(tpc_rad)) /
                     (2*pad_out[i]*pad_center_z))*RadToDeg + 180.;
      }
    }
    break;
  case 2:
    for (G4int i=0; i<NumOfPadTPC; ++i) {
      if (i<NumOfPadTPCIn) {
	pad_in[i]  = 10.+(pad_length_in+pad_gap)*i;
	pad_out[i] = 10.+(pad_length_in+pad_gap)*i+pad_length_in;
	angle[i]   = 360.;
      } else {
	pad_in[i] = 10.+(pad_length_in+pad_gap)*NumOfPadTPCIn +
	  (pad_length_out+pad_gap)*(i-NumOfPadTPCIn);
	pad_out[i] = 10.+(pad_length_in+pad_gap)*NumOfPadTPCIn +
	  (pad_length_out+pad_gap)*(i-NumOfPadTPCIn) + pad_length_out;
      }
    }
    angle[10] = 180. - 155.35;
    angle[11] = 180. - 144.8;
    angle[12] = 180. - 138.;
    angle[13] = 180. - 116.73;
    angle[14] = 180. - 106.;
    angle[15] = 180. - 98.77;
    angle[16] = 180. - 94.29;
    angle[17] = 180. - 89.8;
    angle[18] = 180. - 87.18;
    angle[19] = 180. - 84.16;
    angle[20] = 180. - 81.48;
    angle[21] = 180. - 73.39;
    angle[22] = 180. - 65.51011;
    angle[23] = 180. - 60.19;
    angle[24] = 180. - 56.35239;
    angle[25] = 180. - 52.85;
    angle[26] = 180. - 50.14;
    angle[27] = 180. - 47.17;
    angle[28] = 180. - 41.24;
    angle[29] = 180. - 29.;
    angle[30] = 180. - 23.23;
    angle[31] = 180. - 18.69;
    break;
  case 3:
    //for tracking analysis
    //If you need the dE/dx information, it should be modified.
    //Thin sensitive detector is introduced.
    for (G4int i=0; i<NumOfPadTPC; ++i) {
      G4double pad_radius = tpc::getRadius(i);
      pad_in[i] = pad_radius;
      pad_out[i] = pad_radius + 0.1*mm;
      /*
        double pad_halflength = tpc::getLength(i)/2;
        pad_in[i] = pad_radius-pad_halflength;
        pad_out[i] = pad_radius + pad_halflength;
      */
      if (i<NumOfPadTPCIn) {
	angle[i]   = 360.;
      } else {
	angle[i] = tpc::getsTheta(i);
      }
    }
    break;
  default:
    break;
  }
  // Inner Pads
  const G4double rInnerEA[] = { eff_volume_size[0], eff_volume_size[0] };
  const G4double rOuterEA[] = { eff_volume_size[1], eff_volume_size[1] };
  const G4double zPlaneEA[] = { -eff_volume_size[2], eff_volume_size[2] };
  auto eff_volume = new G4Polyhedra("TpcEffectiveVolumeSolid",
                                  phiStart, phiTotal, numSide, numZPlane,
                                  zPlaneEA, rInnerEA, rOuterEA);
  G4VSolid* pad_solid[NumOfPadTPC];
  pos.set(0, pad_center_z, 0);
  for (G4int i=0; i<NumOfPadTPCIn; ++i) {
    pad_solid[i] = new G4Tubs("TpcPadSolid"+std::to_string(i),
                              pad_in[i]*mm,
                              pad_out[i]*mm,
                              field_cage_size[2],
                              phiStart, phiTotal);
    pad_solid[i] = new G4IntersectionSolid("TpcPadSolid"+std::to_string(i),
                                           pad_solid[i], p10_solid,
                                           nullptr, pos);
    pad_solid[i] = new G4IntersectionSolid("TpcPadSolid"+std::to_string(i),
                                           pad_solid[i], eff_volume,
                                           nullptr, pos);
    pad_lv[i]  = new G4LogicalVolume(pad_solid[i], mat("P10"),
                                     "TpcPadLV"+std::to_string(i));
    new G4PVPlacement(nullptr, -pos, pad_lv[i], "TpcPadPV"+std::to_string(i),
                      p10_lv, true, i, m_check_overlaps);
  }
  // Outer Pads
  for(G4int i=NumOfPadTPCIn; i<NumOfPadTPC; ++i){
    pad_solid[i] = new G4Tubs("TpcPadSolid"+std::to_string(i),
                              pad_in[i]*mm,
                              pad_out[i]*mm,
                              field_cage_size[2],
                              phiStart, phiTotal);
    pad_solid[i] = new G4IntersectionSolid("TpcPadSolid"+std::to_string(i),
                                           pad_solid[i], p10_solid,
                                           nullptr, pos);
    pad_solid[i] = new G4IntersectionSolid("TpcPadSolid"+std::to_string(i),
                                           pad_solid[i], eff_volume,
                                           nullptr, pos);
    pad_lv[i]  = new G4LogicalVolume(pad_solid[i], mat("P10"),
				     "TpcPadLV"+std::to_string(i));
    new G4PVPlacement(nullptr, -pos, pad_lv[i], "TpcPadPV"+std::to_string(i),
                      p10_lv, true, i, m_check_overlaps);
  }
  for (G4int i=0; i<NumOfPadTPC; ++i) {
    pad_lv[i]->SetVisAttributes(ORANGE);
    pad_lv[i]->SetSensitiveDetector(tpc_sd);
  }
  // Dead area
  auto dead_solid = new G4Box("DeadSolid", 5*mm, 250*mm, 0.001*mm);
  auto dead_lv = new G4LogicalVolume(dead_solid, mat("P10"),
                                     "DeadLV");
  auto rotdead1 = new G4RotationMatrix;
  rotdead1->rotateZ(45.*deg);
  new G4PVPlacement(rotdead1, G4ThreeVector(0., 0.*mm, -300.1*mm),
                    dead_lv, "DeadPV1", p10_lv, true, 0, m_check_overlaps);
  auto rotdead2 = new G4RotationMatrix;
  rotdead2->rotateZ(-45.*deg);
  new G4PVPlacement(rotdead2, G4ThreeVector(0., 0.*mm, -300.1*mm),
                    dead_lv, "DeadPV2", p10_lv, true, 1, m_check_overlaps);
  dead_lv->SetVisAttributes(G4Colour::Gray());
}

//_____________________________________________________________________________
void S2SDetectorConstruction::ConstructTargetE90()
{
  using CLHEP::mm;
  using CLHEP::deg;
  // auto target_sd = new TargetSD("TGT");
  // AddNewDetector(target_sd);
  const auto target_pos = gGeom.GetGlobalPosition("Target")*mm;
  const auto target_size = gSize.GetSize("Target")*mm/2.;
  G4VSolid* target_solid;
  auto rot = new G4RotationMatrix;
  switch (m_experiment) {
  case 42: {
    target_solid = new G4Box("Target", target_size.x(),
                             target_size.y(), target_size.z());
  }
    break;
  case 45: case 27: {
    G4double target_r = gSize.Get("Target", G4ThreeVector::X);
    G4double target_z = gSize.Get("Target", G4ThreeVector::Z);
    target_solid = new G4Tubs("TargetSolid", 0.*mm,
                              target_r*mm, target_z*mm, 0., 360*deg);
    rot->rotateX(90.*deg);
    break;
  }
  case 90: {
    if (gConf.Get<G4bool>("TargetVP")) {
      target_solid = new G4Box("TargetSolid",
                               target_size[1],
                               0.001*mm,
                               target_size[2]);
    } else {
      target_solid = new G4Tubs("TargetSolid",
                                target_size[0],
                                target_size[1],
                                target_size[2],
                                0.*deg, 360.*deg);
    }
    rot->rotateX(90.*deg);
    const auto kapton_size = gSize.GetSize("TargetKapton")*mm/2.;
    const auto gfrp_size = gSize.GetSize("TargetGFRP")*mm/2.;
    auto kapton = new G4Tubs("TargetKapton",
                             kapton_size[0],
                             kapton_size[1],
                             kapton_size[2],
                             0.*deg, 360.*deg);
    auto kapton_lv = new G4LogicalVolume(kapton,
                                         mat("Kapton"),
                                         "TargetKaptonLV");
    kapton_lv->SetVisAttributes(G4Colour::Red());
    new G4PVPlacement(rot, target_pos, kapton_lv, "TargetKaptonPV",
                      m_world_lv, true, 0, m_check_overlaps);
    auto gfrp = new G4Tubs("TargetGFRP",
                           gfrp_size[0],
                           gfrp_size[1],
                           gfrp_size[2],
                           0.*deg, 360.*deg);
    auto gfrp_lv = new G4LogicalVolume(gfrp,
                                       mat("G10"),
                                       "TargetGFRPLV");
    gfrp_lv->SetVisAttributes(G4Colour::Green());
    new G4PVPlacement(rot, target_pos, gfrp_lv, "TargetGFRPPV",
                      m_world_lv, true, 0, m_check_overlaps);
  }
    break;
  default:
    G4Exception(FUNC_NAME,
                "Invalid experiment", FatalException,
                ("Found invalid experiment "+std::to_string(m_experiment)
                 +" in "+gConf.Get<G4String>("CONF")).c_str());
    return;
  }
  auto target_lv = new G4LogicalVolume(target_solid,
                                       mat("LD2"),
                                       "TargetLV");
  // target_lv->SetSensitiveDetector(target_sd);
  target_lv->SetVisAttributes(G4Colour::Blue());
  new G4PVPlacement(rot, target_pos,
                    target_lv, "TargetPV",
                    m_world_lv, true, 0, m_check_overlaps);
}
