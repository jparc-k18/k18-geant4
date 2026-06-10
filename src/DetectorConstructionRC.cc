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
#include <G4MultiUnion.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4PVReplica.hh>
#include <G4UniformMagField.hh>
#include <G4Mag_UsualEqRhs.hh>
#include <G4SimpleRunge.hh>
#include <G4VisAttributes.hh>
#include <G4GenericTrap.hh>

#include <cmath>
#include <stdexcept>
#include <vector>

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
  const auto& gGeom = DCGeomMan::GetInstance();
  const auto& gSize = DetSizeMan::GetInstance();
  inline G4Material* mat(const G4String& name) {
    return G4Material::GetMaterial(name, /*warning*/ true);
  }

  constexpr G4int kNumRcModules = 3;
  constexpr G4int kNumRcSegmentsPerModule = NumOfSegRC/kNumRcModules;
  constexpr G4double kRcLayerPitch = -7.0*mm;
  // RC drawing 220707 parts 27-29: top/bottom plates 275.5 x 56 x 3.
  // The plate length is NOT the scint length (scints are 262 long with
  // 5-mm pressers and 1.75-mm play at both ends inside the same cavity).
  constexpr G4double kRcFramePlateLength = 275.5*mm;
  constexpr G4double kRcFramePlateWidth = 56.0*mm;
  constexpr G4double kRcFramePlateThickness = 3.0*mm;
  // Part 26 corner L block: 56 x 30 outline, 7 thick (vertical), arms 11
  // wide; lies flat on the plate outer faces. End-board height check:
  // 130 = 7 + 3 + 110 + 3 + 7.
  constexpr G4double kRcFrameFootThickness = 7.0*mm;
  constexpr G4double kRcFrameFootLength = 30.0*mm;
  constexpr G4double kRcFrameFootLip = 11.0*mm;
  constexpr G4double kRcFrameFootSideLength =
    kRcFrameFootLength - kRcFrameFootLip;   // (19) in drawing 26
  constexpr G4double kRcFrameScintPlateClearance = 0.0*mm;
  // Parts 30/31: Al fiber boards closing both module ends, 110 x 56 x 3.
  constexpr G4double kRcEndBoardThickness = 3.0*mm;
  constexpr G4double kRcEndBoardHeight = 110.0*mm;
  // ---- PD (TD) Al frame, drawing 220707 parts 22-25 ----
  // Assembly local z origin = TD-Y scint mid-plane (DCGEO PDY). All slots
  // are 4 wide in z with a 2-mm wall on the target side; the set screws
  // press TD-Y against the front wall and TD-Z against the back wall
  // ("TD-Y to Z wo hanareru houkou ni oshitsuke"), giving the 0.4-mm
  // surface gap (PDY-PDZ plane distance 3.4 in DCGEO).
  constexpr G4double kPdYFrameOuterHeight = 113.2*mm;
  constexpr G4double kPdYFrameWidth = 7.0*mm;
  constexpr G4double kPdYBottomedFrameDepth = 7.5*mm;  // part 23
  constexpr G4double kPdYThroughFrameDepth = 6.0*mm;   // part 22
  constexpr G4double kPdYFrameAssemblyLength = 275.5*mm;
  constexpr G4double kPdYFrameSlotLength = 92.0*mm;
  constexpr G4double kPdYPitch = 15.0*mm;
  constexpr G4double kPdZPitch = 16.15*mm;
  constexpr G4double kPdFrameSlotOpening = 4.0*mm;
  constexpr G4double kPdFrameSlotFrontWall = 2.0*mm;
  constexpr G4double kPdZFrameLength = 275.5*mm;
  constexpr G4double kPdZFrameSlotLength = 260.5*mm;
  constexpr G4double kPdZFrameDepth = 7.0*mm;          // both rails
  constexpr G4double kPdZTopFrameHeight = 7.0*mm;      // part 25
  constexpr G4double kPdZBottomFrameHeight = 8.0*mm;   // part 24
  constexpr G4double kPdZFrameSlotFloor = 2.0*mm;      // blind slot floor


  inline G4bool has_detector_record(const G4String& name)
  {
    try {
      gGeom.GetDetectorId(name);
      return true;
    } catch(const std::invalid_argument&) {
      return false;
    }
  }

  inline G4String rc_side_name(G4int i)
  {
    return (i == 1) ? "RC-X" : "RC+X";
  }

  inline G4String legacy_rc_name(G4int i)
  {
    return "RC" + std::to_string(i);
  }

  inline G4String legacy_rc_record_prefix(G4int i)
  {
    return (i == 1) ? "RC-M" : "RC-P";
  }

  inline G4String rc_record_prefix(G4int i)
  {
    const auto side = rc_side_name(i);
    if(has_detector_record(side + "-RC1") || has_detector_record(side + "-PDY"))
      return side;

    const auto legacy_side = legacy_rc_record_prefix(i);
    if(has_detector_record(legacy_side + "-RC1") ||
       has_detector_record(legacy_side + "-PDY"))
      return legacy_side;

    return "";
  }

  inline G4String rc_module_record(const G4String& prefix, G4int module)
  {
    return prefix + "-RC" + std::to_string(module + 1);
  }

  inline G4String pd_record_name(const G4String& prefix, const G4String& yz)
  {
    return prefix + "-PD" + yz;
  }

  inline G4String pd_branch_name(G4int i, const G4String& yz)
  {
    return rc_side_name(i) + "-PD" + yz;
  }

  inline void apply_dcgeom_rotation(G4RotationMatrix* rot,
                                    const G4String& record_name)
  {
    rot->rotateY(-gGeom.GetRotAngle2(record_name) * deg);
  }

  inline void apply_rc_side_rotation(G4RotationMatrix* rot,
                                     const G4String& prefix)
  {
    rot->rotateY(-gGeom.GetRotAngle2(rc_module_record(prefix, 0)) * deg);
  }

  // Vertical TD-Y post (parts 22/23). Local z origin of the returned
  // solid is the post center; post_z_offset returns where that center
  // must sit relative to the TD-Y scint mid-plane so that the slot keeps
  // its 2-mm front wall with the scint face on it.
  G4VSolid*
  make_pd_y_frame_solid(const G4String& name,
                        const G4ThreeVector& seg_size,
                        G4bool bottomed_frame,
                        G4double& post_z_offset)
  {
    const G4double depth = bottomed_frame ? kPdYBottomedFrameDepth
                                          : kPdYThroughFrameDepth;
    const G4double slot_front_z = -seg_size.z();
    post_z_offset = slot_front_z - kPdFrameSlotFrontWall + 0.5*depth;
    auto outer_solid = new G4Box(name + "OuterSolid",
                                 0.5*kPdYFrameWidth,
                                 0.5*kPdYFrameOuterHeight,
                                 0.5*depth);
    // The through frame (6 deep) has no back wall: 2 + 4 reaches the back
    // face, so extend the cut past it for a clean subtraction.
    const G4double slot_back_z = bottomed_frame
      ? slot_front_z + kPdFrameSlotOpening
      : post_z_offset + 0.5*depth + 0.5*mm;
    const G4double slot_depth_x =
      bottomed_frame ? kPdYFrameWidth - 0.25*mm
                     : kPdYFrameWidth + 1.*mm;
    auto slot_solid = new G4Box(name + "ScintSlotSolid",
                                0.5*slot_depth_x,
                                0.5*kPdYFrameSlotLength,
                                0.5*(slot_back_z - slot_front_z));
    const G4double slot_center_x = bottomed_frame
      ? -0.5*kPdYFrameWidth + 0.5*slot_depth_x
      : 0.0;
    return new G4SubtractionSolid(
      name + "Solid", outer_solid, slot_solid, nullptr,
      G4ThreeVector(slot_center_x, 0.0,
                    0.5*(slot_front_z + slot_back_z) - post_z_offset));
  }

  // Horizontal TD-Z rail (parts 24/25): 7 deep, 7 (top) / 8 (bottom)
  // tall, with a blind slot from the window-side face leaving a 2-mm
  // floor. rail_z_offset returns the rail-center position relative to
  // the TD-Z scint mid-plane (scint pressed onto the slot back wall).
  G4VSolid*
  make_pd_z_frame_solid(const G4String& name,
                        const G4ThreeVector& seg_size,
                        G4bool bottom_rail,
                        G4double& rail_z_offset)
  {
    const G4double height = bottom_rail ? kPdZBottomFrameHeight
                                        : kPdZTopFrameHeight;
    rail_z_offset = seg_size.z() - kPdFrameSlotOpening
      - kPdFrameSlotFrontWall + 0.5*kPdZFrameDepth;
    auto outer_solid = new G4Box(name + "OuterSolid",
                                 0.5*kPdZFrameLength,
                                 0.5*height,
                                 0.5*kPdZFrameDepth);
    const G4double slot_depth = height - kPdZFrameSlotFloor;
    const G4double inner_face_y = bottom_rail ? 0.5*height : -0.5*height;
    const G4double dir = bottom_rail ? -1. : 1.;  // into the rail
    auto slot_solid = new G4Box(name + "SlotSolid",
                                0.5*kPdZFrameSlotLength,
                                0.5*slot_depth + 0.25*mm,
                                0.5*kPdFrameSlotOpening);
    const G4double slot_center_y =
      inner_face_y + dir*(0.5*slot_depth - 0.25*mm);
    const G4double slot_center_z =
      seg_size.z() - 0.5*kPdFrameSlotOpening - rail_z_offset;
    return new G4SubtractionSolid(
      name + "Solid", outer_solid, slot_solid, nullptr,
      G4ThreeVector(0.0, slot_center_y, slot_center_z));
  }

  G4VSolid*
  make_pd_frame_assembly_solid(const G4String& name,
                               const G4ThreeVector& pd_y_seg_size,
                               const G4ThreeVector& pd_z_seg_size,
                               const G4ThreeVector& pd_z_local_offset)
  {
    auto assembly_solid = new G4MultiUnion(name + "Solid");
    for(G4int iside=0; iside<2; ++iside){
      const G4double sign = (iside == 0) ? -1. : 1.;
      const G4bool bottomed_frame = (iside == 1);
      G4double post_z_offset = 0.;
      auto side_solid =
        make_pd_y_frame_solid(name + (bottomed_frame ? "YBottomed"
                                                     : "YThrough"),
                              pd_y_seg_size, bottomed_frame,
                              post_z_offset);
      const G4ThreeVector side_pos(
        sign*(0.5*kPdYFrameAssemblyLength - 0.5*kPdYFrameWidth),
        0.0,
        post_z_offset);
      assembly_solid->AddNode(*side_solid,
                              G4Transform3D(G4RotationMatrix(), side_pos));
    }

    for(G4int iside=0; iside<2; ++iside){
      const G4bool bottom_rail = (iside == 0);
      const G4double ysign = bottom_rail ? -1. : 1.;
      const G4double rail_height =
        bottom_rail ? kPdZBottomFrameHeight : kPdZTopFrameHeight;
      G4double rail_z_offset = 0.;
      auto rail_solid =
        make_pd_z_frame_solid(name + (bottom_rail ? "ZBottom" : "ZTop"),
                              pd_z_seg_size, bottom_rail, rail_z_offset);
      // Rails are flush with the post ends: 7 + 98.2 + 8 = 113.2.
      const G4ThreeVector rail_pos(
        0.0,
        ysign*(0.5*kPdYFrameOuterHeight - 0.5*rail_height),
        pd_z_local_offset.z() + rail_z_offset);
      assembly_solid->AddNode(*rail_solid,
                              G4Transform3D(G4RotationMatrix(), rail_pos));
    }
    assembly_solid->Voxelize();
    return assembly_solid;
  }

}

//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructRC(G4int i)
{
  const G4String sd_name = rc_side_name(i);
  const auto& half_size = gSize.GetSize("RcSeg")*mm/2;
  const G4String prefix = rc_record_prefix(i);

  if(!prefix.empty() && has_detector_record(rc_module_record(prefix, 0))){
    auto sdRC = new RCSD(sd_name);
    AddNewDetector(sdRC);

    auto segment_solid = new G4Box("RcSegmentSolid",
                                   half_size.x(), half_size.y(), half_size.z());
    auto segment_lv = new G4LogicalVolume(segment_solid,
                                          mat("Scintillator"),
                                          "RcSegmentLV");
    segment_lv->SetVisAttributes(G4Colour::Cyan());
    segment_lv->SetSensitiveDetector(sdRC);

    // RC drawing 220707 parts 27-29: 3-mm top/bottom plates (275.5 long,
    // independent of the 262-mm scint length).
    auto frame_plate_solid = new G4Box("RcAlFramePlateSolid",
                                       0.5*kRcFramePlateLength,
                                       0.5*kRcFramePlateThickness,
                                       0.5*kRcFramePlateWidth);
    auto frame_plate_lv = new G4LogicalVolume(frame_plate_solid, mat("Al"),
                                              "RcAlFramePlateLV");
    frame_plate_lv->SetVisAttributes(G4Colour::Gray());
    // Part 26 corner L block (56 x 30 x 7t, arms 11 wide) lying flat on
    // the plate outer face: 11 x 56 end bar + 11 x 19 side arm.
    auto frame_foot_end_solid = new G4Box("RcAlFrameFootEndSolid",
                                          0.5*kRcFrameFootLip,
                                          0.5*kRcFrameFootThickness,
                                          0.5*kRcFramePlateWidth);
    auto frame_foot_end_lv =
      new G4LogicalVolume(frame_foot_end_solid, mat("Al"),
                          "RcAlFrameFootEndLV");
    frame_foot_end_lv->SetVisAttributes(G4Colour::Gray());
    auto frame_foot_side_solid = new G4Box("RcAlFrameFootSideSolid",
                                           0.5*kRcFrameFootSideLength,
                                           0.5*kRcFrameFootThickness,
                                           0.5*kRcFrameFootLip);
    auto frame_foot_side_lv =
      new G4LogicalVolume(frame_foot_side_solid, mat("Al"),
                          "RcAlFrameFootSideLV");
    frame_foot_side_lv->SetVisAttributes(G4Colour::Gray());
    // Parts 30/31: fiber boards closing both ends of each module.
    auto end_board_solid = new G4Box("RcEndBoardSolid",
                                     0.5*kRcEndBoardThickness,
                                     0.5*kRcEndBoardHeight,
                                     0.5*kRcFramePlateWidth);
    auto end_board_lv = new G4LogicalVolume(end_board_solid, mat("Al"),
                                            "RcEndBoardLV");
    end_board_lv->SetVisAttributes(G4Colour::Gray());
    // // Part 1: readout PCB flush on the +x end board.
    // auto readout_solid = new G4Box("RcReadoutBoardSolid",
    //                                0.5*kRcReadoutBoardThickness,
    //                                0.5*kRcEndBoardHeight,
    //                                0.5*kRcReadoutBoardWidth);
    // auto readout_lv = new G4LogicalVolume(readout_solid, mat("G10"),
    //                                       "RcReadoutBoardLV");
    // readout_lv->SetVisAttributes(G4Colour::Green());

    std::vector<G4ThreeVector> module_positions;
    for(G4int imod=0; imod<kNumRcModules; ++imod){
      module_positions.push_back(
        gGeom.GetGlobalPosition(rc_module_record(prefix, imod)));
    }

    for(G4int imod=0; imod<kNumRcModules; ++imod){
      const auto det_name = rc_module_record(prefix, imod);
      const auto pitch = kRcLayerPitch;
      const auto module_pos = module_positions[imod];
      const G4double plate_center_y =
        half_size.y() + kRcFrameScintPlateClearance
        + 0.5*kRcFramePlateThickness;
      const G4double foot_center_y =
        half_size.y() + kRcFrameScintPlateClearance
        + kRcFramePlateThickness + 0.5*kRcFrameFootThickness;

      for(G4int iplate=0; iplate<2; ++iplate){
        auto frame_rot = new G4RotationMatrix;
        apply_dcgeom_rotation(frame_rot, det_name);
        const G4double sign = (iplate == 0) ? -1. : 1.;
        const G4ThreeVector local_pos(
          0.0,
          sign*plate_center_y,
          0.0);
        const auto global_pos = module_pos + (*frame_rot)*local_pos;
        new G4PVPlacement(
          frame_rot, global_pos, frame_plate_lv,
          sd_name + (iplate == 0 ? "AlBottomPlatePHYS" : "AlTopPlatePHYS"),
          m_world_lv, false, imod, m_check_overlaps);

        for(G4int iend=0; iend<2; ++iend){
          const G4double xsign = (iend == 0) ? -1. : 1.;
          auto foot_end_rot = new G4RotationMatrix;
          apply_dcgeom_rotation(foot_end_rot, det_name);
          const G4ThreeVector foot_end_local_pos(
            xsign*(0.5*kRcFramePlateLength - 0.5*kRcFrameFootLip),
            sign*foot_center_y,
            0.0);
          const auto foot_end_global_pos =
            module_pos + (*foot_end_rot)*foot_end_local_pos;
          new G4PVPlacement(
            foot_end_rot, foot_end_global_pos, frame_foot_end_lv,
            sd_name + (iplate == 0 ? "AlBottomFootEndPHYS"
                                   : "AlTopFootEndPHYS"),
            m_world_lv, false, 2*imod + iend, m_check_overlaps);

          auto foot_side_rot = new G4RotationMatrix;
          apply_dcgeom_rotation(foot_side_rot, det_name);
          const G4double zsign = (iplate == 0) ? 1. : -1.;
          const G4ThreeVector foot_side_local_pos(
            xsign*(0.5*kRcFramePlateLength - kRcFrameFootLip
                   - 0.5*kRcFrameFootSideLength),
            sign*foot_center_y,
            zsign*(0.5*kRcFramePlateWidth - 0.5*kRcFrameFootLip));
          const auto foot_side_global_pos =
            module_pos + (*foot_side_rot)*foot_side_local_pos;
          new G4PVPlacement(
            foot_side_rot, foot_side_global_pos, frame_foot_side_lv,
            sd_name + (iplate == 0 ? "AlBottomFootSidePHYS"
                                   : "AlTopFootSidePHYS"),
            m_world_lv, false, 4*imod + 2*iplate + iend,
            m_check_overlaps);
        }
      }

      // Fiber boards (parts 30/31) on both ends (inside the top/bottom plates);
      // readout PCB (part 1) flush on the +x end board (fiber-readout end, assumed +x local).
      for(G4int iend=0; iend<2; ++iend){
        const G4double xsign = (iend == 0) ? -1. : 1.;
        auto board_rot = new G4RotationMatrix;
        apply_dcgeom_rotation(board_rot, det_name);
        const G4ThreeVector board_local_pos(
          xsign*(0.5*kRcFramePlateLength - 0.5*kRcEndBoardThickness),
          0.0, 0.0);
        new G4PVPlacement(
          board_rot, module_pos + (*board_rot)*board_local_pos,
          end_board_lv, sd_name + "EndBoardPHYS",
          m_world_lv, false, 2*imod + iend, m_check_overlaps);
      }
      // {
      //   auto pcb_rot = new G4RotationMatrix;
      //   apply_dcgeom_rotation(pcb_rot, det_name);
      //   const G4ThreeVector pcb_local_pos(
      //     0.5*kRcFramePlateLength + 0.5*kRcReadoutBoardThickness,
      //     0.0, 0.0);
      //   new G4PVPlacement(
      //     pcb_rot, module_pos + (*pcb_rot)*pcb_local_pos, readout_lv,
      //     sd_name + "ReadoutBoardPHYS",
      //     m_world_lv, false, imod, m_check_overlaps);
      // }


      for(G4int iseg=0; iseg<kNumRcSegmentsPerModule; ++iseg){
        auto rot = new G4RotationMatrix;
        apply_dcgeom_rotation(rot, det_name);
        const G4ThreeVector local_pos(
          0.0, 0.0,
          (-kNumRcSegmentsPerModule/2. + iseg)*pitch + pitch/2.);
        const G4int copy_no = imod*kNumRcSegmentsPerModule + iseg;
        const auto global_pos = module_pos + (*rot)*local_pos;
        new G4PVPlacement(rot, global_pos, segment_lv,
                          sd_name + "SegmentPHYS", m_world_lv, false,
                          copy_no, m_check_overlaps);
      }
    }

    G4cout << "[RC] constructed " << sd_name << " from "
           << prefix << "-RC1..3"
           << " using individual DCGEO module centers" << G4endl;
    return;
  }

  const G4String geom_name = legacy_rc_name(i);
  const auto& ra2 = gGeom.GetRotAngle2(geom_name) * deg;
  const G4double pitch = gGeom.GetWirePitch(geom_name)*mm;
  auto sdRC = new RCSD(sd_name);
  AddNewDetector(sdRC);
  auto mother_solid = new G4Box("RcMotherSolid",
                                half_size.x()*mm,
                                half_size.y()*mm,
                                half_size.z()*NumOfSegRC*mm);
  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       mat("Air"),
                                       "RcMotherLV");
  auto rot = new G4RotationMatrix;
  rot->rotateY(-ra2);
  auto pos = gGeom.GetGlobalPosition(geom_name);
  new G4PVPlacement(rot, pos, mother_lv, "RcMotherPHYS", m_world_lv, false,
                    0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());

  auto segment_solid = new G4Box("RcSegmentSolid", half_size.x(),
                                 half_size.y(), half_size.z());
  auto segment_lv = new G4LogicalVolume(segment_solid,
                                        mat("Scintillator"),
                                        "RcSegmentLV");
  segment_lv->SetVisAttributes(G4Colour::Cyan());
  segment_lv->SetSensitiveDetector(sdRC);
  for(G4int iseg=0; iseg<NumOfSegRC; ++iseg){
    pos = G4ThreeVector(0.0,
                        0.0,
                        (-NumOfSegRC/2 + iseg)*pitch + pitch/2);
    new G4PVPlacement(nullptr, pos, segment_lv,
                      sd_name + "SegmentPHYS", mother_lv, false,
                      iseg, m_check_overlaps);
  }
  G4cout << "[RC] constructed " << sd_name << " from legacy "
         << geom_name << G4endl;
}


//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructPD(G4int i, G4String YZ)
{
  const G4String name = pd_branch_name(i, YZ);
  const G4String prefix = rc_record_prefix(i);
  const G4String detailed_name = prefix.empty() ? "" : pd_record_name(prefix, YZ);

  if(!detailed_name.empty() && has_detector_record(detailed_name)){
    const auto& seg_size = gSize.GetSize("Pd"+YZ+"Seg")*mm/2;
    G4int NumOfSegPD = (YZ == "Y") ? NumOfSegPDY : NumOfSegPDZ;
    const G4double pitch = (YZ == "Y") ? kPdYPitch : kPdZPitch;

    auto sdPD = new PDSD(name);
    AddNewDetector(sdPD);

    auto segment_solid = new G4Box("PdSegmentSolid",
                                   seg_size.x(), seg_size.y(), seg_size.z());
    auto segment_lv = new G4LogicalVolume(segment_solid,
                                          mat("Scintillator"),
                                          "PdSegmentLV");
    segment_lv->SetVisAttributes(YZ == "Y" ? G4Colour::Green()
                                            : G4Colour::Magenta());
    segment_lv->SetSensitiveDetector(sdPD);

    const auto det_pos = gGeom.GetGlobalPosition(detailed_name);

    if(YZ == "Y"){
      const auto& pd_z_seg_size = gSize.GetSize("PdZSeg")*mm/2;
      const G4String pd_z_name = pd_record_name(prefix, "Z");
      const auto pd_z_det_pos = gGeom.GetGlobalPosition(pd_z_name);
      auto frame_rot = new G4RotationMatrix;
      apply_rc_side_rotation(frame_rot, prefix);
      // G4PVPlacement maps a frame-local point l to det_pos + R^-1 * l,
      // so the world offset d must enter the solid as l = R * d.
      const auto pd_z_local_offset =
        (*frame_rot)*(pd_z_det_pos - det_pos);
      auto frame_solid =
        make_pd_frame_assembly_solid(name + "AlFrameAssembly",
                                     seg_size, pd_z_seg_size,
                                     pd_z_local_offset);
      auto frame_lv = new G4LogicalVolume(frame_solid, mat("Al"),
                                          name + "AlFrameAssemblyLV");
      frame_lv->SetVisAttributes(G4Colour::Gray());
      new G4PVPlacement(frame_rot, det_pos, frame_lv,
                        name + "AlFrameAssemblyPHYS",
                        m_world_lv, false, 0, m_check_overlaps);
    }

    for(G4int iseg=0; iseg<NumOfSegPD; ++iseg){
      auto rot = new G4RotationMatrix;
      apply_rc_side_rotation(rot, prefix);
      G4ThreeVector local_pos;
      if(YZ == "Y"){
        local_pos = G4ThreeVector(0.0,
                                  (-NumOfSegPD/2 + iseg)*pitch + pitch/2,
                                  0.0);
      } else {
        local_pos = G4ThreeVector((-NumOfSegPD/2 + iseg)*pitch + pitch/2,
                                  0.0, 0.0);
      }
      const auto global_pos = det_pos + (*rot)*local_pos;
      new G4PVPlacement(rot, global_pos, segment_lv,
                        name + "SegmentPHYS", m_world_lv, false,
                        iseg, m_check_overlaps);
    }

    G4cout << "[PD] constructed " << name << " from "
           << detailed_name << G4endl;
    return;
  }

  const G4String RCname = legacy_rc_name(i);
  const auto& ra2 = gGeom.GetRotAngle2(RCname) * deg;
  const auto& half_size = gSize.GetSize("Pd"+YZ)*mm/2;
  const auto& seg_size = gSize.GetSize("Pd"+YZ+"Seg")*mm/2;
  const auto& Rcseg_size = gSize.GetSize("RcSeg")*mm/2;
  G4double pitch;
  G4int NumOfSegPD = 0;
  auto sdPD = new PDSD(name);
  AddNewDetector(sdPD);
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
  new G4PVPlacement(rot, pos, mother_lv, name + "MotherPHYS", m_world_lv, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());

  auto segment_solid = new G4Box("PdSegmentSolid", seg_size.x(),
                                 seg_size.y(), seg_size.z());
  auto segment_lv = new G4LogicalVolume(segment_solid,
                                        mat("Scintillator"),
                                        "PdSegmentLV");
  segment_lv->SetSensitiveDetector(sdPD);
  for(G4int iseg=0; iseg<NumOfSegPD; ++iseg){
    if(YZ == "Y"){
      segment_lv->SetVisAttributes(G4Colour::Green());
      pos = G4ThreeVector(0.0, (-NumOfSegPD/2 + iseg)*pitch + pitch/2, 0.0);
    }
    else{
      segment_lv->SetVisAttributes(G4Colour::Magenta());
      pos = G4ThreeVector((-NumOfSegPD/2 + iseg)*pitch + pitch/2, 0.0, 0.0);
    }
    new G4PVPlacement(nullptr, pos, segment_lv,
                      name + "SegmentPHYS", mother_lv, false, iseg, m_check_overlaps);
  }
}
