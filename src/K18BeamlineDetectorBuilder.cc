// -*- C++ -*-

#include "K18BeamlineDetectorBuilder.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <G4Box.hh>
#include <G4Colour.hh>
#include <G4Element.hh>
#include <G4ExtrudedSolid.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4PVPlacement.hh>
#include <G4RotationMatrix.hh>
#include <G4SDManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4Tubs.hh>
#include <G4TwoVector.hh>
#include <G4VSolid.hh>
#include <G4VSensitiveDetector.hh>
#include <G4VisAttributes.hh>

#include "ConfMan.hh"
#include "DCSD.hh"
#include "DCGeomMan.hh"
#include "DetSizeMan.hh"
#include "K18BeamlineField.hh"
#include "TOFSD.hh"

namespace
{
const auto& confMan = ConfMan::GetInstance();
const auto& geomMan = DCGeomMan::GetInstance();
const auto& sizeMan = DetSizeMan::GetInstance();

const G4double kDeg2Rad = std::acos(-1.)/180.;
const G4double kRhoD4 = 4.000*m;
const G4double kD4SpaceX = 0.657*m;
const G4double kD4HalfGap = 0.1*m;
const G4double kQ10x = 0.89*m;
const G4double kQ10y = 0.89*m;
const G4double kQ10z = 0.956*m;
const G4double kQ11x = 0.89*m;
const G4double kQ11y = 0.89*m;
const G4double kQ11z = 0.66*m;
const G4double kQ12x = 1.095*m;
const G4double kQ12y = 1.095*m;
const G4double kQ12z = 0.66*m;
const G4double kQ13x = 1.095*m;
const G4double kQ13y = 1.095*m;
const G4double kQ13z = 0.66*m;
const G4double kDriftL1 = 0.242*m;
const G4double kDriftL2 = 0.470*m;
const G4double kDriftL3 = 0.470*m;
const G4double kDriftL4 = 0.240*m;

// Dimensions from the QQDQQ generation scripts, not STEP/drawings:
//   NQ418_v25_native_curves_layered_autorun.py   (Q10)
//   NQ412_v1_fusion_autorun.py                   (Q11/Q12/Q13; yoke half-z=300)
//   K18D4_fusion_script_v32_end_transition.py    (D4 v54)
//   K18_QQDQQ_clean_assembly.py                  (layout, Y-mirror)
// Magnet *centres* stay on the field-layout chain above.
const G4double kQYokeHX = 450.*mm;
const G4double kQYokeHY = 650.*mm;
const G4double kQYokeInner = 390.*mm;
const G4double kQBoreR = 100.*mm;
const G4double kQPoleLimb62 = 62.*mm;
const G4double kQPoleLimb107 = 107.*mm;
const G4double kQPoleBreak = 225.*mm;
const G4double kQPoleTaperZ = 45.*mm;
const G4int kQPoleTaperSlices = 5;
const G4double kQGuardThick = 28.*mm;
const G4double kQGuardStarArm = 170.*mm;
const G4double kQGuardStarNeckX = 29.291*mm;
const G4double kQGuardStarNeckY = 29.412*mm;
const G4double kQGuardStarArcR = 141.214*mm;
const G4double kQGuardStarCx = 170.503*mm;
const G4double kQGuardStarCy = 170.625*mm;
// Beam duct: star = DXF _1__S opening (slightly inset), with the D4 rectangle
// inside the 100 mm half-gap.
const G4double kPipeStarScale = 0.97;
const G4double kPipeD4HalfW = 150.*mm; // radial about R4000
const G4double kPipeD4HalfH = 80.*mm;  // vertical; kD4HalfGap = 100 mm
const G4double kPipeDiskThick = 12.*mm;
const G4double kPipeDiskHoleR = 175.*mm; // covers star arms and rect corners
const G4double kPipeDiskOuterR = 220.*mm;
const G4double kPipeDiskStandoff = 30.*mm; // Q-endguard outer to disk inner face
// NQ418/NQ412 coil hull: one filled pancake + endturn + R-bend per
// quadrant, using the outer contour of the 30 native saddle turns.
const G4double kQCoilHalf = 8.*mm;
const G4double kQCoilClearZ = 1.*mm;
const G4double kQCoilPitchZ = 16.*mm;
const G4double kQCoilSlotR0 = 243.*mm;
const G4double kQCoilSlotPitch = 17.*mm;
const G4double kQCoilSlotCount = 8.;
const G4double kQCoilDInner = 192.*mm;
const G4double kQCoilDOuter = 226.*mm;
const G4double kQCoilCorner = 108.*mm;
const G4double kQCoilDBendR = 50.*mm;
const G4double kQ10IronHZ = 450.*mm;
const G4double kNQ412IronHZ = 300.*mm;
const G4double kQ10GuardTerm = (450. + 72. + 28.)*mm;
const G4double kQ10GuardNon = -(450. + 102. + 28.)*mm;
const G4double kNQ412GuardTerm = (300. + 72. + 28.)*mm;
const G4double kNQ412GuardNon = -(300. + 102. + 28.)*mm;
const G4double kD4YokeR0 = 3040.*mm;
const G4double kD4YokeR1 = 4850.*mm;
const G4double kD4YokeHZ = 505.*mm;
const G4double kD4PoleR0 = 3670.*mm;
const G4double kD4PoleR1 = 4330.*mm;
const G4double kD4SlotR0 = 3520.*mm;
const G4double kD4SlotR1 = 4480.*mm;
const G4double kD4SlotHZ = 83.*mm;
const G4double kD4LegZ0 = 1.*mm;
const G4double kD4LegZ1 = 83.*mm;
const G4double kD4LegZcen = 0.5*(kD4LegZ0 + kD4LegZ1);
const G4int kD4YokeRadialSlices = 36;
const G4double kD4RiserCX = 75.*mm;
const G4double kD4RiserCZ = 183.*mm;
// Fusion v32 pack-wrap of the 0-deg riser (RISER_R_IN/OUT).  Join to the
// pancake is the vertical at y=CX, z=(CZ-RIn)..(CZ-ROut) = 83..1.
// The 4-layer turn centerlines (CX 75..-13.1) describe conductors rather than
// hull edges; the pack-wrap values preserve the correct vertical join.
const G4double kD4RiserRIn = 100.*mm;
const G4double kD4RiserROut = 182.*mm;
const G4double kD4BendLayerR[4] = {110.4*mm, 130.8*mm, 151.2*mm, 171.5*mm};
const G4double kD4BendLayerCX[4] = {75.0*mm, 45.7*mm, 16.4*mm, -13.1*mm};
const G4double kD4ExitYNear = 24.*mm;
const G4double kD4ExitYFar = 120.*mm;
const G4double kD4ArchCS = 230.*mm;
const G4double kD4ArchCZ = 221.*mm;
const G4double kD4ArchAx = 82.*mm; // = ROut-RIn; Fusion ARCH_AX
const G4double kD4EntLegOverlap = 0.5*mm; // Fusion ENTRANCE_LEG_OVERLAP
const G4int kD4CoilPackSlices = 12;
const G4int kD4CoilLayers = 4;
const G4double kD4MainW = 22.5*mm;
const G4double kD4MainCin = kRhoD4 - 377.5*mm;
const G4double kD4MainCout = kRhoD4 + 377.5*mm;
const G4double kD4MainRoff[4] = {-35.25*mm, -11.75*mm, 11.75*mm, 35.25*mm};
const G4double kD4AuxW = 11.*mm;
const G4double kD4AuxCin = kRhoD4 - 449.5*mm;
const G4double kD4AuxCout = kRhoD4 + 449.5*mm;
const G4double kD4AuxRoff[4] = {-18.*mm, -6.*mm, 6.*mm, 18.*mm};

// Fusion exit-face arc ends (v32).  After the clean-assembly Y-mirror these
// sit on the Q11/upstream side: G4 phi = 64° - Fusion phi.
const G4double kD4PhiFaceIn = 72.1465*deg;
const G4double kD4PhiFaceOut = 59.70*deg;
const G4double kD4Alpha = 23.5*deg;

enum class K18Reference { VI, VO };

struct K18PlaneSpec
{
  const char* name;
  G4double half_x;
  G4double half_y;
  G4double half_z;
};

struct K18PassivePlaneSpec
{
  G4String name;
  G4double local_l;
  G4double half_x;
  G4double half_y;
  G4double half_z;
  G4Material* material;
};

G4bool conf_bool_or(const G4String& key, G4bool fallback)
{
  const G4String raw = confMan.Get<G4String>(key);
  return raw.empty() ? fallback : confMan.Get<G4bool>(key);
}

G4double conf_double_or(const G4String& key, G4double fallback)
{
  const G4String raw = confMan.Get<G4String>(key);
  return raw.empty() ? fallback : confMan.Get<G4double>(key);
}

G4String conf_string_or(const G4String& key, const G4String& fallback)
{
  const G4String raw = confMan.Get<G4String>(key);
  return raw.empty() ? fallback : raw;
}

G4double d4_tan_half()
{
  return kRhoD4*std::tan(0.5*K18BeamlineFrame::BendAngle()*kDeg2Rad);
}

G4double target_l()
{
  try {
    return geomMan.GetLocalZ("K18Target");
  } catch(...) {
    const G4double val = confMan.Get<G4double>("K18TargetL");
    return (val != 0.) ? val : 1503.;
  }
}

G4Material* k18_vacuum_material()
{
  const G4double pressure_pa = confMan.Get<G4double>("K18VacuumPressurePa");
  if(pressure_pa <= 0.){
    auto* vacuum = G4Material::GetMaterial("Vacuum", true);
    G4cout << "[K18BeamlineDetectorBuilder] K18 vacuum: legacy Vacuum "
           << "material (set K18VacuumPressurePa > 0 to use pressure gas)"
           << G4endl;
    return vacuum;
  }

  const G4double pressure = pressure_pa*pascal;
  const G4double temperature = NTP_Temperature;
  const G4double density =
    (1.290*mg/cm3)*(pressure/atmosphere)*(CLHEP::STP_Temperature/temperature);

  std::ostringstream oss;
  oss << "K18Vacuum_" << std::setprecision(6) << pressure_pa << "Pa";
  auto name = oss.str();
  std::replace(name.begin(), name.end(), '.', 'p');
  std::replace(name.begin(), name.end(), '+', 'p');
  std::replace(name.begin(), name.end(), '-', 'm');

  if(auto* material = G4Material::GetMaterial(name, false)){
    G4cout << "[K18BeamlineDetectorBuilder] K18 vacuum: reuse " << name
           << " pressure=" << pressure_pa << " Pa density="
           << density/(g/cm3) << " g/cm3" << G4endl;
    return material;
  }

  auto* material =
    new G4Material(name, density, 2, kStateGas, temperature, pressure);
  material->AddElement(G4Element::GetElement("Nitrogen", true), 70.*perCent);
  material->AddElement(G4Element::GetElement("Oxygen", true), 30.*perCent);

  G4cout << "[K18BeamlineDetectorBuilder] K18 vacuum: " << name
         << " pressure=" << pressure_pa << " Pa density="
         << density/(g/cm3) << " g/cm3" << G4endl;
  return material;
}

G4RotationMatrix* make_rotation_from_internal(const K18BeamlineFrame& frame,
                                             const G4RotationMatrix& rot_native)
{
  auto* rot = new G4RotationMatrix;
  const auto x_axis = frame.InternalToS2SVector(
    rot_native*G4ThreeVector(1., 0., 0.)).unit();
  const auto y_axis = frame.InternalToS2SVector(
    rot_native*G4ThreeVector(0., 1., 0.)).unit();
  const auto z_axis = frame.InternalToS2SVector(
    rot_native*G4ThreeVector(0., 0., 1.)).unit();
  rot->rotateAxes(x_axis, y_axis, z_axis);
  // rotateAxes defines the daughter axes; G4PVPlacement stores the inverse.
  rot->invert();
  return rot;
}

G4RotationMatrix* make_axis_rotation(const K18BeamlineFrame& frame,
                                     G4double axis_angle,
                                     G4double roll = 0.)
{
  const auto ex = K18BeamlineFrame::XAxis(axis_angle);
  const auto ey = K18BeamlineFrame::YAxis();
  const auto ez = K18BeamlineFrame::Axis(axis_angle);
  const auto x_axis = frame.InternalToS2SVector(
    std::cos(roll)*ex + std::sin(roll)*ey).unit();
  const auto y_axis = frame.InternalToS2SVector(
    -std::sin(roll)*ex + std::cos(roll)*ey).unit();
  const auto z_axis = frame.InternalToS2SVector(ez).unit();
  auto* rot = new G4RotationMatrix;
  rot->rotateAxes(x_axis, y_axis, z_axis);
  rot->invert();
  return rot;
}

G4ThreeVector magnet_offset(G4double axis_angle, const G4ThreeVector& local)
{
  return local.x()*K18BeamlineFrame::XAxis(axis_angle)
    + local.y()*K18BeamlineFrame::YAxis()
    + local.z()*K18BeamlineFrame::Axis(axis_angle);
}

G4RotationMatrix d4_rotation()
{
  G4RotationMatrix rot;
  rot.rotateZ(-90.*deg);
  return rot;
}

G4ThreeVector on_straight_axis(G4double l)
{
  return G4ThreeVector(l, 0., 0.);
}

G4ThreeVector on_bent_axis(G4double l)
{
  return l*K18BeamlineFrame::Axis(K18BeamlineFrame::BendAngle());
}

G4double axis_angle(K18Reference ref)
{
  return (ref == K18Reference::VI) ? 0. : K18BeamlineFrame::BendAngle();
}

G4ThreeVector native_from_reference(K18Reference ref, G4double l)
{
  if(ref == K18Reference::VI)
    return K18BeamlineFrame::Vin() + l*K18BeamlineFrame::Axis(0.);
  return K18BeamlineFrame::Vout()
    + l*K18BeamlineFrame::Axis(K18BeamlineFrame::BendAngle());
}

G4ThreeVector native_from_dcgeo(K18Reference ref, const G4String& name)
{
  const auto global = geomMan.GetGlobalPosition(name)*mm;
  const auto angle = axis_angle(ref);
  return native_from_reference(ref, geomMan.GetLocalZ(name)*mm)
    + global.x()*K18BeamlineFrame::XAxis(angle)
    + global.y()*K18BeamlineFrame::YAxis();
}

std::vector<G4double> bh2_segment_widths()
{
  std::vector<G4double> widths;
  try {
    for(G4int i=0; i<8; ++i)
      widths.push_back(sizeMan.Get("Bh2SegWidth", i)*mm);
  } catch(...) {
    widths.clear();
  }
  if(widths.empty())
    widths.push_back(sizeMan.GetSize("Bh2Seg").x()*mm);
  return widths;
}

std::vector<G4double> bh1_segment_widths()
{
  std::vector<G4double> widths;
  try {
    for(G4int i=0; i<11; ++i)
      widths.push_back(sizeMan.Get("Bh1SegWidth", i)*mm);
  } catch(...) {
    widths.clear();
  }
  if(widths.empty())
    widths = {30.*mm, 20.*mm, 16.*mm, 12.*mm, 8.*mm, 8.*mm,
              8.*mm, 12.*mm, 16.*mm, 20.*mm, 30.*mm};
  return widths;
}

G4RotationMatrix* make_local_stereo_rotation(G4double tilt_angle)
{
  if(std::abs(tilt_angle) < 1.e-9)
    return nullptr;
  auto* rot = new G4RotationMatrix;
  rot->rotateZ(tilt_angle*deg);
  // The pointer-form G4PVPlacement constructor takes the frame rotation,
  // i.e. the inverse of the intended daughter-object rotation.  Keep the
  // plane-local x axis aligned with the DCGEO readout coordinate s.
  rot->invert();
  return rot;
}

void place_box(const G4String& name,
               G4double half_x,
               G4double half_y,
               G4double half_z,
               const G4ThreeVector& native_center,
               G4double axis_angle,
               const K18BeamlineFrame& frame,
               G4LogicalVolume* world_lv,
               G4bool check_overlaps,
               G4VisAttributes* vis,
               G4Material* material,
               G4VSensitiveDetector* sd = nullptr,
               G4int copy_no = 0)
{
  auto* solid = new G4Box(name + "Solid", half_x, half_y, half_z);
  auto* lv = new G4LogicalVolume(solid, material, name + "LV");
  lv->SetVisAttributes(vis);
  if(sd)
    lv->SetSensitiveDetector(sd);
  new G4PVPlacement(make_axis_rotation(frame, axis_angle),
                    frame.InternalToS2SPoint(native_center),
                    lv, name + "PV", world_lv, false, copy_no, check_overlaps);
}

void place_tubs(const G4String& name,
                G4double rmin,
                G4double rmax,
                G4double half_z,
                G4double start_phi,
                G4double delta_phi,
                const G4ThreeVector& native_center,
                const G4RotationMatrix& native_rotation,
                const K18BeamlineFrame& frame,
                G4LogicalVolume* world_lv,
                G4bool check_overlaps,
                G4VisAttributes* vis,
                G4Material* material)
{
  auto* solid = new G4Tubs(name + "Solid", rmin, rmax, half_z,
                           start_phi, delta_phi);
  auto* lv = new G4LogicalVolume(solid, material, name + "LV");
  lv->SetVisAttributes(vis);
  new G4PVPlacement(make_rotation_from_internal(frame, native_rotation),
                    frame.InternalToS2SPoint(native_center),
                    lv, name + "PV", world_lv, false, 0, check_overlaps);
}

void place_window(const G4String& name,
                  G4double radius,
                  G4double thickness,
                  const G4ThreeVector& native_center,
                  G4double beam_axis_angle,
                  const K18BeamlineFrame& frame,
                  G4LogicalVolume* world_lv,
                  G4bool check_overlaps,
                  G4VisAttributes* vis,
                  G4Material* material)
{
  auto* solid = new G4Tubs(name + "Solid", 0., radius, 0.5*thickness,
                           0.*deg, 360.*deg);
  auto* lv = new G4LogicalVolume(solid, material, name + "LV");
  lv->SetVisAttributes(vis);
  new G4PVPlacement(make_axis_rotation(frame, beam_axis_angle),
                    frame.InternalToS2SPoint(native_center),
                    lv, name + "PV", world_lv, false, 0, check_overlaps);
}

void place_solid_on_axis(const G4String& pv_name,
                         G4LogicalVolume* lv,
                         const G4ThreeVector& native_center,
                         G4double axis_angle,
                         const K18BeamlineFrame& frame,
                         G4LogicalVolume* world_lv,
                         G4bool check_overlaps,
                         G4int copy_no = 0,
                         G4double roll = 0.,
                         const G4ThreeVector& local_offset = G4ThreeVector())
{
  new G4PVPlacement(
    make_axis_rotation(frame, axis_angle, roll),
    frame.InternalToS2SPoint(native_center
                             + magnet_offset(axis_angle, local_offset)),
    lv, pv_name, world_lv, false, copy_no, check_overlaps);
}

// Q R-bend: G4Tubs quarter in the (u, z) plane, axis along magnet radial.
// Fusion _bend_tangent_chain_local: r = min(bend_reach, |arc_z - leg_z|).
void place_q_rbend(const G4String& pv_name,
                   G4LogicalVolume* lv,
                   G4double radius,
                   G4double layer_offset,
                   G4double r_arc,
                   G4double axis_deg,
                   G4double side,
                   G4double s_z,
                   G4double z_arc,
                   G4double roll,
                   const G4ThreeVector& native_center,
                   G4double axis_angle,
                   const K18BeamlineFrame& frame,
                   G4LogicalVolume* world_lv,
                   G4bool check_overlaps,
                   G4int copy_no)
{
  const G4double a = axis_deg + roll;
  const G4double ca = std::cos(a);
  const G4double sa = std::sin(a);
  const G4double u_c = side*(layer_offset + r_arc);
  const G4double z_c = z_arc - s_z*r_arc;
  const G4ThreeVector local(radius*ca + u_c*(-sa),
                            radius*sa + u_c*ca,
                            z_c);
  const auto ex = K18BeamlineFrame::XAxis(axis_angle);
  const auto ey = K18BeamlineFrame::YAxis();
  const auto ez = K18BeamlineFrame::Axis(axis_angle);
  const auto er = ca*ex + sa*ey;
  const auto et = -sa*ex + ca*ey;
  G4RotationMatrix native;
  native.rotateAxes(et.unit(), ez.unit(), er.unit());
  new G4PVPlacement(
    make_rotation_from_internal(frame, native),
    frame.InternalToS2SPoint(native_center
                             + magnet_offset(axis_angle, local)),
    lv, pv_name, world_lv, false, copy_no, check_overlaps);
}

G4double d4_bend_phi()
{
  return K18BeamlineFrame::BendAngle()*deg;
}

// Clean-assembly Y-mirror: Fusion phi f maps to G4 phi = 64° - f, so the
// 23.5° face is upstream (near phi=0) and the radial face is downstream.
G4double d4_from_fusion(G4double fusion_phi)
{
  return d4_bend_phi() - fusion_phi;
}

G4Tubs* make_phi_tubs(const G4String& name,
                      G4double rmin, G4double rmax, G4double half_z,
                      G4double phi0, G4double phi1)
{
  G4double dphi = phi1 - phi0;
  if(dphi < 0.){
    std::swap(phi0, phi1);
    dphi = -dphi;
  }
  return new G4Tubs(name, rmin, rmax, half_z, phi0, dphi);
}

void place_d4_tubs_at(const G4String& pv_name,
                      G4VSolid* solid,
                      const G4ThreeVector& native_center,
                      const G4RotationMatrix& native_rotation,
                      const K18BeamlineFrame& frame,
                      G4LogicalVolume* world_lv,
                      G4bool check_overlaps,
                      G4VisAttributes* vis,
                      G4Material* material,
                      G4int copy_no,
                      const G4ThreeVector& local_offset);

G4VSolid* make_d4_exit_hull(const G4String& name,
                            G4double z0, G4double z1,
                            G4double r_near, G4double r_far,
                            G4double half_z, G4double zsgn, G4int n = 48)
{
  std::vector<G4TwoVector> poly;
  poly.reserve(2*n + 4);
  auto add = [&](G4double x, G4double z){
    poly.emplace_back(x, zsgn*z);
  };
  add(0., z0);
  add(0., z1);
  const G4double cin = z1 + r_near;
  for(G4int i=1; i<=n; ++i){
    const G4double a = 270.*deg + (180.*deg - 270.*deg)*i/n;
    add(r_near*std::cos(a), cin + r_near*std::sin(a));
  }
  const G4double coutz = z0 + r_far;
  for(G4int i=0; i<n; ++i){
    const G4double a = 180.*deg + (270.*deg - 180.*deg)*i/n;
    add(r_far*std::cos(a), coutz + r_far*std::sin(a));
  }
  return new G4ExtrudedSolid(name, poly, half_z,
                             G4TwoVector(), 1., G4TwoVector(), 1.);
}

void append_xy_arc(std::vector<G4TwoVector>& poly,
                   G4double cx, G4double cy,
                   G4double x0, G4double y0,
                   G4double x1, G4double y1,
                   G4int n = 24)
{
  const G4double a0 = std::atan2(y0 - cy, x0 - cx);
  G4double da = std::atan2(y1 - cy, x1 - cx) - a0;
  const G4double pi = std::acos(-1.);
  while(da > pi) da -= 2.*pi;
  while(da < -pi) da += 2.*pi;
  const G4double r = std::hypot(x0 - cx, y0 - cy);
  for(G4int i = 1; i <= n; ++i){
    const G4double a = a0 + da*i/n;
    poly.emplace_back(cx + r*std::cos(a), cy + r*std::sin(a));
  }
}

// DXF _1__S right half-plate (x>=0): 450x1300 with the native star opening.
G4VSolid* make_endguard_half(const G4String& name, G4bool right_half)
{
  const G4double sx = right_half ? 1. : -1.;
  std::vector<G4TwoVector> poly;
  poly.reserve(48);
  poly.emplace_back(0., kQYokeHY);
  poly.emplace_back(sx*kQYokeHX, kQYokeHY);
  poly.emplace_back(sx*kQYokeHX, -kQYokeHY);
  poly.emplace_back(0., -kQYokeHY);
  poly.emplace_back(0., -kQGuardStarArm);
  poly.emplace_back(sx*kQGuardStarNeckX, -kQGuardStarArm);
  append_xy_arc(poly, sx*kQGuardStarCx, -kQGuardStarCy,
                sx*kQGuardStarNeckX, -kQGuardStarArm,
                sx*kQGuardStarArm, -kQGuardStarNeckY);
  poly.emplace_back(sx*kQGuardStarArm, kQGuardStarNeckY);
  append_xy_arc(poly, sx*kQGuardStarCx, kQGuardStarCy,
                sx*kQGuardStarArm, kQGuardStarNeckY,
                sx*kQGuardStarNeckX, kQGuardStarArm);
  poly.emplace_back(0., kQGuardStarArm);
  if(!right_half)
    std::reverse(poly.begin(), poly.end());
  return new G4ExtrudedSolid(name, poly, 0.5*kQGuardThick,
                             G4TwoVector(), 1., G4TwoVector(), 1.);
}

std::vector<G4TwoVector> star_opening_polygon(G4double scale, G4int narc = 10)
{
  const G4double arm = scale*kQGuardStarArm;
  const G4double nx = scale*kQGuardStarNeckX;
  const G4double ny = scale*kQGuardStarNeckY;
  const G4double cx = scale*kQGuardStarCx;
  const G4double cy = scale*kQGuardStarCy;
  std::vector<G4TwoVector> poly;
  poly.reserve(8 + 4*narc);
  poly.emplace_back(nx, arm);
  append_xy_arc(poly, cx, cy, nx, arm, arm, ny, narc);
  poly.emplace_back(arm, -ny);
  append_xy_arc(poly, cx, -cy, arm, -ny, nx, -arm, narc);
  poly.emplace_back(-nx, -arm);
  append_xy_arc(poly, -cx, -cy, -nx, -arm, -arm, -ny, narc);
  poly.emplace_back(-arm, ny);
  append_xy_arc(poly, -cx, cy, -arm, ny, -nx, arm, narc);
  return poly;
}

std::vector<G4TwoVector> rect_opening_polygon(G4double hx, G4double hy)
{
  return {
    G4TwoVector(-hx, -hy), G4TwoVector(-hx, hy),
    G4TwoVector(hx, hy), G4TwoVector(hx, -hy)
  };
}

G4ThreeVector on_beam_axis(G4double l, G4double axis_angle)
{
  return (std::abs(axis_angle) < 0.5)
    ? on_straight_axis(l) : on_bent_axis(l);
}

void place_extruded_pipe(const G4String& pv_name,
                         const std::vector<G4TwoVector>& poly,
                         G4double half_z,
                         const G4ThreeVector& native_center,
                         G4double axis_angle,
                         const K18BeamlineFrame& frame,
                         G4LogicalVolume* world_lv,
                         G4bool check_overlaps,
                         G4VisAttributes* vis,
                         G4Material* material,
                         G4int copy_no)
{
  auto* solid = new G4ExtrudedSolid(
    pv_name + "Solid" + std::to_string(copy_no), poly, half_z,
    G4TwoVector(), 1., G4TwoVector(), 1.);
  auto* lv = new G4LogicalVolume(
    solid, material, pv_name + "LV" + std::to_string(copy_no));
  lv->SetVisAttributes(vis);
  place_solid_on_axis(pv_name, lv, native_center, axis_angle, frame, world_lv,
                      check_overlaps, copy_no);
}

void place_pipe_span(const G4String& pv_name,
                     const std::vector<G4TwoVector>& poly,
                     G4double l0, G4double l1, G4double axis_angle,
                     const K18BeamlineFrame& frame,
                     G4LogicalVolume* world_lv,
                     G4bool check_overlaps,
                     G4VisAttributes* vis,
                     G4Material* material,
                     G4int& copy)
{
  if(l1 < l0) std::swap(l0, l1);
  const G4double hz = 0.5*(l1 - l0);
  if(hz < 0.2*mm) return;
  place_extruded_pipe(pv_name, poly, hz, on_beam_axis(0.5*(l0 + l1), axis_angle),
                      axis_angle, frame, world_lv, check_overlaps, vis, material,
                      copy++);
}

void place_pipe_disk(const G4String& pv_name,
                     G4double l, G4double axis_angle,
                     const K18BeamlineFrame& frame,
                     G4LogicalVolume* world_lv,
                     G4bool check_overlaps,
                     G4VisAttributes* hole_vis,
                     G4VisAttributes* flange_vis,
                     G4Material* vacuum,
                     G4Material* flange,
                     G4int& copy)
{
  const auto center = on_beam_axis(l, axis_angle);
  auto* hole = new G4Tubs(pv_name + "Hole" + std::to_string(copy),
                          0., kPipeDiskHoleR, 0.5*kPipeDiskThick,
                          0., 360.*deg);
  auto* hole_lv = new G4LogicalVolume(hole, vacuum,
                                      pv_name + "HoleLV" + std::to_string(copy));
  hole_lv->SetVisAttributes(hole_vis);
  place_solid_on_axis(pv_name, hole_lv, center, axis_angle, frame, world_lv,
                      check_overlaps, copy++);
  auto* ring = new G4Tubs(pv_name + "Flange" + std::to_string(copy),
                          kPipeDiskHoleR, kPipeDiskOuterR, 0.5*kPipeDiskThick,
                          0., 360.*deg);
  auto* ring_lv = new G4LogicalVolume(ring, flange,
                                      pv_name + "FlangeLV" + std::to_string(copy));
  ring_lv->SetVisAttributes(flange_vis);
  place_solid_on_axis("K18PipeFlangePV", ring_lv, center, axis_angle, frame,
                      world_lv, check_overlaps, copy++);
}

G4double d4_fusion_face_phi(G4double radius)
{
  const G4double r0 = kD4YokeR0;
  const G4double r1 = kD4YokeR1;
  const G4double p0x = r0*std::cos(kD4PhiFaceIn);
  const G4double p0y = r0*std::sin(kD4PhiFaceIn);
  const G4double p1x = r1*std::cos(kD4PhiFaceOut);
  const G4double p1y = r1*std::sin(kD4PhiFaceOut);
  const G4double dx = p1x - p0x;
  const G4double dy = p1y - p0y;
  const G4double aa = dx*dx + dy*dy;
  const G4double bb = 2.*(p0x*dx + p0y*dy);
  const G4double cc = p0x*p0x + p0y*p0y - radius*radius;
  const G4double disc = std::max(0., bb*bb - 4.*aa*cc);
  const G4double sdisc = std::sqrt(disc);
  G4double best_d2 = 1.e99;
  G4double best_phi = kD4PhiFaceIn;
  for(G4double t : {(-bb - sdisc)/(2.*aa), (-bb + sdisc)/(2.*aa)}){
    const G4double x = p0x + t*dx;
    const G4double y = p0y + t*dy;
    const G4double d2 = (x - p0x)*(x - p0x) + (y - p0y)*(y - p0y);
    if(d2 < best_d2){
      best_d2 = d2;
      best_phi = std::atan2(y, x);
    }
  }
  return best_phi;
}

// NQ418 block29 / NQ412 block115 NE pole: CORE_POLE_FACE_NE plus 62- or
// 107-family C/D limbs to the ±390 mm window.  Placed as four copies; the
// 45 mm ends use the 107-family solid, the central span uses 62-family.
G4VSolid* make_q_pole_quadrant(const G4String& name, G4double limb,
                              G4double half_z)
{
  static const G4double kPoleFace[][2] = {
    {33.090691, 151.100000}, {33.622333, 148.710678},
    {34.255807, 145.960656}, {34.852762, 143.460654},
    {35.470892, 140.960652}, {36.111343, 138.460649},
    {36.843085, 135.710678}, {37.534529, 133.210678},
    {38.252422, 130.710678}, {39.074516, 127.960636},
    {39.853139, 125.460632}, {40.663423, 122.960628},
    {41.507340, 120.460623}, {42.477030, 117.710678},
    {43.398755, 115.210678}, {44.361369, 112.710678},
    {45.470834, 109.960595}, {46.528687, 107.460586},
    {47.636934, 104.960576}, {48.799262, 102.460565},
    {50.145081, 99.710678}, {51.434679, 97.210678},
    {52.792358, 94.710678}, {54.371169, 91.960502},
    {55.890600, 89.460482}, {57.497396, 86.960459},
    {59.199317, 84.460434}, {61.191513, 81.710678},
    {63.122803, 79.210678}, {65.179974, 76.710678},
    {67.603842, 73.960293}, {69.919895, 71.510405},
    {71.510405, 69.919895}, {73.960293, 67.603842},
    {76.710678, 65.179974}, {79.210678, 63.122803},
    {81.710678, 61.191513}, {84.460434, 59.199317},
    {86.960459, 57.497396}, {89.460482, 55.890600},
    {91.960502, 54.371169}, {94.710678, 52.792358},
    {97.210678, 51.434679}, {99.710678, 50.145081},
    {102.460565, 48.799262}, {104.960576, 47.636934},
    {107.460586, 46.528687}, {109.960595, 45.470834},
    {112.710678, 44.361369}, {115.210678, 43.398755},
    {117.710678, 42.477030}, {120.460623, 41.507340},
    {122.960628, 40.663423}, {125.460632, 39.853139},
    {127.960636, 39.074516}, {130.710678, 38.252422},
    {133.210678, 37.534529}, {135.710678, 36.843085},
    {138.460649, 36.111343}, {140.960652, 35.470892},
    {143.460654, 34.852762}, {145.960656, 34.255807},
    {148.710678, 33.622333}, {151.210678, 33.066448},
  };
  std::vector<G4TwoVector> poly;
  poly.reserve(72);
  for(const auto& p : kPoleFace)
    poly.emplace_back(p[0]*mm, p[1]*mm);
  const G4bool end107 = limb > 80.*mm;
  poly.emplace_back(173.0*mm, end107 ? 54.99345*mm : 35.35*mm);
  poly.emplace_back(kQPoleBreak, limb);
  poly.emplace_back(kQYokeInner, limb);
  poly.emplace_back(kQYokeInner, kQYokeInner);
  poly.emplace_back(limb, kQYokeInner);
  poly.emplace_back(limb, kQPoleBreak);
  poly.emplace_back(end107 ? 54.99345*mm : 35.35*mm, 173.0*mm);
  return new G4ExtrudedSolid(name, poly, half_z,
                             G4TwoVector(), 1., G4TwoVector(), 1.);
}

G4LogicalVolume* iron_lv(const G4String& name, G4VSolid* solid,
                         G4Material* iron, G4VisAttributes* vis)
{
  auto* lv = new G4LogicalVolume(solid, iron, name);
  lv->SetVisAttributes(vis);
  return lv;
}

void place_quad_hardware(const G4String& name,
                         G4double iron_hz,
                         G4double guard_up,
                         G4double guard_dn,
                         const G4ThreeVector& native_center,
                         G4double axis_angle,
                         const K18BeamlineFrame& frame,
                         G4LogicalVolume* world_lv,
                         G4bool check_overlaps,
                         G4VisAttributes* iron_vis,
                         G4VisAttributes* coil_vis,
                         G4VisAttributes* guard_vis,
                         G4Material* iron,
                         G4Material* copper)
{
  // Clean-assembly physical master: 4 clamp plates + 4 pole quadrants.
  // These remain separate because G4 BooleanProcessor cannot form a stable
  // polyhedron for the combined tree (K18Q*YokeU3).
  const G4double hy_tb = 0.5*(kQYokeHY - kQYokeInner);
  const G4double y_tb = 0.5*(kQYokeHY + kQYokeInner);
  const G4double hx_lr = 0.5*(kQYokeHX - kQYokeInner);
  const G4double x_lr = 0.5*(kQYokeHX + kQYokeInner);
  auto* clamp_tb = iron_lv(
    name + "ClampTBLV",
    new G4Box(name + "ClampTB", kQYokeHX, hy_tb, iron_hz), iron, iron_vis);
  auto* clamp_lr = iron_lv(
    name + "ClampLRLV",
    new G4Box(name + "ClampLR", hx_lr, kQYokeInner, iron_hz), iron, iron_vis);
  G4int copy = 0;
  place_solid_on_axis(name + "PV", clamp_tb, native_center, axis_angle,
                      frame, world_lv, check_overlaps, copy++, 0.,
                      G4ThreeVector(0., y_tb, 0.));
  place_solid_on_axis(name + "PV", clamp_tb, native_center, axis_angle,
                      frame, world_lv, check_overlaps, copy++, 0.,
                      G4ThreeVector(0., -y_tb, 0.));
  place_solid_on_axis(name + "PV", clamp_lr, native_center, axis_angle,
                      frame, world_lv, check_overlaps, copy++, 0.,
                      G4ThreeVector(x_lr, 0., 0.));
  place_solid_on_axis(name + "PV", clamp_lr, native_center, axis_angle,
                      frame, world_lv, check_overlaps, copy++, 0.,
                      G4ThreeVector(-x_lr, 0., 0.));

  const G4double central_hz = iron_hz - kQPoleTaperZ;
  const G4double slice_hz = 0.5*kQPoleTaperZ/kQPoleTaperSlices;
  auto* pole62 = iron_lv(
    name + "Pole62LV",
    make_q_pole_quadrant(name + "Pole62", kQPoleLimb62, central_hz),
    iron, iron_vis);
  for(G4int i=0; i<4; ++i){
    const G4double roll = i*90.*deg;
    place_solid_on_axis(name + "PV", pole62, native_center, axis_angle,
                        frame, world_lv, check_overlaps, copy++, roll);
    for(G4int s=0; s<kQPoleTaperSlices; ++s){
      const G4double frac = (s + 0.5)/kQPoleTaperSlices;
      const G4double limb = kQPoleLimb62
        + frac*(kQPoleLimb107 - kQPoleLimb62);
      const G4String slice_name = name + "PoleTaper" + std::to_string(s);
      auto* pole_slice = iron_lv(
        slice_name + "LV",
        make_q_pole_quadrant(slice_name, limb, slice_hz),
        iron, iron_vis);
      const G4double z_pos = central_hz + (s + 0.5)*kQPoleTaperZ/kQPoleTaperSlices;
      place_solid_on_axis(name + "PV", pole_slice, native_center, axis_angle,
                          frame, world_lv, check_overlaps, copy++, roll,
                          G4ThreeVector(0., 0., z_pos));
      place_solid_on_axis(name + "PV", pole_slice, native_center, axis_angle,
                          frame, world_lv, check_overlaps, copy++, roll,
                          G4ThreeVector(0., 0., -z_pos));
    }
  }

  // One filled hull per quadrant.  Z-legs stop at the R-bend tangent
  // (z=z_leg, y=off_in..off_out).  R-bend is a quarter about the 108 mm
  // corner so both joins are coplanar.  Endturn z-span is the R-bend
  // arrival.  0.5 mm overlap at each joint.
  const G4double kJointOl = 0.5*mm;
  const G4double r_in = kQCoilDInner - kQCoilHalf;
  const G4double r_out = kQCoilSlotR0 + kQCoilSlotCount*kQCoilSlotPitch
    + kQCoilHalf;
  const G4double off_in = 9.*mm - kQCoilHalf;
  const G4double off_out = 43.*mm + kQCoilHalf;
  const G4double pack_hr = 0.5*(r_out - r_in);
  const G4double pack_ho = 0.5*(off_out - off_in);
  const G4double r_mid = 0.5*(r_in + r_out);
  const G4double off_mid = 0.5*(off_in + off_out);
  const G4double r_end_in = std::max(r_in - kQCoilCorner, 0.1*mm);
  const G4double r_end_out = r_out - kQCoilCorner;
  const G4double r_arc_in = std::max(kQCoilCorner - off_out, 0.1*mm);
  const G4double r_arc_out = kQCoilCorner - off_in;
  const G4double r_arc_place = kQCoilCorner - off_mid;
  const G4double zend_hi = iron_hz + kQCoilHalf + kQCoilClearZ
    + 2.*kQCoilPitchZ;
  const G4double z_leg = zend_hi + kQCoilHalf - r_arc_out;
  const G4double z_end_c = z_leg + 0.5*(r_arc_in + r_arc_out);
  const G4double z_end_h = 0.5*(r_arc_out - r_arc_in) + kJointOl;
  const G4double z_arc_place = z_leg + r_arc_place;

  auto* coil_legx = new G4LogicalVolume(
    new G4Box(name + "CoilHullLegX", pack_hr, pack_ho, z_leg + kJointOl),
    copper, name + "CoilHullLegXLV");
  auto* coil_legy = new G4LogicalVolume(
    new G4Box(name + "CoilHullLegY", pack_ho, pack_hr, z_leg + kJointOl),
    copper, name + "CoilHullLegYLV");
  auto* coil_end = new G4LogicalVolume(
    new G4Tubs(name + "CoilHullEnd", r_end_in, r_end_out, z_end_h,
               0.*deg, 90.*deg),
    copper, name + "CoilHullEndLV");
  auto* coil_rb_sf = new G4LogicalVolume(
    new G4Tubs(name + "CoilHullRbSF", r_arc_in, r_arc_out, pack_hr,
               90.*deg, 90.*deg),
    copper, name + "CoilHullRbSFLV");
  auto* coil_rb_ef = new G4LogicalVolume(
    new G4Tubs(name + "CoilHullRbEF", r_arc_in, r_arc_out, pack_hr,
               0.*deg, 90.*deg),
    copper, name + "CoilHullRbEFLV");
  auto* coil_rb_eb = new G4LogicalVolume(
    new G4Tubs(name + "CoilHullRbEB", r_arc_in, r_arc_out, pack_hr,
               -90.*deg, 90.*deg),
    copper, name + "CoilHullRbEBLV");
  auto* coil_rb_sb = new G4LogicalVolume(
    new G4Tubs(name + "CoilHullRbSB", r_arc_in, r_arc_out, pack_hr,
               180.*deg, 90.*deg),
    copper, name + "CoilHullRbSBLV");
  for(auto* lv : {coil_legx, coil_legy, coil_end,
                  coil_rb_sf, coil_rb_ef, coil_rb_eb, coil_rb_sb})
    lv->SetVisAttributes(coil_vis);

  G4int ccopy = 0;
  for(G4int i=0; i<4; ++i){
    const G4double roll = i*90.*deg;
    const G4double c = std::cos(roll);
    const G4double s = std::sin(roll);
    auto rz = [c, s](G4double x, G4double y, G4double z){
      return G4ThreeVector(c*x - s*y, s*x + c*y, z);
    };
    place_solid_on_axis(name + "CoilPV", coil_legx, native_center, axis_angle,
                        frame, world_lv, check_overlaps, ccopy++, roll,
                        rz(r_mid, off_mid, 0.));
    place_solid_on_axis(name + "CoilPV", coil_legy, native_center, axis_angle,
                        frame, world_lv, check_overlaps, ccopy++, roll,
                        rz(off_mid, r_mid, 0.));
    place_solid_on_axis(name + "CoilPV", coil_end, native_center, axis_angle,
                        frame, world_lv, check_overlaps, ccopy++, roll,
                        rz(kQCoilCorner, kQCoilCorner, z_end_c));
    place_solid_on_axis(name + "CoilPV", coil_end, native_center, axis_angle,
                        frame, world_lv, check_overlaps, ccopy++, roll,
                        rz(kQCoilCorner, kQCoilCorner, -z_end_c));
    place_q_rbend(name + "CoilPV", coil_rb_sf, r_mid, off_mid, r_arc_place,
                  0., +1., +1., z_arc_place, roll, native_center, axis_angle,
                  frame, world_lv, check_overlaps, ccopy++);
    place_q_rbend(name + "CoilPV", coil_rb_ef, r_mid, off_mid, r_arc_place,
                  90.*deg, -1., +1., z_arc_place, roll, native_center,
                  axis_angle, frame, world_lv, check_overlaps, ccopy++);
    place_q_rbend(name + "CoilPV", coil_rb_eb, r_mid, off_mid, r_arc_place,
                  90.*deg, -1., -1., -z_arc_place, roll, native_center,
                  axis_angle, frame, world_lv, check_overlaps, ccopy++);
    place_q_rbend(name + "CoilPV", coil_rb_sb, r_mid, off_mid, r_arc_place,
                  0., +1., -1., -z_arc_place, roll, native_center, axis_angle,
                  frame, world_lv, check_overlaps, ccopy++);
  }

  // _1__S 900x1300 half plates with the native 4-neck / 4-arc star opening.
  auto* guard_r = iron_lv(
    name + "GuardRightLV", make_endguard_half(name + "GuardRight", true),
    iron, guard_vis);
  auto* guard_l = iron_lv(
    name + "GuardLeftLV", make_endguard_half(name + "GuardLeft", false),
    iron, guard_vis);
  const G4double z_up = guard_up - 0.5*kQGuardThick;
  const G4double z_dn = guard_dn + 0.5*kQGuardThick;
  G4int gcopy = 0;
  for(G4double z : {z_up, z_dn}){
    place_solid_on_axis(name + "GuardPV", guard_r, native_center, axis_angle,
                        frame, world_lv, check_overlaps, gcopy++, 0.,
                        G4ThreeVector(0., 0., z));
    place_solid_on_axis(name + "GuardPV", guard_l, native_center, axis_angle,
                        frame, world_lv, check_overlaps, gcopy++, 0.,
                        G4ThreeVector(0., 0., z));
  }
}

void place_d4_tubs_at(const G4String& pv_name,
                      G4VSolid* solid,
                      const G4ThreeVector& native_center,
                      const G4RotationMatrix& native_rotation,
                      const K18BeamlineFrame& frame,
                      G4LogicalVolume* world_lv,
                      G4bool check_overlaps,
                      G4VisAttributes* vis,
                      G4Material* material,
                      G4int copy_no,
                      const G4ThreeVector& local_offset)
{
  auto* lv = new G4LogicalVolume(solid, material, pv_name + "LV" + std::to_string(copy_no));
  lv->SetVisAttributes(vis);
  const auto native = native_center + native_rotation*local_offset;
  new G4PVPlacement(make_rotation_from_internal(frame, native_rotation),
                    frame.InternalToS2SPoint(native),
                    lv, pv_name, world_lv, false, copy_no, check_overlaps);
}

void place_d4_tubs_at(const G4String& pv_name,
                      G4VSolid* solid,
                      const G4ThreeVector& native_center,
                      const G4RotationMatrix& native_rotation,
                      const K18BeamlineFrame& frame,
                      G4LogicalVolume* world_lv,
                      G4bool check_overlaps,
                      G4VisAttributes* vis,
                      G4Material* material,
                      G4int copy_no)
{
  place_d4_tubs_at(pv_name, solid, native_center, native_rotation, frame,
                   world_lv, check_overlaps, vis, material, copy_no,
                   G4ThreeVector());
}

// 90-deg riser envelope in the face-normal / z plane.  Fusion sketches this
// on a constant-radius plane and extrudes in r.  G4Tubs axis = radial;
// solid X = inside (decreasing G4 phi), Y = +z, Z = -radial.
void place_d4_riser(const G4String& pv_name,
                    G4double r0,
                    G4double r1,
                    G4double phi_face,
                    G4double cx,
                    G4double cz,
                    G4double r_path,
                    G4double half_thick,
                    const G4ThreeVector& native_center,
                    const G4RotationMatrix& native_rotation,
                    const K18BeamlineFrame& frame,
                    G4LogicalVolume* world_lv,
                    G4bool check_overlaps,
                    G4VisAttributes* vis,
                    G4Material* copper,
                    G4int copy_no)
{
  const G4double rmid = 0.5*(r0 + r1);
  const G4double hz = 0.5*(r1 - r0);
  const G4double c = std::cos(phi_face);
  const G4double s = std::sin(phi_face);
  const G4ThreeVector radial(c, s, 0.);
  const G4ThreeVector inside(s, -c, 0.); // -phi-hat
  const G4ThreeVector zhat(0., 0., 1.);
  G4RotationMatrix riser_in_tubs;
  riser_in_tubs.rotateAxes(inside.unit(), zhat.unit(), (-radial).unit());
  const G4RotationMatrix tot = native_rotation*riser_in_tubs;
  const G4ThreeVector c_tubs = rmid*radial + cx*inside + cz*zhat;
  const G4ThreeVector native = native_center + native_rotation*c_tubs;
  const G4double phi0 = (cz > 0.) ? 180.*deg : 90.*deg;
  auto* solid = new G4Tubs(pv_name + "Solid" + std::to_string(copy_no),
                           r_path - half_thick, r_path + half_thick, hz,
                           phi0, 90.*deg);
  auto* lv = new G4LogicalVolume(solid, copper,
                                 pv_name + "LV" + std::to_string(copy_no));
  lv->SetVisAttributes(vis);
  new G4PVPlacement(make_rotation_from_internal(frame, tot),
                    frame.InternalToS2SPoint(native),
                    lv, pv_name, world_lv, false, copy_no, check_overlaps);
}

void place_d4_native(const G4String& pv_name,
                     G4VSolid* solid,
                     const G4ThreeVector& c_tubs,
                     const G4RotationMatrix& rot_in_tubs,
                     const G4ThreeVector& native_center,
                     const G4RotationMatrix& native_rotation,
                     const K18BeamlineFrame& frame,
                     G4LogicalVolume* world_lv,
                     G4bool check_overlaps,
                     G4VisAttributes* vis,
                     G4Material* material,
                     G4int copy_no)
{
  auto* lv = new G4LogicalVolume(solid, material,
                                 pv_name + "LV" + std::to_string(copy_no));
  lv->SetVisAttributes(vis);
  const G4RotationMatrix tot = native_rotation*rot_in_tubs;
  const G4ThreeVector native = native_center + native_rotation*c_tubs;
  new G4PVPlacement(make_rotation_from_internal(frame, tot),
                    frame.InternalToS2SPoint(native),
                    lv, pv_name, world_lv, false, copy_no, check_overlaps);
}

void place_d4_hardware(const G4ThreeVector& native_center,
                       const G4RotationMatrix& native_rotation,
                       const K18BeamlineFrame& frame,
                       G4LogicalVolume* world_lv,
                       G4bool check_overlaps,
                       G4VisAttributes* iron_vis,
                       G4VisAttributes* coil_vis,
                       G4Material* iron,
                       G4Material* copper)
{
  // v32 H-frame: inner/outer uprights plus top/bottom lids.  Each radial
  // band uses face_angle_at_radius() so the 23.5° shell is a staircase of
  // G4Tubs instead of one inner-phi / one outer-phi step at R4000.
  const G4double phi_down = d4_from_fusion(0.);
  const G4double lid_hz = 0.5*(kD4YokeHZ - kD4SlotHZ);
  const G4double lid_z = kD4SlotHZ + lid_hz;
  const G4double dr = (kD4YokeR1 - kD4YokeR0)/kD4YokeRadialSlices;
  G4int copy = 0;
  for(G4int ib=0; ib<kD4YokeRadialSlices; ++ib){
    const G4double r0 = kD4YokeR0 + ib*dr;
    const G4double r1 = r0 + dr;
    const G4double rmid = 0.5*(r0 + r1);
    const G4double phi_up = d4_from_fusion(d4_fusion_face_phi(rmid));
    const G4bool in_slot = (r0 + 0.5*mm >= kD4SlotR0 && r1 - 0.5*mm <= kD4SlotR1);
    if(in_slot){
      auto* lid = make_phi_tubs("K18D4Lid" + std::to_string(copy),
                                r0, r1, lid_hz, phi_up, phi_down);
      for(G4double z : {lid_z, -lid_z}){
        place_d4_tubs_at("K18D4PV", lid, native_center, native_rotation,
                         frame, world_lv, check_overlaps, iron_vis, iron, copy++,
                         G4ThreeVector(0., 0., z));
      }
    }else{
      place_d4_tubs_at("K18D4PV",
        make_phi_tubs("K18D4Yoke" + std::to_string(copy),
                      r0, r1, kD4YokeHZ, phi_up, phi_down),
        native_center, native_rotation, frame, world_lv, check_overlaps,
        iron_vis, iron, copy++);
    }
  }

  // Fusion pole_end_taper_cutter removes iron at the ends (122 -> 83 mm
  // ruled chamfer).  The central |z|<83 region remains the beam aperture.

  // Fusion turn_grid hull: one inner pack and one outer pack.  Pancake and
  // end-bends follow the pack-wrap outline (not the 4-layer turn CLs).
  // The radial ArchFace is limited to the winding pack and leaves the bore.
  G4int ccopy = 0;
  const G4ThreeVector zhat(0., 0., 1.);
  const G4double c_dn = std::cos(phi_down);
  const G4double s_dn = std::sin(phi_down);
  const G4ThreeVector ent_radial(c_dn, s_dn, 0.);
  const G4ThreeVector ent_inside(s_dn, -c_dn, 0.); // -phi-hat
  const G4ThreeVector ent_n_out = -ent_inside;
  const G4ThreeVector ent_origin = kRhoD4*ent_radial;

  const G4double phi_fi = d4_from_fusion(kD4PhiFaceIn);
  const G4double phi_fo = d4_from_fusion(kD4PhiFaceOut);
  const G4ThreeVector p_fi(kD4YokeR0*std::cos(phi_fi),
                           kD4YokeR0*std::sin(phi_fi), 0.);
  const G4ThreeVector p_fo(kD4YokeR1*std::cos(phi_fo),
                           kD4YokeR1*std::sin(phi_fo), 0.);
  G4ThreeVector s_hat = (p_fo - p_fi).unit();
  G4ThreeVector n_out(-s_hat.y(), s_hat.x(), 0.);
  const G4double phi_beam = d4_from_fusion(d4_fusion_face_phi(kRhoD4));
  const G4ThreeVector origin(kRhoD4*std::cos(phi_beam),
                             kRhoD4*std::sin(phi_beam), 0.);
  const G4ThreeVector into_magnet(-origin.y(), origin.x(), 0.);
  if(n_out.dot(into_magnet) > 0.) n_out = -n_out;
  const G4ThreeVector inside = -n_out;
  const G4ThreeVector tubs_z = inside.cross(zhat).unit();

  auto face_s_at = [&](G4double radius){
    const G4double b = 2.*origin.dot(s_hat);
    const G4double c = origin.mag2() - radius*radius;
    const G4double disc = std::max(0., b*b - 4.*c);
    const G4double sdisc = std::sqrt(disc);
    G4double best = 0.;
    G4double best_pen = 1.e99;
    for(G4double s : {0.5*(-b - sdisc), 0.5*(-b + sdisc)}){
      const G4ThreeVector p = origin + s*s_hat;
      const G4double t = (p - p_fi).dot(s_hat);
      const G4double pen = (t < 0. ? -t : 0.)
        + (t > (p_fo - p_fi).mag() ? t - (p_fo - p_fi).mag() : 0.);
      if(pen < best_pen){
        best_pen = pen;
        best = s;
      }
    }
    return best;
  };

  auto place_exit_envelope = [&](G4double s_c, G4double pack_hs,
                                 G4double zsgn){
    G4RotationMatrix rin;
    rin.rotateAxes(inside.unit(), zhat.unit(), tubs_z);
    const G4ThreeVector c_tubs = origin + s_c*s_hat;
    auto* solid = make_d4_exit_hull(
      "K18D4ExitEnv" + std::to_string(ccopy),
      kD4LegZ0, kD4LegZ1, kD4ExitYNear, kD4ExitYFar, pack_hs, zsgn);
    place_d4_native("K18D4CoilPV", solid, c_tubs, rin,
                    native_center, native_rotation, frame, world_lv,
                    check_overlaps, coil_vis, copper, ccopy++);
  };

  auto place_ent_envelope = [&](G4double r0p, G4double r1p, G4double zsgn){
    // Concentric pack wrap about Fusion (RISER_CX, RISER_CZ).  Same solid
    // as place_d4_riser: join is the pancake face y=75, z=1..83.
    const G4double r_path = 0.5*(kD4RiserRIn + kD4RiserROut);
    const G4double half_th = 0.5*(kD4RiserROut - kD4RiserRIn);
    place_d4_riser("K18D4CoilPV", r0p, r1p, phi_down,
                   kD4RiserCX, zsgn*kD4RiserCZ, r_path, half_th,
                   native_center, native_rotation, frame, world_lv,
                   check_overlaps, coil_vis, copper, ccopy++);
  };

  auto place_face_box = [&](const G4ThreeVector& org,
                            const G4ThreeVector& sh,
                            const G4ThreeVector& no,
                            G4double s_c, G4double y_out, G4double z_c,
                            G4double hs, G4double hy, G4double hz){
    G4RotationMatrix rin;
    rin.rotateAxes(sh.unit(), no.unit(), sh.cross(no).unit());
    const G4ThreeVector c_tubs = org + s_c*sh + y_out*no + z_c*zhat;
    auto* solid = new G4Box("K18D4FaceBox" + std::to_string(ccopy),
                            hs, hy, hz);
    place_d4_native("K18D4CoilPV", solid, c_tubs, rin,
                    native_center, native_rotation, frame, world_lv,
                    check_overlaps, coil_vis, copper, ccopy++);
  };

  auto place_shoulder = [&](const G4ThreeVector& org,
                            const G4ThreeVector& sh,
                            const G4ThreeVector& no,
                            G4double s_center, G4double y_out, G4double zsh,
                            G4double r_sh, G4double half_w, G4double half_h,
                            G4bool inner_side, G4double zsgn){
    G4RotationMatrix rin;
    rin.rotateAxes(sh.unit(), zhat.unit(), sh.cross(zhat).unit());
    const G4ThreeVector c_tubs =
      org + s_center*sh + y_out*no + zsgn*zsh*zhat;
    const G4double phi0 = inner_side
      ? ((zsgn > 0.) ? 90.*deg : 180.*deg)
      : ((zsgn > 0.) ? 0.*deg : 270.*deg);
    const G4double rmin = std::max(r_sh - half_w, 0.1*mm);
    auto* solid = new G4Tubs("K18D4ArchSh" + std::to_string(ccopy),
                             rmin, r_sh + half_w, half_h, phi0, 90.*deg);
    place_d4_native("K18D4CoilPV", solid, c_tubs, rin,
                    native_center, native_rotation, frame, world_lv,
                    check_overlaps, coil_vis, copper, ccopy++);
  };

  auto place_staple = [&](const G4ThreeVector& org,
                          const G4ThreeVector& sh,
                          const G4ThreeVector& no,
                          G4double s_i, G4double s_o, G4double y_out,
                          G4double z_bot, G4double zsh,
                          G4double half_w, G4double half_h, G4double zsgn){
    if(zsh > z_bot + 0.5*mm){
      const G4double hz = 0.5*(zsh - z_bot);
      const G4double zc = 0.5*(zsh + z_bot);
      place_face_box(org, sh, no, s_i, y_out, zsgn*zc, half_w, half_h, hz);
      place_face_box(org, sh, no, s_o, y_out, zsgn*zc, half_w, half_h, hz);
    }
    const G4double r_in = std::abs(s_i) - kD4ArchCS;
    const G4double r_out = std::abs(s_o) - kD4ArchCS;
    if(r_in <= 2.*half_w || r_out <= 2.*half_w)
      return;
    place_shoulder(org, sh, no, -kD4ArchCS, y_out, zsh, r_in,
                   half_w, half_h, true, zsgn);
    place_shoulder(org, sh, no, +kD4ArchCS, y_out, zsh, r_out,
                   half_w, half_h, false, zsgn);
    const G4double z_bar = zsh + 0.5*(r_in + r_out);
    place_face_box(org, sh, no, 0., y_out, zsgn*z_bar,
                   kD4ArchCS, half_h, half_w);
  };

  // Main+Aux radial hull.  Pancake legs + pack-wrap bends from measured
  // outline (entrance concentric r=100-182; exit compact R=24/120).
  G4double in0 = 1.e99;
  G4double in1 = -1.e99;
  G4double out0 = 1.e99;
  G4double out1 = -1.e99;
  auto accum_pack = [&](const G4double* roffs, G4int nroff,
                        G4double w, G4double c_in, G4double c_out){
    for(G4int i=0; i<nroff; ++i){
      in0 = std::min(in0, c_in + roffs[i] - 0.5*w);
      in1 = std::max(in1, c_in + roffs[i] + 0.5*w);
      out0 = std::min(out0, c_out - roffs[i] - 0.5*w);
      out1 = std::max(out1, c_out - roffs[i] + 0.5*w);
    }
  };
  accum_pack(kD4MainRoff, 4, kD4MainW, kD4MainCin, kD4MainCout);
  accum_pack(kD4AuxRoff, 4, kD4AuxW, kD4AuxCin, kD4AuxCout);

  const G4double leg_hz = 0.5*(kD4LegZ1 - kD4LegZ0);
  const G4double leg_z = kD4LegZcen;
  auto place_pack_legs = [&](G4double r0, G4double r1){
    const G4double drp = (r1 - r0)/kD4CoilPackSlices;
    for(G4int is=0; is<kD4CoilPackSlices; ++is){
      const G4double ra = r0 + is*drp;
      const G4double rb = ra + drp;
      const G4double rmid = 0.5*(ra + rb);
      const G4double phi_up = d4_from_fusion(d4_fusion_face_phi(rmid));
      const G4double phi_dn = phi_down
        - (kD4RiserCX + kD4EntLegOverlap)/rmid;
      for(G4double z : {leg_z, -leg_z}){
        place_d4_tubs_at("K18D4CoilPV",
          make_phi_tubs("K18D4Pack" + std::to_string(ccopy),
                        ra, rb, leg_hz, phi_up, phi_dn),
          native_center, native_rotation, frame, world_lv, check_overlaps,
          coil_vis, copper, ccopy++, G4ThreeVector(0., 0., z));
      }
    }
  };
  place_pack_legs(in0, in1);
  place_pack_legs(out0, out1);

  const G4double s_in = face_s_at(0.5*(in0 + in1));
  const G4double s_out = face_s_at(0.5*(out0 + out1));
  const G4double hs_in = 0.5*std::abs(face_s_at(in1) - face_s_at(in0));
  const G4double hs_out = 0.5*std::abs(face_s_at(out1) - face_s_at(out0));
  const G4double s_ent_in = 0.5*(in0 + in1) - kRhoD4;
  const G4double s_ent_out = 0.5*(out0 + out1) - kRhoD4;
  const G4double hw_ent_in = 0.5*(in1 - in0);

  const G4double ent_r_path = 0.5*(kD4RiserRIn + kD4RiserROut);
  const G4double y_ent = ent_r_path - kD4RiserCX; // 66 mm outside of U-turn
  const G4double hy_ent = 0.5*kD4ArchAx;         // 41 mm = half of 82 mm wrap
  const G4double y_ex = 0.5*(kD4ExitYNear + kD4ExitYFar); // 72 mm
  const G4double hy_ex = 0.5*(kD4ExitYFar - kD4ExitYNear); // 48 mm
  const G4double z_ex_top = kD4LegZ0 + kD4ExitYFar;

  for(G4double zsgn : {1., -1.}){
    place_ent_envelope(in0, in1, zsgn);
    place_ent_envelope(out0, out1, zsgn);
    place_staple(ent_origin, ent_radial, ent_n_out,
                 s_ent_in, s_ent_out, y_ent, kD4RiserCZ, kD4ArchCZ,
                 hw_ent_in, hy_ent, zsgn);
    place_exit_envelope(s_in, hs_in, zsgn);
    place_exit_envelope(s_out, hs_out, zsgn);
    place_staple(origin, s_hat, n_out,
                 s_in, s_out, y_ex, z_ex_top, kD4ArchCZ,
                 hs_in, hy_ex, zsgn);
  }
  G4cout << "[K18D4] hull entrance: center=("
         << kD4RiserCX << "," << kD4RiserCZ << ") r="
         << kD4RiserRIn << ".." << kD4RiserROut
         << " join y=" << kD4RiserCX
         << " z=" << (kD4RiserCZ - kD4RiserRIn)
         << ".." << (kD4RiserCZ - kD4RiserROut)
         << "  exit join y=0 z=" << kD4LegZ0 << ".." << kD4LegZ1
         << " R=" << kD4ExitYNear << "/" << kD4ExitYFar
         << "  copies=" << ccopy << G4endl;
  G4cout << "[K18D4] turn CLs (not hull) CX={";
  for(G4int i=0; i<kD4CoilLayers; ++i){
    if(i) G4cout << ",";
    G4cout << kD4BendLayerCX[i];
  }
  G4cout << "} R={";
  for(G4int i=0; i<kD4CoilLayers; ++i){
    if(i) G4cout << ",";
    G4cout << kD4BendLayerR[i];
  }
  G4cout << "}" << G4endl;
}


void place_k18_chamber(const G4String& chamber_name,
                       const std::vector<K18PlaneSpec>& planes,
                       const std::vector<K18PassivePlaneSpec>& passive_planes,
                       K18Reference reference,
                       G4double half_x,
                       G4double half_y,
                       G4double half_z,
                       const K18BeamlineFrame& frame,
                       G4LogicalVolume* world_lv,
                       G4bool check_overlaps,
                       G4VisAttributes* chamber_vis,
                       G4VisAttributes* plane_vis,
                       G4Material* chamber_material,
                       G4Material* plane_material,
                       DCSD* sd,
                       G4double longitudinal_shift = 0.)
{
  if(planes.empty())
    return;

  G4double min_l = geomMan.GetLocalZ(planes.front().name)*mm;
  G4double max_l = min_l;
  for(const auto& plane : planes){
    const G4double l = geomMan.GetLocalZ(plane.name)*mm;
    min_l = std::min(min_l, l);
    max_l = std::max(max_l, l);
  }

  const G4double center_l = 0.5*(min_l + max_l);
  const G4double actual_half_z = std::max(half_z, 0.5*(max_l - min_l) + 1.*mm);

  for(const auto& plane : passive_planes){
    if(!plane.material)
      throw std::runtime_error("K18 passive plane has no material: " +
                               std::string(plane.name));
    if(plane.half_x > half_x || plane.half_y > half_y ||
       std::abs(plane.local_l - center_l) + plane.half_z > actual_half_z)
      throw std::runtime_error("K18 passive plane is outside chamber mother: " +
                               std::string(plane.name));
  }

  auto* chamber_solid = new G4Box(chamber_name + "Solid",
                                  half_x, half_y, actual_half_z);
  auto* chamber_lv = new G4LogicalVolume(chamber_solid, chamber_material,
                                         chamber_name + "LV");
  chamber_lv->SetVisAttributes(chamber_vis);
  new G4PVPlacement(make_axis_rotation(frame, axis_angle(reference)),
                    frame.InternalToS2SPoint(native_from_reference(reference,
                                                                   center_l + longitudinal_shift)),
                    chamber_lv, chamber_name + "PV", world_lv, false, 0,
                    check_overlaps);

  for(const auto& plane : planes){
    const G4double l = geomMan.GetLocalZ(plane.name)*mm;
    auto* plane_solid = new G4Box(G4String(plane.name) + "Solid",
                                  plane.half_x, plane.half_y, plane.half_z);
    auto* plane_lv = new G4LogicalVolume(plane_solid, plane_material,
                                         G4String(plane.name) + "LV");
    plane_lv->SetVisAttributes(plane_vis);
    if(sd)
      plane_lv->SetSensitiveDetector(sd);
    new G4PVPlacement(make_local_stereo_rotation(geomMan.GetTiltAngle(plane.name)),
                      G4ThreeVector(0., 0., l - center_l),
                      plane_lv, G4String(plane.name) + "PV", chamber_lv,
                      false, geomMan.GetDetectorId(plane.name), check_overlaps);
  }

  for(const auto& plane : passive_planes){
    auto* solid = new G4Box(plane.name + "Solid",
                            plane.half_x, plane.half_y, plane.half_z);
    auto* lv = new G4LogicalVolume(solid, plane.material, plane.name + "LV");
    lv->SetVisAttributes(plane_vis);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0., 0.,
                                    plane.local_l - center_l),
                      lv, plane.name + "PV", chamber_lv, false, 0,
                      check_overlaps);
  }
}

std::vector<K18PassivePlaneSpec>
make_bcout_films(const G4String& chamber_name,
                 const std::vector<std::pair<G4String, G4String>>& layer_pairs,
                 G4double half_x,
                 G4double half_y,
                 G4double window_thickness,
                 G4double cathode_thickness,
                 G4double anode_cathode_gap,
                 G4double aluminum_coating_thickness,
                 G4double carbon_coating_thickness,
                 G4Material* mylar,
                 G4Material* aluminum,
                 G4Material* carbon)
{
  std::vector<K18PassivePlaneSpec> planes;
  if(layer_pairs.empty() || window_thickness <= 0. ||
     cathode_thickness <= 0. || anode_cathode_gap <= 0. || !mylar)
    return planes;

  const G4double first_l = geomMan.GetLocalZ(layer_pairs.front().first)*mm;
  const G4double last_l = geomMan.GetLocalZ(layer_pairs.back().second)*mm;
  const G4double chamber_center_l = 0.5*(first_l + last_l);
  // MWDC-629 assembly drawing: full beam-axis envelope is 85 mm.
  const G4double window_offset = 42.5*mm;
  for(G4int side=0; side<2; ++side){
    const G4double sign = side == 0 ? -1. : 1.;
    const G4double window_l = chamber_center_l + sign*window_offset;
    planes.push_back({chamber_name + "Window" + std::to_string(side),
                      window_l, half_x, half_y,
                      0.5*window_thickness, mylar});
    if(aluminum_coating_thickness > 0. && aluminum){
      for(G4int face=0; face<2; ++face){
        const G4double face_sign = face == 0 ? -1. : 1.;
        planes.push_back({chamber_name + "WindowAl" +
                            std::to_string(side) + "Face" +
                            std::to_string(face),
                          window_l + face_sign*0.5*(window_thickness +
                                                    aluminum_coating_thickness),
                          half_x, half_y,
                          0.5*aluminum_coating_thickness, aluminum});
      }
    }
  }

  // Each X/X', U/U', or V/V' doublet is bounded by two 12-um
  // carbon-coated Mylar cathodes.
  G4int index = 0;
  for(const auto& pair : layer_pairs){
    const G4double first = geomMan.GetLocalZ(pair.first)*mm;
    const G4double second = geomMan.GetLocalZ(pair.second)*mm;
    // The drawing places the planar cathodes 2 mm from the adjacent anode
    // wire plane.  The u/v stereo angle belongs to the wires, not the face.
    for(const auto cathode_l : {first - anode_cathode_gap,
                                second + anode_cathode_gap}){
      const auto cathode_name = chamber_name + "Cathode" +
        std::to_string(index++);
      planes.push_back({cathode_name, cathode_l, half_x, half_y,
                        0.5*cathode_thickness, mylar});
      if(carbon_coating_thickness > 0. && carbon){
        for(G4int face=0; face<2; ++face){
          const G4double face_sign = face == 0 ? -1. : 1.;
          planes.push_back({cathode_name + "CarbonFace" +
                              std::to_string(face),
                            cathode_l + face_sign*0.5*(cathode_thickness +
                                                       carbon_coating_thickness),
                            half_x, half_y,
                            0.5*carbon_coating_thickness, carbon});
        }
      }
    }
  }
  return planes;
}

void place_k18_bft_realistic(const K18BeamlineFrame& frame,
                             G4LogicalVolume* world_lv,
                             G4bool check_overlaps,
                             G4VisAttributes* plane_vis,
                             G4VisAttributes* window_vis,
                             G4Material* fibre_material,
                             G4Material* proxy_material,
                             G4Material* black_sheet_material,
                             DCSD* sd,
                             G4double fibre_diameter,
                             G4double black_sheet_thickness)
{
  const std::vector<G4String> names{"BFT-X", "BFT-XP"};
  const G4double half_x = 80.*mm;
  const G4double half_y = 40.*mm;
  const G4double radius = 0.5*fibre_diameter;
  if(fibre_diameter <= 0. || fibre_diameter > 2.*mm)
    throw std::runtime_error("K18 BFT fibre diameter must be in (0,2] mm");
  if(black_sheet_thickness < 0.)
    throw std::runtime_error("K18 BFT black-sheet thickness must be non-negative");
  const G4double center_l = 0.5*(geomMan.GetLocalZ(names[0]) +
                                 geomMan.GetLocalZ(names[1]))*mm;
  const G4double x_offset_0 = geomMan.CalcWirePosition(names[0], 80)*mm;
  const G4double x_offset_1 = geomMan.CalcWirePosition(names[1], 80)*mm;
  const G4double dx = std::abs(x_offset_1 - x_offset_0);
  const G4double nonoverlap_separation =
    std::sqrt(std::max(0., fibre_diameter*fibre_diameter - dx*dx));
  const G4double physical_separation = std::max(
    std::abs(geomMan.GetLocalZ(names[1]) - geomMan.GetLocalZ(names[0]))*mm,
    nonoverlap_separation + 1.e-3*mm);

  // In the K18 internal frame the beam axis is native X, the measured BFT
  // coordinate is native Y, and native Z is vertical.  An unrotated G4Tubs
  // therefore already has its cylinder axis along the physical fibres.
  G4RotationMatrix fibre_rotation;
  for(G4int layer=0; layer<2; ++layer){
    const G4double l = center_l + (layer == 0 ? -0.5 : 0.5)*physical_separation;
    auto* solid = new G4Tubs("K18" + names[layer] + "FibreSolid",
                             0., radius, half_y, 0., CLHEP::twopi);
    auto* lv = new G4LogicalVolume(solid, fibre_material,
                                   "K18" + names[layer] + "FibreLV");
    lv->SetVisAttributes(plane_vis);
    G4int placed = 0;
    for(G4int wire=1; wire<=240; ++wire){
      const G4double x = geomMan.CalcWirePosition(names[layer], wire)*mm;
      if(std::abs(x) > half_x - radius)
        continue;
      const auto native_center = native_from_reference(K18Reference::VI, l)
        + x*K18BeamlineFrame::XAxis(0.);
      new G4PVPlacement(make_rotation_from_internal(frame, fibre_rotation),
                        frame.InternalToS2SPoint(native_center), lv,
                        "K18" + names[layer] + "FibrePV", world_lv, false,
                        wire, check_overlaps);
      ++placed;
    }
    if(placed == 0)
      throw std::runtime_error("K18 BFT realistic plane has no fibres");
    G4cout << "[K18BeamlineDetectorBuilder] " << names[layer]
           << " circular fibres: diameter=" << fibre_diameter/mm
           << " mm, placed=" << placed << ", physical_l=" << l/mm
           << " mm" << G4endl;
  }

  // Infinitesimal vacuum scoring planes provide a continuous BFT x coordinate
  // without double-counting scintillator.  They sit immediately outside the
  // fibre doublet and carry the two DCGEO detector IDs.
  const G4double proxy_half_z = 0.5e-3*mm;
  const G4double fibre_outer = 0.5*physical_separation + radius;
  for(G4int layer=0; layer<2; ++layer){
    const G4double sign = layer == 0 ? -1. : 1.;
    const G4double l = center_l + sign*(fibre_outer + proxy_half_z);
    place_box("K18" + names[layer] + "Proxy", half_x, half_y, proxy_half_z,
              native_from_reference(K18Reference::VI, l), 0., frame,
              world_lv, check_overlaps, plane_vis, proxy_material, sd,
              geomMan.GetDetectorId(names[layer]));
  }

  // The installed BFT is wrapped in an ordinary black PVC sheet.  Only the
  // upstream and downstream skins cross the beam; side wrapping is irrelevant
  // to the material budget.  The adopted 0.15-mm value is a common commercial
  // black-PVC sheet thickness, fixed rather than scanned.
  if(black_sheet_thickness > 0.){
    const G4double proxy_outer = fibre_outer + 2.*proxy_half_z;
    for(G4int side=0; side<2; ++side){
      const G4double sign = side == 0 ? -1. : 1.;
      const G4double l = center_l + sign*(proxy_outer +
                                           0.5*black_sheet_thickness);
      place_box("K18BFTBlackSheet" + std::to_string(side),
                half_x, half_y, 0.5*black_sheet_thickness,
                native_from_reference(K18Reference::VI, l), 0., frame,
                world_lv, check_overlaps, window_vis, black_sheet_material);
    }
  }
}

void place_bh1(const K18BeamlineFrame& frame,
               G4LogicalVolume* world_lv,
               G4bool check_overlaps,
               G4VisAttributes* vis,
               G4Material* material,
               G4VSensitiveDetector* sd)
{
  const auto widths = bh1_segment_widths();
  const G4double half_height = 0.5*66.*mm;
  const G4double half_thick = 0.5*5.*mm;
  const G4double overlap = 2.*mm;
  const G4double stagger = 3.*mm;
  const auto center = native_from_dcgeo(K18Reference::VI, "BH1");
  const auto ex = K18BeamlineFrame::XAxis(axis_angle(K18Reference::VI));
  const auto ez = K18BeamlineFrame::Axis(axis_angle(K18Reference::VI));

  std::vector<G4double> local_x(widths.size(), 0.);
  std::vector<G4double> local_z(widths.size(), 0.);
  const G4int center_seg = 5;
  local_z[center_seg] = -stagger;

  for(G4int i=center_seg - 1; i>=0; --i){
    local_x[i] = local_x[i + 1] - 0.5*(widths[i + 1] + widths[i]) + overlap;
    local_z[i] = (i%2 == 0) ? stagger : -stagger;
  }
  for(std::size_t i=center_seg + 1; i<widths.size(); ++i){
    local_x[i] = local_x[i - 1] + 0.5*(widths[i - 1] + widths[i]) - overlap;
    local_z[i] = (i%2 == 0) ? stagger : -stagger;
  }

  for(std::size_t i=0; i<widths.size(); ++i){
    place_box("K18BH1-" + std::to_string(i),
              0.5*widths[i], half_height, half_thick,
              center + local_x[i]*ex + local_z[i]*ez,
              axis_angle(K18Reference::VI), frame, world_lv, check_overlaps,
              vis, material, sd, static_cast<G4int>(i));
  }
}

void place_segmented_bh2(const K18BeamlineFrame& frame,
                         G4LogicalVolume* world_lv,
                         G4bool check_overlaps,
                         G4VisAttributes* vis,
                         G4Material* material,
                         G4VSensitiveDetector* sd)
{
  const auto size = sizeMan.GetSize("Bh2Seg")*mm;
  const auto widths = bh2_segment_widths();
  G4double total_width = 0.;
  for(const auto width : widths)
    total_width += width;

  const G4double longitudinal_shift =
    conf_double_or("K18BH2LongitudinalShift", 0.)*mm;
  const auto center = native_from_dcgeo(K18Reference::VO, "BH2")
    + longitudinal_shift*K18BeamlineFrame::Axis(axis_angle(K18Reference::VO));
  G4double x = -0.5*total_width;
  for(std::size_t i=0; i<widths.size(); ++i){
    const auto width = widths[i];
    const auto local_x = x + 0.5*width;
    x += width;
    place_box("K18BH2-" + std::to_string(i + 1),
              0.5*width, 0.5*size.y(), 0.5*size.z(),
              center + local_x*K18BeamlineFrame::XAxis(axis_angle(K18Reference::VO)),
              axis_angle(K18Reference::VO), frame, world_lv, check_overlaps,
              vis, material, sd, static_cast<G4int>(i + 1));
  }
}
}

K18BeamlineDetectorBuilder::K18BeamlineDetectorBuilder() = default;
K18BeamlineDetectorBuilder::~K18BeamlineDetectorBuilder() = default;

void
K18BeamlineDetectorBuilder::Construct(G4LogicalVolume* world_lv,
                                      G4bool check_overlaps)
{
  if(!world_lv)
    return;

  const G4bool all_vacuum = conf_bool_or("K18AllVacuum", false);
  auto* vacuum = all_vacuum
    ? G4Material::GetMaterial("Vacuum", true) : k18_vacuum_material();
  auto* air = G4Material::GetMaterial("Air", true);
  auto* iron = G4Material::GetMaterial("Fe", true);
  auto* copper = G4Material::GetMaterial("Cu", true);
  auto* bcout_gas =
    G4Material::GetMaterial("Ar76IsoButane20Methylal4Gas", true);
  auto* scintillator = G4Material::GetMaterial("Scintillator", true);
  auto* mylar = G4Material::GetMaterial("Mylar", true);
  auto* aluminum = G4Material::GetMaterial("Al", true);
  auto* carbon = G4Material::GetMaterial("C", true);
  auto* pvc = G4Material::GetMaterial("PVC", true);
  const G4bool vacuum_windows =
    !all_vacuum && conf_bool_or("K18VacuumWindows", false);
  const G4bool beam_pipe = conf_bool_or("K18BeamPipe", false);
  const G4double window_thickness =
    conf_double_or("K18VacuumWindowThickness", 0.5)*mm;
  const G4double window_radius =
    conf_double_or("K18VacuumWindowRadius", 100.)*mm;
  const G4bool realistic_bft =
    !all_vacuum && conf_bool_or("K18BFTRealisticMaterial", true);
  const G4double bft_fibre_diameter =
    conf_double_or("K18BFTFibreDiameter", 1.0)*mm;
  const G4bool bft_black_sheet =
    realistic_bft && conf_bool_or("K18BFTBlackSheet", true);
  const G4double bft_black_sheet_thickness = bft_black_sheet
    ? conf_double_or("K18BFTBlackSheetThickness", 0.15)*mm : 0.;
  const G4bool bcout_films =
    !all_vacuum && conf_bool_or("K18BcOutFilms", true);
  const G4double bcout_window_thickness =
    conf_double_or("K18BcOutWindowMylarThickness", 0.0125)*mm;
  const G4double bcout_cathode_thickness =
    conf_double_or("K18BcOutCathodeMylarThickness", 0.012)*mm;
  const G4double bcout_anode_cathode_gap =
    conf_double_or("K18BcOutAnodeCathodeGap", 2.0)*mm;
  // Both Al and evaporated carbon coatings are 5 nm per face by default.
  // Keep both values explicit and configurable.
  const G4double bcout_aluminum_coating_thickness =
    conf_double_or("K18BcOutAlCoatingThickness", 0.000005)*mm;
  const G4double bcout_carbon_coating_thickness =
    conf_double_or("K18BcOutCarbonCoatingThickness", 0.000005)*mm;
  const G4String window_material_name =
    conf_string_or("K18VacuumWindowMaterial", "SUS316L");
  auto* window_material = G4Material::GetMaterial(window_material_name, true);
  const G4String bh2_material_name =
    conf_string_or("K18BH2Material", "Scintillator");
  auto* configured_bh2_material =
    G4Material::GetMaterial(bh2_material_name, true);
  if(bcout_films && (bcout_window_thickness <= 0. ||
                     bcout_cathode_thickness <= 0. ||
                     bcout_anode_cathode_gap <= 0. ||
                     bcout_aluminum_coating_thickness < 0. ||
                     bcout_carbon_coating_thickness < 0.))
    throw std::runtime_error(
      "K18 BcOut film thicknesses must be positive when films are enabled");

  // Hall paint: Q yokes orange, D yokes dark green. Coils stay copper.
  auto* q_iron_vis = new G4VisAttributes(G4Colour(1.00, 0.42, 0.00, 0.82));
  q_iron_vis->SetForceSolid(true);
  auto* d_iron_vis = new G4VisAttributes(G4Colour(0.04, 0.32, 0.12, 0.85));
  d_iron_vis->SetForceSolid(true);
  auto* coil_vis = new G4VisAttributes(G4Colour(0.72, 0.45, 0.18, 0.70));
  coil_vis->SetForceSolid(true);
  auto* guard_vis = new G4VisAttributes(G4Colour(0.90, 0.32, 0.00, 0.90));
  guard_vis->SetForceSolid(true);
  auto* duct_vis = new G4VisAttributes(G4Colour(0.30, 0.65, 1.00, 0.30));
  duct_vis->SetForceWireframe(true);
  auto* chamber_vis = new G4VisAttributes(G4Colour(0.05, 0.70, 0.30, 0.24));
  chamber_vis->SetForceWireframe(true);
  auto* bc_plane_vis = new G4VisAttributes(G4Colour(0.10, 0.90, 0.20, 0.55));
  bc_plane_vis->SetForceSolid(true);
  auto* bft_plane_vis = new G4VisAttributes(G4Colour(1.00, 0.58, 0.05, 0.70));
  bft_plane_vis->SetForceSolid(true);
  auto* bh1_vis = new G4VisAttributes(G4Colour(0.95, 0.55, 0.15, 0.75));
  bh1_vis->SetForceSolid(true);
  auto* bh2_vis = new G4VisAttributes(G4Colour(0.95, 0.15, 0.15, 0.75));
  bh2_vis->SetForceSolid(true);
  auto* window_vis = new G4VisAttributes(G4Colour(0.65, 0.65, 0.70, 0.80));
  window_vis->SetForceSolid(true);
  auto* pipe_vis = new G4VisAttributes(G4Colour(0.20, 0.45, 0.95, 0.38));
  pipe_vis->SetForceSolid(true);

  auto* bc_sd = new DCSD("K18BC");
  G4SDManager::GetSDMpointer()->AddNewDetector(bc_sd);
  auto* bft_sd = new DCSD("K18BFT");
  G4SDManager::GetSDMpointer()->AddNewDetector(bft_sd);
  auto* bh1_sd = new TOFSD("K18BH1");
  G4SDManager::GetSDMpointer()->AddNewDetector(bh1_sd);
  auto* bh2_sd = new TOFSD("K18BH2");
  G4SDManager::GetSDMpointer()->AddNewDetector(bh2_sd);

  const auto target_global = geomMan.GetGlobalPosition("Target")*mm;
  const K18BeamlineFrame frame(target_global, target_l());

  const auto d4rot = d4_rotation();

  const G4double q10_l =
    -(d4_tan_half() + kDriftL2 + kQ11z + kDriftL1 + 0.5*kQ10z);
  const G4double vc1_l =
    -(d4_tan_half() + kDriftL2 + kQ11z + 0.5*kDriftL1);
  const G4double q11_l =
    -(d4_tan_half() + kDriftL2 + 0.5*kQ11z);
  const G4double vc2_l = -(d4_tan_half() + 0.5*kDriftL2);
  const G4double vc3_l = d4_tan_half() + 0.5*kDriftL3;
  const G4double q12_l = d4_tan_half() + kDriftL3 + 0.5*kQ12z;
  const G4double vc4_l = d4_tan_half() + kDriftL3 + kQ12z + 0.5*kDriftL4;
  const G4double q13_l =
    d4_tan_half() + kDriftL3 + kQ12z + kDriftL4 + 0.5*kQ13z;
  const G4double bend = K18BeamlineFrame::BendAngle();

  const G4double vi_l = K18BeamlineFrame::VinL();
  const G4double vo_l = K18BeamlineFrame::VoutL();
  const G4double q11_dn_outer = q11_l + kNQ412GuardTerm;
  const G4double d4_in_l = -d4_tan_half();
  const G4double d4_out_l = d4_tan_half();
  const G4double q12_up_outer = q12_l - kNQ412GuardTerm;
  const G4double disk_h = 0.5*kPipeDiskThick;
  // Switch plates sit just outside the Q11/Q12 terminal faces: a short
  // star stub, then the disk, then rectangle to D4.  Mid-gap disks crowd
  // the D4 aperture; a flush disk sits on the endguard.
  const G4double disk_in_l = q11_dn_outer + kPipeDiskStandoff + disk_h;
  const G4double disk_out_l = q12_up_outer - kPipeDiskStandoff - disk_h;
  const G4double win_h = (vacuum_windows && window_thickness > 0.)
    ? 0.5*window_thickness : 0.;

  // Connecting vacuum ducts stop at the magnet faces; the iron yokes own the
  // K18Q*PV and K18D4PV physical-volume names.
  if(beam_pipe){
    const auto star = star_opening_polygon(kPipeStarScale);
    const auto rect = rect_opening_polygon(kPipeD4HalfW, kPipeD4HalfH);
    G4int pcopy = 0;
    place_pipe_span("K18PipePV", star, vi_l + win_h, disk_in_l - disk_h, 0.,
                    frame, world_lv, check_overlaps, pipe_vis, vacuum, pcopy);
    place_pipe_disk("K18PipePV", disk_in_l, 0.,
                    frame, world_lv, check_overlaps, pipe_vis, window_vis,
                    vacuum, window_material, pcopy);
    place_pipe_span("K18PipePV", rect, disk_in_l + disk_h, d4_in_l, 0.,
                    frame, world_lv, check_overlaps, pipe_vis, vacuum, pcopy);
    const G4double phi_up = d4_from_fusion(d4_fusion_face_phi(kRhoD4));
    const G4double phi_dn = d4_from_fusion(0.);
    place_d4_tubs_at("K18PipePV",
      make_phi_tubs("K18PipeD4", kRhoD4 - kPipeD4HalfW, kRhoD4 + kPipeD4HalfW,
                    kPipeD4HalfH, phi_up, phi_dn),
      G4ThreeVector(-d4_tan_half(), kRhoD4, 0.), d4rot, frame, world_lv,
      check_overlaps, pipe_vis, vacuum, pcopy++);
    place_pipe_span("K18PipePV", rect, d4_out_l, disk_out_l - disk_h, bend,
                    frame, world_lv, check_overlaps, pipe_vis, vacuum, pcopy);
    place_pipe_disk("K18PipePV", disk_out_l, bend,
                    frame, world_lv, check_overlaps, pipe_vis, window_vis,
                    vacuum, window_material, pcopy);
    place_pipe_span("K18PipePV", star, disk_out_l + disk_h, vo_l - win_h, bend,
                    frame, world_lv, check_overlaps, pipe_vis, vacuum, pcopy);
    G4cout << "[K18BeamPipe] star scale=" << kPipeStarScale
           << " D4 rect ±" << kPipeD4HalfW/mm << " x ±" << kPipeD4HalfH/mm
           << " mm  face gap=" << kPipeDiskStandoff/mm
           << " mm  disk_in=" << (disk_in_l - q11_dn_outer)/mm
           << " mm after Q11 guard (D4 "
           << (d4_in_l - disk_in_l)/mm
           << " mm)  disk_out=" << (q12_up_outer - disk_out_l)/mm
           << " mm before Q12 guard (D4 "
           << (disk_out_l - d4_out_l)/mm
           << " mm)  hole R=" << kPipeDiskHoleR/mm
           << " mm, thick=" << kPipeDiskThick/mm
           << " mm  copies=" << pcopy << G4endl;
  }else{
    place_window("K18VC1", kQBoreR, kDriftL1,
                 on_straight_axis(vc1_l), 0.,
                 frame, world_lv, check_overlaps, duct_vis, vacuum);
    place_window("K18VC2", kQBoreR, kDriftL2,
                 on_straight_axis(vc2_l), 0.,
                 frame, world_lv, check_overlaps, duct_vis, vacuum);
    place_window("K18VC3", kQBoreR, kDriftL3,
                 on_bent_axis(vc3_l), bend,
                 frame, world_lv, check_overlaps, duct_vis, vacuum);
    place_window("K18VC4", kQBoreR, kDriftL4,
                 on_bent_axis(vc4_l), bend,
                 frame, world_lv, check_overlaps, duct_vis, vacuum);
  }

  if(!all_vacuum){
    place_quad_hardware("K18Q10", kQ10IronHZ, kQ10GuardTerm, kQ10GuardNon,
                        on_straight_axis(q10_l), 0.,
                        frame, world_lv, check_overlaps,
                        q_iron_vis, coil_vis, guard_vis, iron, copper);
    place_quad_hardware("K18Q11", kNQ412IronHZ, kNQ412GuardTerm, kNQ412GuardNon,
                        on_straight_axis(q11_l), 0.,
                        frame, world_lv, check_overlaps,
                        q_iron_vis, coil_vis, guard_vis, iron, copper);
    place_d4_hardware(G4ThreeVector(-d4_tan_half(), kRhoD4, 0.), d4rot,
                      frame, world_lv, check_overlaps,
                      d_iron_vis, coil_vis, iron, copper);
    // Q12/Q13: D4-facing (terminal) endguard is upstream along the exit axis.
    place_quad_hardware("K18Q12", kNQ412IronHZ,
                        -kNQ412GuardNon, -kNQ412GuardTerm,
                        on_bent_axis(q12_l), bend,
                        frame, world_lv, check_overlaps,
                        q_iron_vis, coil_vis, guard_vis, iron, copper);
    place_quad_hardware("K18Q13", kNQ412IronHZ,
                        -kNQ412GuardNon, -kNQ412GuardTerm,
                        on_bent_axis(q13_l), bend,
                        frame, world_lv, check_overlaps,
                        q_iron_vis, coil_vis, guard_vis, iron, copper);
  }

  (void)kQ10x; (void)kQ10y; (void)kQ11x; (void)kQ11y;
  (void)kQ12x; (void)kQ12y; (void)kQ13x; (void)kQ13y;

  if(vacuum_windows){
    if(window_thickness <= 0. || window_radius <= 0.){
      G4cerr << "[K18BeamlineDetectorBuilder] K18VacuumWindows requested but "
             << "thickness/radius is non-positive; windows are not placed"
             << G4endl;
    }else{
      // VI/VO are the QQDQQ optical ends (design 272/270 mm from Q10/Q13 EFB).
      // Point-to-point optics makes <x|θ>≈0, so MSC at these planes is
      // suppressed in reconstructed p (Yamamoto D / Ichikawa D).
      place_window("K18VIWindow", window_radius, window_thickness,
                   on_straight_axis(vi_l), 0.,
                   frame, world_lv, check_overlaps, window_vis, window_material);
      place_window("K18VOWindow", window_radius, window_thickness,
                   on_bent_axis(vo_l), bend,
                   frame, world_lv, check_overlaps, window_vis, window_material);
    }
  }

  auto* bft_mother_material = all_vacuum ? vacuum : air;
  auto* bcout_mother_material = all_vacuum ? vacuum : bcout_gas;
  auto* bh_material = all_vacuum ? vacuum : scintillator;
  auto* bh2_material = all_vacuum ? vacuum : configured_bh2_material;
  auto* bft_plane_material = all_vacuum ? vacuum : scintillator;
  auto* bc_plane_material = all_vacuum ? vacuum : bcout_gas;

  place_bh1(frame, world_lv, check_overlaps, bh1_vis, bh_material, bh1_sd);

  if(realistic_bft){
    place_k18_bft_realistic(frame, world_lv, check_overlaps,
                            bft_plane_vis, window_vis, scintillator, vacuum,
                            pvc, bft_sd, bft_fibre_diameter,
                            bft_black_sheet_thickness);
  }else{
    place_k18_chamber("K18BFT",
                      {{"BFT-X", 80.*mm, 40.*mm, 0.05*mm},
                       {"BFT-XP", 80.*mm, 40.*mm, 0.05*mm}}, {},
                      K18Reference::VI, 85.*mm, 45.*mm, 1.*mm,
                      frame, world_lv, check_overlaps, chamber_vis,
                      bft_plane_vis, bft_mother_material,
                      bft_plane_material, bft_sd);
  }

  const auto bc3_films = bcout_films ? make_bcout_films(
    "K18BC3", {{"BC3-X1", "BC3-X2"}, {"BC3-V1", "BC3-V2"},
                {"BC3-U1", "BC3-U2"}},
    100.*mm, 75.*mm, bcout_window_thickness,
    bcout_cathode_thickness, bcout_anode_cathode_gap,
    bcout_aluminum_coating_thickness,
    bcout_carbon_coating_thickness, mylar, aluminum, carbon)
    : std::vector<K18PassivePlaneSpec>{};
  const auto bc4_films = bcout_films ? make_bcout_films(
    "K18BC4", {{"BC4-U1", "BC4-U2"}, {"BC4-V1", "BC4-V2"},
                {"BC4-X1", "BC4-X2"}},
    125.*mm, 75.*mm, bcout_window_thickness,
    bcout_cathode_thickness, bcout_anode_cathode_gap,
    bcout_aluminum_coating_thickness,
    bcout_carbon_coating_thickness, mylar, aluminum, carbon)
    : std::vector<K18PassivePlaneSpec>{};

  place_k18_chamber("K18BC3",
                    {{"BC3-X1", 100.*mm, 75.*mm, 0.005*mm},
                     {"BC3-X2", 100.*mm, 75.*mm, 0.005*mm},
                     {"BC3-V1", 100.*mm, 75.*mm, 0.005*mm},
                     {"BC3-V2", 100.*mm, 75.*mm, 0.005*mm},
                     {"BC3-U1", 100.*mm, 75.*mm, 0.005*mm},
                     {"BC3-U2", 100.*mm, 75.*mm, 0.005*mm}}, bc3_films,
                    K18Reference::VO, 175.*mm, 150.*mm, 43.*mm,
                    frame, world_lv, check_overlaps, chamber_vis,
                    bc_plane_vis, bcout_mother_material,
                    bc_plane_material, bc_sd);

  place_k18_chamber("K18BC4",
                    {{"BC4-U1", 125.*mm, 75.*mm, 0.005*mm},
                     {"BC4-U2", 125.*mm, 75.*mm, 0.005*mm},
                     {"BC4-V1", 125.*mm, 75.*mm, 0.005*mm},
                     {"BC4-V2", 125.*mm, 75.*mm, 0.005*mm},
                     {"BC4-X1", 125.*mm, 75.*mm, 0.005*mm},
                     {"BC4-X2", 125.*mm, 75.*mm, 0.005*mm}}, bc4_films,
                    K18Reference::VO, 175.*mm, 150.*mm, 43.*mm,
                    frame, world_lv, check_overlaps, chamber_vis,
                    bc_plane_vis, bcout_mother_material,
                    bc_plane_material, bc_sd,
                    conf_double_or("K18BC4LongitudinalShift", 0.)*mm);

  const G4bool bh2_enabled = conf_bool_or("K18BH2Enabled", true);
  if(bh2_enabled)
    place_segmented_bh2(frame, world_lv, check_overlaps, bh2_vis, bh2_material,
                        bh2_sd);

  G4cout << "[K18BeamlineDetectorBuilder] material mode="
         << (all_vacuum ? "all-vacuum" : "standard")
         << "; BH1/BH2/BFT/BC3/BC4="
         << (all_vacuum ? vacuum->GetName() : "nominal")
         << "; VI/VO " << window_material_name << " windows=";
  if(vacuum_windows && window_thickness > 0. && window_radius > 0.)
    G4cout << "on (thickness=" << window_thickness/mm
           << " mm, radius=" << window_radius/mm
           << " mm at VI/VO)";
  else
    G4cout << "off";
  G4cout << G4endl;
  G4cout << "[K18BeamlineDetectorBuilder] K18BeamPipe="
         << (beam_pipe ? "on (star Q / rect D4 / circular switch disks)" : "off (circular VC1-4)")
         << G4endl;
  G4cout << "[K18BeamlineDetectorBuilder] BH2 material="
         << (bh2_enabled ? bh2_material->GetName() : "removed") << G4endl;

  G4cout << "[K18BeamlineDetectorBuilder] BFT material: "
         << (realistic_bft ? "1-mm circular fibres" : "legacy 0.1-mm slabs")
         << "; black PVC skins="
         << (bft_black_sheet ? "on" : "off")
         << " (each=" << bft_black_sheet_thickness/mm << " mm); "
         << "BcOut films=" << (bcout_films ? "on" : "off")
         << " (window=" << bcout_window_thickness/mm
         << " mm x2/chamber, cathode=" << bcout_cathode_thickness/mm
         << " mm x6/chamber, anode-cathode gap="
         << bcout_anode_cathode_gap/mm << " mm, Al/C coatings="
         << bcout_aluminum_coating_thickness/mm << "/"
         << bcout_carbon_coating_thickness/mm << " mm/face)" << G4endl;

  G4cout << "[K18BeamlineDetectorBuilder] layout shifts along K1.8 axis: "
         << "BC4=" << conf_double_or("K18BC4LongitudinalShift", 0.)
         << " mm, BH2=" << conf_double_or("K18BH2LongitudinalShift", 0.)
         << " mm" << G4endl;

  G4cout << "[K18BeamlineDetectorBuilder] K1.8 QQDQQ from clean-assembly "
         << "scripts (NQ418 v25 block29, NQ412 v1 block115, D4 v32/v54): "
         << "4 clamp plates + 4 poles (5-slice 62->107 end taper), _1__S star "
         << "endguards, Q coil hull (filled pancake+endturn+R-bend per quadrant), D4 H-frame "
         << "960x166 with " << kD4YokeRadialSlices
         << "-slice 23.5 deg face, open slot (no aperture-filling taper), "
         << "D4 Main+Aux pack hull (one pancake, one bend envelope), "
         << "face-oriented 23.5 deg upstream after Y-mirror"
         << (all_vacuum ? " (omitted, K18AllVacuum)" : "; K18Q*PV / K18D4PV")
         << G4endl;
  G4cout << "[K18BeamlineDetectorBuilder] S2S centres mm: Q10="
         << frame.InternalToS2SPoint(on_straight_axis(q10_l))/mm
         << " Q11=" << frame.InternalToS2SPoint(on_straight_axis(q11_l))/mm
         << " D4=" << frame.InternalToS2SPoint(G4ThreeVector(-d4_tan_half(),
                                                            kRhoD4, 0.))/mm
         << " Q12=" << frame.InternalToS2SPoint(on_bent_axis(q12_l))/mm
         << " Q13=" << frame.InternalToS2SPoint(on_bent_axis(q13_l))/mm
         << G4endl;
  G4cout << "[K18BeamlineDetectorBuilder] S2S VI/VO/disks mm: VI="
         << frame.InternalToS2SPoint(on_straight_axis(vi_l))/mm
         << " disk_in=" << frame.InternalToS2SPoint(on_straight_axis(disk_in_l))/mm
         << " disk_out=" << frame.InternalToS2SPoint(on_bent_axis(disk_out_l))/mm
         << " VO=" << frame.InternalToS2SPoint(on_bent_axis(vo_l))/mm
         << " (VI=" << vi_l/mm << " VO=" << vo_l/mm << " mm internal)"
         << G4endl;
  G4cout << "[K18BeamlineDetectorBuilder] K1.8 native Q10-Q13, D4,"
         << " VC1-VC4, BH1, BFT, BC3/BC4, and BH2 geometry added to G4S2S world"
         << G4endl;
}
