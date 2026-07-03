// -*- C++ -*-

#include "S2SDetectorConstruction.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <sstream>
#include <utility>
#include <vector>

#include <G4Box.hh>
#include <G4GenericTrap.hh>
#include <G4IntersectionSolid.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4PVPlacement.hh>
#include <G4Polyhedra.hh>
#include <G4QuadrangularFacet.hh>
#include <G4RotationMatrix.hh>
#include <G4SubtractionSolid.hh>
#include <G4SystemOfUnits.hh>
#include <G4TessellatedSolid.hh>
#include <G4TriangularFacet.hh>
#include <G4Trd.hh>
#include <G4Torus.hh>
#include <G4Tubs.hh>
#include <G4TwoVector.hh>
#include <G4UnionSolid.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>
#include <G4ios.hh>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "GeSD.hh"

namespace
{
const auto& confMan = ConfMan::GetInstance();
const auto& geomMan = DCGeomMan::GetInstance();

struct Axes
{
  G4ThreeVector x;
  G4ThreeVector y;
  G4ThreeVector z;
};

struct CloverSpec
{
  G4int sl;
  std::array<G4String, 4> ge_names;
  G4ThreeVector u_axis;
  G4ThreeVector detector_axis;
  G4double bgo_window_to_ge_window_tip;
};

struct SingleSpec
{
  G4int sl;
  G4String name;
  G4ThreeVector detector_axis;
  G4ThreeVector v0_axis;
  G4double tip_to_ge_center;
  G4double tip_to_inner_wall;
};

struct LocalPlacement
{
  G4LogicalVolume* lv;
  G4String phys_name;
  G4ThreeVector pos;
  G4RotationMatrix rot;
  G4int copy_number;
};

const G4double kCloverGePitch = 41.*mm;
const G4double kCloverGeHalfPitch = 0.5*kCloverGePitch;
const G4double kCloverGeLength = 71.*mm;
const G4double kGeShieldTipToGeFront = 20.*mm;
const G4double kCloverGeTemplateWarnShift = 0.05*mm;

G4Material*
mat(const G4String& name)
{
  return G4Material::GetMaterial(name, /*warning*/ true);
}

G4String
name_with_id(const G4String& base, G4int id)
{
  std::ostringstream os;
  os << base << id;
  return os.str();
}

G4VisAttributes*
make_solid_vis(const G4Colour& color)
{
  auto att = new G4VisAttributes(color);
  att->SetForceSolid(true);
  att->SetForceAuxEdgeVisible(true);
  return att;
}

Axes
make_axes(G4ThreeVector x_axis, G4ThreeVector z_axis)
{
  Axes axes;
  axes.z = z_axis.unit();
  axes.x = x_axis.unit();
  axes.y = axes.z.cross(axes.x).unit();
  axes.x = axes.y.cross(axes.z).unit();
  return axes;
}

G4RotationMatrix
make_rotation(const Axes& axes)
{
  G4RotationMatrix rot;
  rot.rotateAxes(axes.x, axes.y, axes.z);
  return rot;
}

std::array<G4TwoVector, 4>
hex_side_quad(G4double inner_apothem, G4double outer_apothem, G4int side)
{
  const G4double theta = (side*60. + 30.)*deg;
  const G4double half_inner = inner_apothem*std::tan(30.*deg);
  const G4double half_outer = outer_apothem*std::tan(30.*deg);
  const G4TwoVector normal(std::cos(theta), std::sin(theta));
  const G4TwoVector tangent(-std::sin(theta), std::cos(theta));

  return {{
    inner_apothem*normal - half_inner*tangent,
    inner_apothem*normal + half_inner*tangent,
    outer_apothem*normal + half_outer*tangent,
    outer_apothem*normal - half_outer*tangent
  }};
}

G4VSolid*
make_hex_side_panel_solid(const G4String& name, G4double length,
                          G4double inner_front, G4double outer_front,
                          G4double inner_back, G4double outer_back,
                          G4int side)
{
  const auto front = hex_side_quad(inner_front, outer_front, side);
  const auto back = hex_side_quad(inner_back, outer_back, side);
  std::vector<G4TwoVector> vertices;
  vertices.reserve(8);
  vertices.insert(vertices.end(), front.begin(), front.end());
  vertices.insert(vertices.end(), back.begin(), back.end());
  return new G4GenericTrap(name, length/2., vertices);
}

G4double
clover_crystal_local_z_rotation(const G4ThreeVector& local_pos)
{
  const G4bool xp = local_pos.x() > 0.;
  const G4bool yp = local_pos.y() > 0.;
  if(!xp && !yp) return 0.*deg;
  if(!xp &&  yp) return 270.*deg;
  if( xp &&  yp) return 180.*deg;
  return 90.*deg;
}

G4ThreeVector
to_local(const G4ThreeVector& global, const G4ThreeVector& origin,
         const Axes& axes)
{
  const G4ThreeVector d = global - origin;
  return G4ThreeVector(d.dot(axes.x), d.dot(axes.y), d.dot(axes.z));
}

G4ThreeVector
to_global(const G4ThreeVector& local, const G4ThreeVector& origin,
          const Axes& axes)
{
  return origin + local.x()*axes.x + local.y()*axes.y + local.z()*axes.z;
}

std::array<G4String, 4>
surface_names(const CloverSpec& spec)
{
  std::array<G4String, 4> names;
  for(G4int i=0; i<4; ++i)
    names[i] = spec.ge_names[i] + "-SURF";
  return names;
}

G4VSolid*
make_clover_ge_crystal_solid(const G4String& prefix)
{
  auto ge_base = new G4Tubs(prefix + "GeBase",
                            0.*mm, 50.*mm/2, kCloverGeLength/2,
                            0.*deg, 360.*deg);
  auto trd4 = new G4Box(prefix + "CutBox", 45.*mm/2, 45.*mm/2, 80.*mm/2);
  auto subtract2 = new G4SubtractionSolid(prefix + "GeCutY",
                                          ge_base, trd4, nullptr,
                                          G4ThreeVector(0, 43.*mm, 0));
  auto subtract3 = new G4SubtractionSolid(prefix + "GeCutYX",
                                          subtract2, trd4, nullptr,
                                          G4ThreeVector(43.*mm, 0, 0));
  auto trd5 = new G4Trd(prefix + "TaperCut",
                        50.*mm/2, 41.*mm/2,
                        50.*mm/2, 41.*mm/2,
                        36.1*mm/2);
  auto subtract4 = new G4SubtractionSolid(prefix + "GeCutTaperY",
                                          subtract3, trd5, nullptr,
                                          G4ThreeVector(0, -45.5*mm, -17.05*mm));
  return new G4SubtractionSolid(prefix + "GeCrystal",
                                subtract4, trd5, nullptr,
                                G4ThreeVector(-45.5*mm, 0, -17.05*mm));
}

G4ThreeVector
ideal_clover_ge_local_position(const CloverSpec& spec,
                               const G4ThreeVector& measured_local)
{
  const G4double sx = measured_local.x() >= 0. ? 1. : -1.;
  const G4double sy = measured_local.y() >= 0. ? 1. : -1.;
  return G4ThreeVector(sx*kCloverGeHalfPitch,
                       sy*kCloverGeHalfPitch,
                       spec.bgo_window_to_ge_window_tip
                       + kGeShieldTipToGeFront
                       + 0.5*kCloverGeLength);
}

G4LogicalVolume*
make_clover_ge_crystal_lv(const G4String& prefix, GeSD* ge_sd,
                          const G4Colour& color)
{
  auto lv = new G4LogicalVolume(make_clover_ge_crystal_solid(prefix),
                                mat("Ge"),
                                prefix + "LV");
  lv->SetSensitiveDetector(ge_sd);
  lv->SetVisAttributes(make_solid_vis(color));
  return lv;
}

std::vector<G4ThreeVector>
chamfered_square_vertices(G4double half_width, G4double chamfer,
                          G4double z)
{
  const G4double c = std::min(chamfer, 0.9*half_width);
  if(c <= 0.)
    return {
      G4ThreeVector( half_width,  half_width, z),
      G4ThreeVector(-half_width,  half_width, z),
      G4ThreeVector(-half_width, -half_width, z),
      G4ThreeVector( half_width, -half_width, z)
    };

  constexpr G4int kArcSegments = 4;
  std::vector<G4ThreeVector> vertices;
  vertices.reserve(4*kArcSegments + 4);
  auto add_arc = [&](G4double cx, G4double cy,
                     G4double a0, G4double a1,
                     G4int first_step, G4int last_step) {
    for(G4int i=first_step; i<=last_step; ++i){
      const G4double t = static_cast<G4double>(i)/kArcSegments;
      const G4double angle = a0 + (a1 - a0)*t;
      vertices.emplace_back(cx + c*std::cos(angle),
                            cy + c*std::sin(angle), z);
    }
  };

  vertices.emplace_back( half_width - c,  half_width,     z);
  vertices.emplace_back(-half_width + c,  half_width,     z);
  add_arc(-half_width + c,  half_width - c,  90.*deg, 180.*deg, 1, kArcSegments);
  vertices.emplace_back(-half_width,     -half_width + c, z);
  add_arc(-half_width + c, -half_width + c, 180.*deg, 270.*deg, 1, kArcSegments);
  vertices.emplace_back( half_width - c, -half_width,     z);
  add_arc( half_width - c, -half_width + c, 270.*deg, 360.*deg, 1, kArcSegments);
  vertices.emplace_back( half_width,      half_width - c, z);
  add_arc( half_width - c,  half_width - c,   0.*deg,  90.*deg, 1, kArcSegments - 1);
  return vertices;
}

G4VSolid*
make_chamfered_square_envelope_solid(
  const G4String& name,
  const std::vector<std::pair<G4double, G4double>>& z_half_width,
  G4double chamfer)
{
  auto solid = new G4TessellatedSolid(name);
  std::vector<std::vector<G4ThreeVector>> planes;
  planes.reserve(z_half_width.size());
  for(const auto& section : z_half_width)
    planes.push_back(chamfered_square_vertices(section.second, chamfer,
                                               section.first));

  for(std::size_t iz=0; iz+1<planes.size(); ++iz){
    for(std::size_t iv=0; iv<planes[iz].size(); ++iv){
      const std::size_t jv = (iv + 1) % planes[iz].size();
      solid->AddFacet(new G4QuadrangularFacet(planes[iz][iv],
                                              planes[iz][jv],
                                              planes[iz + 1][jv],
                                              planes[iz + 1][iv],
                                              ABSOLUTE));
    }
  }

  const auto& front = planes.front();
  const G4ThreeVector front_center(0., 0., z_half_width.front().first);
  for(std::size_t iv=0; iv<front.size(); ++iv){
    const std::size_t jv = (iv + 1) % front.size();
    solid->AddFacet(new G4TriangularFacet(front_center, front[jv],
                                          front[iv], ABSOLUTE));
  }

  const auto& back = planes.back();
  const G4ThreeVector back_center(0., 0., z_half_width.back().first);
  for(std::size_t iv=0; iv<back.size(); ++iv){
    const std::size_t jv = (iv + 1) % back.size();
    solid->AddFacet(new G4TriangularFacet(back_center, back[iv],
                                          back[jv], ABSOLUTE));
  }

  solid->SetSolidClosed(true);
  return solid;
}

std::vector<LocalPlacement>
make_clover_shield_placements(const G4String& prefix, GeSD* bgo_sd,
                              G4int sl, G4int bgo_copy_offset)
{
  const G4int num_sides = 4;
  const G4double bgo_phi_start = 45.*deg;
  // CloverGeBGO_size.pdf gives the whole shield placement reference, while
  // outer.pdf gives the physical SUS outer-wall dimensions.  The slot-dependent
  // BGO-tip to Ge-front offset is applied when placing this unit, not baked into
  // the BGO shield material dimensions.
  const G4double front_case_z = 6.*mm;
  const G4double case_taper_z = 47.6*mm;
  const G4double lower_taper_end_z = front_case_z + case_taper_z;
  const G4double body_end_z = 188.*mm;
  const G4double rear_flange_z = 8.*mm;
  const G4double rear_flange_end_z = body_end_z + rear_flange_z;
  const G4double bgo_crystal_z = 200.*mm;
  const G4double bgo_crystal_end_z = front_case_z + bgo_crystal_z;
  const G4double mother_end_z =
    bgo_crystal_end_z > rear_flange_end_z ? bgo_crystal_end_z : rear_flange_end_z;
  const G4double case_tip_outer_side = 117.1*mm;
  const G4double case_tip_outer_ro = 0.5*case_tip_outer_side;
  const G4double front_opening_ro = 0.5*92.0*mm;
  const G4double rear_flange_inner_ro = 177.4*mm/2.;
  const G4double case_wall_thickness = 3.*mm;
  const G4double upper_taper_angle = 83.*deg;
  const G4double case_front_end_inner_side = 113.0*mm;
  const G4double case_lower_taper_inner_side = 141.8*mm;
  const G4double case_lower_taper_outer_side_default =
    case_lower_taper_inner_side
    + 2.*case_wall_thickness/std::sin(upper_taper_angle);
  const G4double case_lower_taper_outer_side =
    case_lower_taper_outer_side_default;
  const G4double case_front_end_outer_side_default =
    case_tip_outer_side
    + (case_lower_taper_outer_side - case_tip_outer_side)
    * front_case_z/lower_taper_end_z;
  const G4double case_front_end_outer_side =
    case_front_end_outer_side_default;
  const G4double case_front_end_outer_ro =
    0.5*case_front_end_outer_side;
  const G4double case_lower_taper_outer_ro =
    0.5*case_lower_taper_outer_side;
  const G4double case_rear_outer_ro = 0.5*181.4*mm;
  const G4double case_front_end_inner_ro = 0.5*case_front_end_inner_side;
  const G4double case_lower_taper_inner_ro = 0.5*case_lower_taper_inner_side;
  const G4double case_rear_inner_side_default =
    2.*case_rear_outer_ro - 2.*case_wall_thickness/std::sin(upper_taper_angle);
  const G4double case_rear_inner_ro = 0.5*case_rear_inner_side_default;
  const G4double rear_flange_outer_ro = 0.5*234.0*mm;
  const G4double protect_box_outer_side = 240.0*mm;
  const G4double protect_box_outer_ro = 0.5*protect_box_outer_side;
  const G4double protect_box_wall_thickness = 6.0*mm;
  const G4double protect_box_inner_ro =
    protect_box_outer_ro - protect_box_wall_thickness;
  const G4double protect_front_clearance_ro =
    std::max(protect_box_inner_ro, rear_flange_outer_ro + 0.05*mm);
  const G4double protect_z_start = body_end_z;
  const G4double protect_z_length = 255.0*mm;
  const G4double bgo_front_inner_ro = 92.0*mm/2.;
  const G4double bgo_front_outer_ro = 0.5*112.0*mm;
  const G4double bgo_rear_inner_ro = 69.4*mm;
  const G4double bgo_rear_outer_ro = 176.8*mm/2.;
  const G4double bgo_crystal_center_half_width = 53.1*mm/2.;
  const G4double bgo_bar_gap = 0.10*mm;
  const G4double outer_edge_miter = 0.20*mm;
  const G4double front_edge_round = 0.20*mm;
  const G4double mag_shield_thickness = 9.0*mm;
  const G4double mag_shield_z_start = 98.0*mm;
  const G4double mag_shield_z_end = body_end_z;
  auto case_side_outer_ro_nominal = [&](G4double z) {
    if(z <= front_case_z){
      return case_tip_outer_ro
        + (case_front_end_outer_ro - case_tip_outer_ro)
        * z/front_case_z;
    }
    if(z <= lower_taper_end_z){
      return case_front_end_outer_ro
        + (case_lower_taper_outer_ro - case_front_end_outer_ro)
        * (z - front_case_z)/(lower_taper_end_z - front_case_z);
    }
    return case_lower_taper_outer_ro
      + (case_rear_outer_ro - case_lower_taper_outer_ro)
      * (z - lower_taper_end_z)/(body_end_z - lower_taper_end_z);
  };
  auto case_side_inner_ro = [&](G4double z) {
    if(z <= front_case_z){
      return front_opening_ro
        + (case_front_end_inner_ro - front_opening_ro)
        * z/front_case_z;
    }
    if(z <= lower_taper_end_z){
      return case_front_end_inner_ro
        + (case_lower_taper_inner_ro - case_front_end_inner_ro)
        * (z - front_case_z)/(lower_taper_end_z - front_case_z);
    }
    return case_lower_taper_inner_ro
      + (case_rear_inner_ro - case_lower_taper_inner_ro)
      * (z - lower_taper_end_z)/(body_end_z - lower_taper_end_z);
  };
  auto case_side_outer_ro = [&](G4double z) {
    const G4double ro = case_side_outer_ro_nominal(z);
    if(front_edge_round <= 0.*mm || z >= front_edge_round)
      return ro;
    return std::max(case_side_inner_ro(z) + 0.1*mm,
                    ro - (front_edge_round - z));
  };
  auto bgo_outer_ro = [&](G4double z) {
    const G4double bz = std::min(std::max(z - front_case_z, 0.*mm),
                                 bgo_crystal_z);
    return bgo_front_outer_ro
      + (bgo_rear_outer_ro - bgo_front_outer_ro)*bz/bgo_crystal_z;
  };
  auto envelope_outer_ro = [&](G4double z) {
    if(z < body_end_z)
      return case_side_outer_ro(z);
    if(z <= rear_flange_end_z)
      return rear_flange_outer_ro;
    return bgo_outer_ro(z);
  };
  std::vector<std::pair<G4double, G4double>> envelope_z = {{
    {0.*mm, envelope_outer_ro(0.*mm)},
    {front_case_z, envelope_outer_ro(front_case_z)},
    {lower_taper_end_z, envelope_outer_ro(lower_taper_end_z)},
    {body_end_z, case_side_outer_ro(body_end_z)},
    {body_end_z, envelope_outer_ro(body_end_z)},
    {rear_flange_end_z, envelope_outer_ro(rear_flange_end_z)},
    {mother_end_z, envelope_outer_ro(mother_end_z)}
  }};
  if(front_edge_round > 0.*mm){
    envelope_z.insert(envelope_z.begin() + 1,
                      {front_edge_round,
                       envelope_outer_ro(front_edge_round)});
  }
  auto make_envelope = [&](const G4String& name) {
    return make_chamfered_square_envelope_solid(name, envelope_z,
                                                outer_edge_miter);
  };
  auto clip_to_envelope =
    [&](const G4String& name, G4VSolid* solid,
        const G4RotationMatrix& child_rot,
        const G4ThreeVector& child_offset) -> G4VSolid* {
      const G4RotationMatrix inv_rot = child_rot.inverse();
      const G4ThreeVector env_offset = -(inv_rot*child_offset);
      return new G4IntersectionSolid(name, solid, make_envelope(name + "Envelope"),
                                     G4Transform3D(inv_rot, env_offset));
    };

  std::vector<LocalPlacement> placements;
  auto bgo_copy_number = [&](G4int side, G4int bar) {
    const G4int local = 3*side + bar;
    // The analyzer BtoG map numbers each clover's 12 BGO pieces by the
    // installed crystal view.  In the Geant local frame SL1/3/5 are rotated by
    // 180 degrees relative to that view, while SL2/4 already match it.
    const G4int rotated = (sl % 2 == 1) ? (local + 6) % 12 : local;
    return bgo_copy_offset + rotated;
  };

  {
    std::array<G4double, 5> z = {{0.*mm, front_edge_round,
                                  front_case_z, lower_taper_end_z,
                                  body_end_z}};
    G4int num_z = 5;
    if(front_edge_round <= 0.*mm){
      z = {{0.*mm, front_case_z, lower_taper_end_z, body_end_z, body_end_z}};
      num_z = 4;
    }
    std::array<G4double, 5> ri;
    std::array<G4double, 5> ro;
    for(G4int i=0; i<num_z; ++i){
      ro[i] = case_side_outer_ro(z[i]);
      ri[i] = case_side_inner_ro(z[i]);
    }
    auto raw_solid = new G4Polyhedra(prefix + "BGOOuterWallRawSolid",
                                     bgo_phi_start, 360.*deg,
                                     num_sides, num_z,
                                     z.data(), ri.data(), ro.data());
    G4RotationMatrix identity;
    auto solid = clip_to_envelope(prefix + "BGOOuterWallSolid",
                                  raw_solid, identity, G4ThreeVector());
    auto lv = new G4LogicalVolume(solid, mat("SUS316L"),
                                  prefix + "BGOOuterWallLV");
    lv->SetVisAttributes(make_solid_vis(G4Colour(1.0, 0.0, 1.0,
                                                 0.85)));
    placements.push_back({lv, prefix + "BGOOuterWallPHYS",
                          G4ThreeVector(), G4RotationMatrix(), 0});
  }

  {
    const G4ThreeVector bgo_pos(0, 0, front_case_z + 0.5*bgo_crystal_z);
    auto add_bgo_plane_vertices =
      [&](std::vector<G4TwoVector>& vertices,
          G4int bar, G4double inner_ro, G4double outer_ro) {
        const G4double center = bgo_crystal_center_half_width;
        const G4double gap = bgo_bar_gap;
        if(bar == 0){
          vertices.emplace_back(-inner_ro + gap, inner_ro);
          vertices.emplace_back(-outer_ro + gap, outer_ro);
          vertices.emplace_back(-center - 0.5*gap, outer_ro);
          vertices.emplace_back(-center - 0.5*gap, inner_ro);
        } else if(bar == 1){
          vertices.emplace_back(-center + 0.5*gap, inner_ro);
          vertices.emplace_back(-center + 0.5*gap, outer_ro);
          vertices.emplace_back( center - 0.5*gap, outer_ro);
          vertices.emplace_back( center - 0.5*gap, inner_ro);
        } else {
          vertices.emplace_back( center + 0.5*gap, inner_ro);
          vertices.emplace_back( center + 0.5*gap, outer_ro);
          vertices.emplace_back( outer_ro - gap, outer_ro);
          vertices.emplace_back( inner_ro - gap, inner_ro);
        }
      };

    for(G4int side=0; side<num_sides; ++side){
      G4RotationMatrix side_rot;
      side_rot.rotateZ(-side*90.*deg);
      for(G4int bar=0; bar<3; ++bar){
        std::vector<G4TwoVector> vertices;
        vertices.reserve(8);
        add_bgo_plane_vertices(vertices, bar,
                               bgo_front_inner_ro, bgo_front_outer_ro);
        add_bgo_plane_vertices(vertices, bar,
                               bgo_rear_inner_ro, bgo_rear_outer_ro);
        const G4String name =
          name_with_id(name_with_id(prefix + "BGOCrystalSide", side) + "Bar", bar);
        auto raw_solid = new G4GenericTrap(name + "RawSolid",
                                           0.5*bgo_crystal_z,
                                           vertices);
        auto solid = clip_to_envelope(name + "Solid",
                                      raw_solid, side_rot, bgo_pos);
        auto lv = new G4LogicalVolume(solid, mat("BGO"), name + "LV");
        if(bgo_sd)
          lv->SetSensitiveDetector(bgo_sd);
        lv->SetVisAttributes(make_solid_vis(G4Colour(0.0, 0.85, 0.95,
                                                     0.40)));
        placements.push_back({lv, name + "PHYS", bgo_pos, side_rot,
                              bgo_copy_number(side, bar)});
      }
    }
  }

  {
    G4double z[2] = {0.*mm, rear_flange_z};
    G4double ri[2] = {rear_flange_inner_ro, rear_flange_inner_ro};
    G4double ro[2] = {rear_flange_outer_ro, rear_flange_outer_ro};
    auto raw_solid = new G4Polyhedra(prefix + "BGORearFlangeRawSolid",
                                     bgo_phi_start, 360.*deg,
                                     num_sides, 2, z, ri, ro);
    G4RotationMatrix identity;
    auto solid = clip_to_envelope(prefix + "BGORearFlangeSolid",
                                  raw_solid, identity,
                                  G4ThreeVector(0, 0, body_end_z));
    auto lv = new G4LogicalVolume(solid, mat("SUS316L"),
                                  prefix + "BGORearFlangeLV");
    lv->SetVisAttributes(make_solid_vis(G4Colour(1.0, 0.0, 1.0,
                                                 0.85)));
    placements.push_back({lv, prefix + "BGORearFlangePHYS",
                          G4ThreeVector(0, 0, body_end_z),
                          G4RotationMatrix(), 0});
  }

  if(protect_box_wall_thickness > 0.1*mm && protect_z_length > 0.1*mm){
    std::array<G4double, 4> z = {{0.*mm, protect_z_length,
                                  protect_z_length, protect_z_length}};
    std::array<G4double, 4> ri = {{protect_box_inner_ro, protect_box_inner_ro,
                                   protect_box_inner_ro, protect_box_inner_ro}};
    std::array<G4double, 4> ro = {{protect_box_outer_ro, protect_box_outer_ro,
                                   protect_box_outer_ro, protect_box_outer_ro}};
    G4int num_z = 2;
    const G4double flange_relief_z =
      std::clamp(rear_flange_end_z - protect_z_start,
                 0.*mm, protect_z_length);
    if(flange_relief_z > 0.1*mm
       && protect_front_clearance_ro > protect_box_inner_ro + 0.1*mm){
      const G4double relief_step_z =
        std::min(protect_z_length, flange_relief_z + 0.05*mm);
      z = {{0.*mm, flange_relief_z, relief_step_z, protect_z_length}};
      ri = {{protect_front_clearance_ro, protect_front_clearance_ro,
             protect_box_inner_ro, protect_box_inner_ro}};
      num_z = 4;
    }
    auto solid = new G4Polyhedra(prefix + "BGOProtectBoxSolid",
                                 bgo_phi_start, 360.*deg,
                                 num_sides, num_z, z.data(), ri.data(), ro.data());
    auto lv = new G4LogicalVolume(solid, mat("SUS316L"),
                                  prefix + "BGOProtectBoxLV");
    lv->SetVisAttributes(make_solid_vis(G4Colour(0.95, 0.78, 0.10,
                                                 0.85)));
    placements.push_back({lv, prefix + "BGOProtectBoxPHYS",
                          G4ThreeVector(0, 0, protect_z_start),
                          G4RotationMatrix(), 0});
  }

  if(mag_shield_thickness > 0.*mm
     && mag_shield_z_end > mag_shield_z_start + 0.1*mm){
    G4double z[2] = {mag_shield_z_start, mag_shield_z_end};
    G4double ri[2] = {
      case_side_outer_ro(mag_shield_z_start) + 0.05*mm,
      case_side_outer_ro(mag_shield_z_end) + 0.05*mm
    };
    G4double ro[2] = {
      ri[0] + mag_shield_thickness,
      ri[1] + mag_shield_thickness
    };
    auto solid = new G4Polyhedra(prefix + "MagShieldRawSolid",
                                 bgo_phi_start, 360.*deg,
                                 num_sides, 2, z, ri, ro);
    auto lv = new G4LogicalVolume(solid, mat("Fe"),
                                  prefix + "MagShieldLV");
    lv->SetVisAttributes(make_solid_vis(G4Colour(0.20, 0.20, 0.20,
                                                 0.85)));
    placements.push_back({lv, prefix + "MagShieldPHYS",
                          G4ThreeVector(), G4RotationMatrix(), 0});
  }

  return placements;
}

struct SingleGeDimensions
{
  G4double diameter = 0.;
  G4double length = 0.;
};

const G4double kSingleGeReferenceRadius = (67.8/2.0)*mm;
const G4double kSingleGeReferenceLength = 76.1*mm;
const G4double kSingleGeReferenceEdge = 15.*mm;
const G4double kSingleGeReferenceAnodeRadius = 4.5*mm;
const G4double kSingleGeReferenceAnodeLength =
  0.7*kSingleGeReferenceLength;
const G4double kSingleGeLegacyOriginOffset =
  0.5*kSingleGeReferenceEdge;

SingleGeDimensions
single_ge_dimensions()
{
  return {2.0*kSingleGeReferenceRadius, kSingleGeReferenceLength};
}

G4ThreeVector
single_ge_placement_offset()
{
  return G4ThreeVector(0., 0., kSingleGeLegacyOriginOffset);
}

G4VSolid*
make_single_ge_solid(const G4String& prefix, const SingleGeDimensions&)
{
  auto solid_ge1 =
    new G4Tubs(prefix + "GeDet1",
               0.*mm, kSingleGeReferenceRadius,
               0.5*(kSingleGeReferenceLength - kSingleGeReferenceEdge),
               0.*deg, 360.*deg);

  auto solid_ge2 =
    new G4Tubs(prefix + "GeDet2",
               0.*mm, kSingleGeReferenceAnodeRadius,
               0.5*kSingleGeReferenceAnodeLength,
               0.*deg, 360.*deg);

  auto solid_ge3 =
    new G4Tubs(prefix + "GeDet3",
               0.*mm, kSingleGeReferenceRadius - kSingleGeReferenceEdge,
               0.05*mm + 0.5*kSingleGeReferenceEdge,
               0.*deg, 360.*deg);

  auto solid_ge4 =
    new G4Torus(prefix + "GeDet4",
                0.*mm,
                kSingleGeReferenceEdge - 0.001*mm,
                kSingleGeReferenceRadius - kSingleGeReferenceEdge,
                0.*deg,
                360.*deg);

  auto solid_ge_uni1 =
    new G4UnionSolid(prefix + "GeDet_uni1",
                     solid_ge1, solid_ge4, nullptr,
                     G4ThreeVector(0., 0.,
                                   -0.5*(kSingleGeReferenceLength
                                          - kSingleGeReferenceEdge)));
  auto solid_ge_uni2 =
    new G4UnionSolid(prefix + "GeDet_uni2",
                     solid_ge_uni1, solid_ge3, nullptr,
                     G4ThreeVector(0., 0.,
                                   -0.5*kSingleGeReferenceLength + 0.1*mm));
  return new G4SubtractionSolid(
    prefix + "SingleGe",
    solid_ge_uni2, solid_ge2, nullptr,
    G4ThreeVector(0., 0.,
                  0.5*(0.3*kSingleGeReferenceLength
                       - kSingleGeReferenceEdge) + 0.005*mm));
}

void
construct_single_unit(const G4String& prefix, GeSD* ge_sd, GeSD* bgo_sd,
                      const SingleSpec& spec, G4int ge_copy_number,
                      G4int bgo_copy_offset,
                      const G4ThreeVector& center, const Axes& axes,
                      G4LogicalVolume* world_lv, G4bool check_overlaps)
{
  constexpr G4int num_sides = 6;
  const G4double bgo_depth = 150.*mm;       // BGO crystal drawing
  const G4double bgo_taper_depth = 50.*mm;  // 19 mm body tapers to 9 mm tip
  const G4double bgo_front_offset = 3.*mm;  // Stainless tip thickness before the BGO crystal.
  const G4double bgo_inner_apothem = 45.6*mm;
  const G4double bgo_tip_outer_apothem = bgo_inner_apothem + 9.0*mm;
  const G4double bgo_body_outer_apothem = bgo_inner_apothem + 19.0*mm;
  const G4double case_depth = 152.5*mm;
  const G4double case_taper_depth = 50.5*mm;
  const G4double case_wall = 2.0*mm;
  const G4double case_tip_outer_apothem = bgo_tip_outer_apothem + case_wall;
  const G4double case_body_inner_apothem = 65.5*mm;
  const G4double case_body_outer_apothem = 67.5*mm;
  const G4double inner_wall_depth = 134.*mm;
  const G4double inner_wall_inner_radius = 43.5*mm;
  const G4double inner_wall_outer_radius = 44.5*mm;
  const G4double tip_z = -spec.tip_to_ge_center;
  const G4RotationMatrix unit_rot = make_rotation(axes);
  auto place = [&](G4LogicalVolume* lv, const G4String& phys_name,
                   const G4ThreeVector& local_pos, G4int copy_number) {
    new G4PVPlacement(G4Transform3D(unit_rot,
                                    to_global(local_pos, center, axes)),
                      lv, phys_name, world_lv,
                      false, copy_number, check_overlaps);
  };

  const G4double case_clearance = 0.05*mm; // Avoid coincident case/BGO side boundaries in navigation.
  for(G4int side=0; side<num_sides; ++side){
    auto taper_solid =
      make_hex_side_panel_solid(name_with_id(prefix + "SingleBGOCaseTaper", side),
                                case_taper_depth,
                                bgo_tip_outer_apothem + case_clearance,
                                case_tip_outer_apothem,
                                case_body_inner_apothem,
                                case_body_outer_apothem,
                                side);
    auto taper_lv = new G4LogicalVolume(taper_solid, mat("SUS316L"),
                                        name_with_id(prefix + "SingleBGOCaseTaperLV", side));
    taper_lv->SetVisAttributes(make_solid_vis(G4Colour(0.75, 0.75, 0.75,
                                                       0.45)));
    place(taper_lv, name_with_id(prefix + "SingleBGOCaseTaperPHYS", side),
          G4ThreeVector(0., 0., tip_z + 0.5*case_taper_depth),
          ge_copy_number);

    auto body_solid =
      make_hex_side_panel_solid(name_with_id(prefix + "SingleBGOCaseBody", side),
                                case_depth - case_taper_depth,
                                case_body_inner_apothem,
                                case_body_outer_apothem,
                                case_body_inner_apothem,
                                case_body_outer_apothem,
                                side);
    auto body_lv = new G4LogicalVolume(body_solid, mat("SUS316L"),
                                       name_with_id(prefix + "SingleBGOCaseBodyLV", side));
    body_lv->SetVisAttributes(make_solid_vis(G4Colour(0.75, 0.75, 0.75,
                                                      0.45)));
    place(body_lv, name_with_id(prefix + "SingleBGOCaseBodyPHYS", side),
          G4ThreeVector(0., 0.,
                        tip_z + case_taper_depth
                        + 0.5*(case_depth - case_taper_depth)),
          ge_copy_number);
  }

  for(G4int side=0; side<num_sides; ++side){
    auto taper_solid =
      make_hex_side_panel_solid(name_with_id(prefix + "SingleBGOCrystalTaper", side),
                                bgo_taper_depth,
                                bgo_inner_apothem, bgo_tip_outer_apothem,
                                bgo_inner_apothem, bgo_body_outer_apothem,
                                side);
    auto taper_lv = new G4LogicalVolume(taper_solid, mat("BGO"),
                                        name_with_id(prefix + "SingleBGOCrystalTaperLV", side));
    if(bgo_sd)
      taper_lv->SetSensitiveDetector(bgo_sd);
    taper_lv->SetVisAttributes(make_solid_vis(G4Colour(0.0, 0.85, 0.95,
                                                       0.45)));
    place(taper_lv, name_with_id(prefix + "SingleBGOCrystalTaperPHYS", side),
          G4ThreeVector(0., 0., tip_z + bgo_front_offset + 0.5*bgo_taper_depth),
          bgo_copy_offset + side);

    auto body_solid =
      make_hex_side_panel_solid(name_with_id(prefix + "SingleBGOCrystalBody", side),
                                bgo_depth - bgo_taper_depth,
                                bgo_inner_apothem, bgo_body_outer_apothem,
                                bgo_inner_apothem, bgo_body_outer_apothem,
                                side);
    auto body_lv = new G4LogicalVolume(body_solid, mat("BGO"),
                                       name_with_id(prefix + "SingleBGOCrystalBodyLV", side));
    if(bgo_sd)
      body_lv->SetSensitiveDetector(bgo_sd);
    body_lv->SetVisAttributes(make_solid_vis(G4Colour(0.0, 0.85, 0.95,
                                                      0.45)));
    place(body_lv, name_with_id(prefix + "SingleBGOCrystalBodyPHYS", side),
          G4ThreeVector(0., 0.,
                        tip_z + bgo_front_offset + bgo_taper_depth
                        + 0.5*(bgo_depth - bgo_taper_depth)),
          bgo_copy_offset + side);
  }

  {
    const G4double inner_wall_z = tip_z + spec.tip_to_inner_wall;
    auto ring = new G4Tubs(prefix + "SingleInnerWall",
                           inner_wall_inner_radius, inner_wall_outer_radius,
                           0.5*inner_wall_depth, 0.*deg, 360.*deg);
    auto lv = new G4LogicalVolume(ring, mat("Al"),
                                  prefix + "SingleInnerWallLV");
    lv->SetVisAttributes(make_solid_vis(G4Colour(0.95, 0.95, 0.15,
                                                 0.5)));
    place(lv, prefix + "SingleInnerWallTipPHYS",
          G4ThreeVector(0., 0., inner_wall_z + 0.5*inner_wall_depth),
          ge_copy_number);
  }

  const SingleGeDimensions ge_dim = single_ge_dimensions();
  auto ge_lv = new G4LogicalVolume(make_single_ge_solid(prefix, ge_dim),
                                   mat("Ge"),
                                   prefix + "SingleGeLV");
  ge_lv->SetSensitiveDetector(ge_sd);
  ge_lv->SetVisAttributes(make_solid_vis(G4Colour(1., 0., 1.)));
  place(ge_lv, prefix + "SingleGePHYS", single_ge_placement_offset(),
        ge_copy_number);
}

G4ThreeVector
average_position(const std::array<G4String, 4>& names)
{
  G4ThreeVector sum;
  for(const auto& name : names)
    sum += geomMan.GetGlobalPosition(name);
  return sum / 4.;
}

std::array<G4int, 4>
clover_u_signs(G4int sl)
{
  switch(sl){
  case 1:
    return {{-1, -1, +1, +1}};
  case 2:
    return {{+1, +1, -1, -1}};
  case 3:
    return {{+1, -1, -1, +1}};
  case 4:
    return {{-1, +1, +1, -1}};
  case 5:
    return {{+1, -1, -1, +1}};
  default:
    return {{+1, -1, -1, +1}};
  }
}

G4ThreeVector
derive_clover_u_axis(const CloverSpec& spec, const G4ThreeVector& center)
{
  const auto signs = clover_u_signs(spec.sl);
  G4ThreeVector u;
  for(G4int i=0; i<4; ++i)
    u += signs[i]*(geomMan.GetGlobalPosition(spec.ge_names[i]) - center);
  return u.mag2() > 0. ? u.unit() : spec.u_axis.unit();
}

G4ThreeVector
derive_clover_detector_axis(const CloverSpec& spec,
                            const G4ThreeVector& center,
                            const G4ThreeVector& surface_center)
{
  const G4ThreeVector axis = center - surface_center;
  return axis.mag2() > 0. ? axis.unit() : spec.detector_axis.unit();
}

const std::array<CloverSpec, 5> kClovers = {{
  {1,
   {{"HBXX-SL1-GE1", "HBXX-SL1-GE2", "HBXX-SL1-GE3", "HBXX-SL1-GE4"}},
   G4ThreeVector(0.8642463751760788, 0.5030687855502549, 0.0),
   G4ThreeVector(-0.5030687855502549, 0.8642463751760788, 0.0),
   9.0*mm},
  {2,
   {{"HBXX-SL2-GE1", "HBXX-SL2-GE2", "HBXX-SL2-GE3", "HBXX-SL2-GE4"}},
   G4ThreeVector(0.8578769252742173, -0.5138552141236429, 0.0),
   G4ThreeVector(-0.5138552141236429, -0.8578769252742173, 0.0),
   10.0*mm},
  {3,
   {{"HBXX-SL3-GE1", "HBXX-SL3-GE2", "HBXX-SL3-GE3", "HBXX-SL3-GE4"}},
   G4ThreeVector(-0.8582726594949202, 0.5131939613474781, 0.0),
   G4ThreeVector(0.5131939613474781, 0.8582726594949202, 0.0),
   9.0*mm},
  {4,
   {{"HBXX-SL4-GE1", "HBXX-SL4-GE2", "HBXX-SL4-GE3", "HBXX-SL4-GE4"}},
   G4ThreeVector(-0.8559155419457903, -0.5171156399255817, 0.0),
   G4ThreeVector(0.5171156399255817, -0.8559155419457903, 0.0),
   9.5*mm},
  {5,
   {{"HBXX-SL5-GE1", "HBXX-SL5-GE2", "HBXX-SL5-GE3", "HBXX-SL5-GE4"}},
   G4ThreeVector(0.8454312550, 0.3116994020, -0.4336928350),
   G4ThreeVector(-0.5140863940, 0.6950568980, -0.5026043050),
   10.0*mm}
}};

const std::array<G4String, 3> kSingleNames = {{
  "HBXX-SL6", "HBXX-SL7", "HBXX-SL8"
}};

const std::array<SingleSpec, 3> kSingles = {{
  {6, "HBXX-SL6",
   G4ThreeVector(-0.469441273, -0.512392198, -0.719082142),
   G4ThreeVector(-0.707539279, 0.705494961, -0.040804760),
   55.0*mm, 11.0*mm},
  {7, "HBXX-SL7",
   G4ThreeVector(0.484431345, 0.540590068, -0.687814401),
   G4ThreeVector(0.725790910, -0.687303406, -0.029010047),
   57.0*mm, 13.0*mm},
  {8, "HBXX-SL8",
   G4ThreeVector(0.502069619, -0.579126960, -0.642291259),
   G4ThreeVector(0.737215653, 0.674887553, -0.032247045),
   57.0*mm, 13.0*mm}
}};
}


//_____________________________________________________________________________
void
S2SDetectorConstruction::ConstructHBXXGe()
{
  if(m_experiment != 63)
    return;
  try {
    geomMan.GetDetectorId(kClovers.front().ge_names.front());
  } catch(const std::exception& e) {
    G4cout << "[HBXXGe] skipped: " << e.what() << G4endl;
    return;
  }

  auto ge_sd = new GeSD("Ge");
  AddNewDetector(ge_sd);
  auto bgo_sd = new GeSD("BGO");
  AddNewDetector(bgo_sd);

  for(const auto& spec : kClovers){
    const G4ThreeVector center = average_position(spec.ge_names);
    const auto surf_names = surface_names(spec);
    const G4ThreeVector surface_center = average_position(surf_names);
    G4ThreeVector u_axis = derive_clover_u_axis(spec, center);
    G4ThreeVector detector_axis =
      derive_clover_detector_axis(spec, center, surface_center);
    const Axes axes = make_axes(u_axis, detector_axis);
    // GE*-SURF is the active Ge crystal front.  The clover drawing gives the
    // BGO shield tip to Ge shield tip gap per slot, then a 20 mm setback from
    // the Ge shield tip to the active Ge crystal front.
    const G4double bgo_tip_to_ge_front =
      spec.bgo_window_to_ge_window_tip + kGeShieldTipToGeFront;
    G4ThreeVector unit_origin = surface_center - bgo_tip_to_ge_front*axes.z;
    const G4String prefix = name_with_id("HBXXSL", spec.sl);
    const auto unit_rot = make_rotation(axes);
    auto place_local =
      [&](G4LogicalVolume* lv, const G4String& phys_name,
          const G4ThreeVector& local_pos, const G4RotationMatrix& local_rot,
          G4int copy_number) {
        const G4RotationMatrix global_rot = unit_rot*local_rot;
        const G4ThreeVector global_pos = unit_origin + unit_rot*local_pos;
        new G4PVPlacement(G4Transform3D(global_rot, global_pos),
                          lv, phys_name, m_world_lv,
                          false, copy_number, m_check_overlaps);
      };

    for(const auto& placement :
        make_clover_shield_placements(prefix, bgo_sd, spec.sl,
                                      (spec.sl - 1)*12)){
      place_local(placement.lv, placement.phys_name, placement.pos,
                  placement.rot, placement.copy_number);
    }

    for(G4int i=0; i<4; ++i){
      const G4int seg = (spec.sl - 1)*4 + i;
      const G4ThreeVector global_pos = geomMan.GetGlobalPosition(spec.ge_names[i]);
      const G4ThreeVector measured_local_pos =
        to_local(global_pos, unit_origin, axes);
      const G4ThreeVector local_pos =
        ideal_clover_ge_local_position(spec, measured_local_pos);
      const G4double template_shift = (local_pos - measured_local_pos).mag();
      if(template_shift > kCloverGeTemplateWarnShift){
        G4cout << "[HBXXGe] warning: " << spec.ge_names[i]
               << " DCGEO center differs from rigid clover template by "
               << template_shift/mm << " mm"
               << " measured_local=" << measured_local_pos/mm
               << " template_local=" << local_pos/mm
               << G4endl;
      }
      const G4double rot_z = clover_crystal_local_z_rotation(local_pos);
      auto ge_lv = make_clover_ge_crystal_lv(name_with_id(prefix + "Ge", i),
                                             ge_sd,
                                             G4Colour(0.95 - 0.15*i,
                                                      0.20 + 0.18*i,
                                                      0.25 + 0.12*i));
      G4RotationMatrix local_rot;
      local_rot.rotateZ(rot_z);
      place_local(ge_lv, "HBXXGePHYS", local_pos, local_rot, seg);
    }

    G4cout << "[HBXXGe] clover SL" << spec.sl
           << " center=" << center/mm
           << " surface=" << surface_center/mm
           << " bgo_tip_to_ge_shield_tip="
           << spec.bgo_window_to_ge_window_tip/mm
           << " ge_shield_tip_to_ge_front="
           << kGeShieldTipToGeFront/mm
           << " front_origin=" << unit_origin/mm
           << " axis=" << axes.z
           << " u=" << axes.x
           << " mm" << G4endl;
  }

  for(G4int i=0; i<static_cast<G4int>(kSingles.size()); ++i){
    const auto& spec = kSingles[i];
    const G4int seg = 20 + i;
    const G4ThreeVector center = geomMan.GetGlobalPosition(spec.name);
    const Axes axes = make_axes(spec.v0_axis, spec.detector_axis);
    construct_single_unit(name_with_id("HBXXSingle", seg),
                          ge_sd, bgo_sd, spec, seg, 60 + i*6,
                          center, axes,
                          m_world_lv, m_check_overlaps);
    G4cout << "[HBXXGe] single " << kSingleNames[i]
           << " seg=" << seg
           << " center=" << center/mm
           << " axis=" << spec.detector_axis
           << " tip=" << (center - spec.tip_to_ge_center*axes.z)/mm
           << " single_ge_diameter=" << single_ge_dimensions().diameter/mm
           << " mm single_ge_length=" << single_ge_dimensions().length/mm
           << " mm" << G4endl;
  }
}
