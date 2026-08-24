// -*- C++ -*-

#include "K18BeamlineField.hh"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

#include <G4SystemOfUnits.hh>
#include <G4ios.hh>
#include <G4ThreeVector.hh>

#include "ConfMan.hh"
#include "K18FieldMap.hh"
#include "S2SField.hh"

namespace
{
const G4double kDeg2Rad = std::acos(-1.)/180.;
const G4double kRad2Deg = 180./std::acos(-1.);

// K1.8 native magnet constants. Magnet sizes/fields follow the standalone
// K1.8 beamline field geometry.  VI/VO drift lengths follow the configured
// detector-coordinate geometry.
const G4double kRhoD4       = 4.000*m;
const G4double kBendAngleD4 = 64.000;
const G4double kD4SpaceX    = 0.657*m;
const G4double kD4HalfGap   = 0.1*m;
const G4double kD4Length    = 4.468*m;
const G4double kD4B0        = 11.5137;
const G4double kAlphaD4     = 23.5;

const G4double kQ10x  = 0.89*m;
const G4double kQ10y  = 0.89*m;
const G4double kQ10z  = 0.956*m;
const G4double kQ10B0 = -7.30734;
const G4double kQ10a0 = 0.1*m;

const G4double kQ11x  = 0.89*m;
const G4double kQ11y  = 0.89*m;
const G4double kQ11z  = 0.66*m;
const G4double kQ11B0 = 5.60272;
const G4double kQ11a0 = 0.1*m;

const G4double kQ12x  = 1.095*m;
const G4double kQ12y  = 1.095*m;
const G4double kQ12z  = 0.66*m;
const G4double kQ12B0 = -5.99934;
const G4double kQ12a0 = 0.1*m;

const G4double kQ13x  = 1.095*m;
const G4double kQ13y  = 1.095*m;
const G4double kQ13z  = 0.66*m;
const G4double kQ13B0 = 5.00944;
const G4double kQ13a0 = 0.1*m;

const G4double kDriftL0 = 0.272*m; // VI to Q10 upstream EFB (magnet surface +0.102)
const G4double kDriftL1 = 0.242*m;
const G4double kDriftL2 = 0.470*m;
const G4double kDriftL3 = 0.470*m;
const G4double kDriftL4 = 0.240*m;
const G4double kDriftL5 = 0.270*m; // Q13 downstream EFB to VO (magnet surface +0.100)

inline G4double d4_tan_half()
{
  return kRhoD4*std::tan(0.5*kBendAngleD4*kDeg2Rad);
}

inline G4ThreeVector internal_on_axis(G4double angle_deg, G4double l)
{
  return l*K18BeamlineFrame::Axis(angle_deg);
}

G4double conf_double_or(const G4String& key, G4double fallback)
{
  const G4String raw = ConfMan::GetInstance().Get<G4String>(key);
  return raw.empty() ? fallback : ConfMan::GetInstance().Get<G4double>(key);
}

G4bool field_map_enabled(const G4String& value)
{
  return !value.empty() && value != "none" && value != "skip" && value != "0";
}

std::shared_ptr<const K18FieldMap> configured_k18_field_map()
{
  static const std::shared_ptr<const K18FieldMap> map = [](){
    const auto& conf = ConfMan::GetInstance();
    const G4String file_name = conf.Get<G4String>("K18FLDMAP");
    if(!field_map_enabled(file_name))
      return std::shared_ptr<const K18FieldMap>();

    const G4double value_nmr = conf_double_or("K18FLDNMR", 1.);
    const G4double value_calc = conf_double_or("K18FLDCALC", 1.);
    auto loaded = std::make_shared<K18FieldMap>(file_name,
                                                value_nmr, value_calc);
    if(!loaded->Initialize())
      throw std::runtime_error("failed to initialize K18FLDMAP: " + file_name);
    return std::shared_ptr<const K18FieldMap>(loaded);
  }();
  return map;
}

G4String strip_k18_prefix(const G4String& name)
{
  const G4String prefix("K18");
  if(name.compare(0, prefix.size(), prefix) == 0)
    return name.substr(prefix.size());
  return name;
}

G4double magnet_scale(const G4String& name, const G4String& scale_prefix)
{
  G4double global_scale = conf_double_or("K18GlobalScale", 1.0);
  G4double local_scale = conf_double_or(name + "Scale", 1.0);
  if(!scale_prefix.empty()){
    global_scale = conf_double_or(scale_prefix + "GlobalScale", global_scale);
    local_scale = conf_double_or(scale_prefix + strip_k18_prefix(name)
                                 + "Scale", local_scale);
    local_scale = conf_double_or(scale_prefix + name + "Scale", local_scale);
  }
  const G4double scale = global_scale*local_scale;
  if(std::abs(scale - 1.0) > 1.e-12){
    G4cout << "[K18BeamlineField] ";
    if(!scale_prefix.empty())
      G4cout << scale_prefix << " ";
    G4cout << name << " scale=" << scale << G4endl;
  }
  return scale;
}
}

K18BeamlineFrame::K18BeamlineFrame(const G4ThreeVector& target_global,
                                   G4double target_l)
  : m_target_global(target_global),
    m_target_internal(Vout() + target_l*Axis(BendAngle()))
{
}

G4double
K18BeamlineFrame::BendAngle()
{
  return kBendAngleD4;
}

G4ThreeVector
K18BeamlineFrame::Axis(G4double angle_deg)
{
  const G4double a = angle_deg*kDeg2Rad;
  return G4ThreeVector(std::cos(a), std::sin(a), 0.);
}

G4ThreeVector
K18BeamlineFrame::XAxis(G4double angle_deg)
{
  const G4double a = angle_deg*kDeg2Rad;
  return G4ThreeVector(-std::sin(a), std::cos(a), 0.);
}

G4ThreeVector
K18BeamlineFrame::YAxis()
{
  return G4ThreeVector(0., 0., 1.);
}

G4double
K18BeamlineFrame::VinL()
{
  return -(d4_tan_half() + kDriftL2 + kQ11z + kDriftL1 + kQ10z + kDriftL0);
}

G4ThreeVector
K18BeamlineFrame::Vin()
{
  return G4ThreeVector(VinL(), 0., 0.);
}

G4double
K18BeamlineFrame::Q10UpstreamOffset()
{
  return kDriftL0;
}

G4double
K18BeamlineFrame::VoutL()
{
  return d4_tan_half() + kDriftL3 + kQ12z + kDriftL4 + kQ13z + kDriftL5;
}

G4ThreeVector
K18BeamlineFrame::Vout()
{
  return internal_on_axis(BendAngle(), VoutL());
}

G4ThreeVector
K18BeamlineFrame::InternalToS2SPoint(const G4ThreeVector& point) const
{
  return m_target_global + InternalToS2SVector(point - m_target_internal);
}

G4ThreeVector
K18BeamlineFrame::S2SToInternalPoint(const G4ThreeVector& point) const
{
  return m_target_internal + S2SToInternalVector(point - m_target_global);
}

G4ThreeVector
K18BeamlineFrame::InternalToS2SVector(const G4ThreeVector& vec) const
{
  const auto ex = XAxis(BendAngle());
  const auto ey = YAxis();
  const auto ez = Axis(BendAngle());
  return G4ThreeVector(vec.dot(ex), vec.dot(ey), vec.dot(ez));
}

G4ThreeVector
K18BeamlineFrame::S2SToInternalVector(const G4ThreeVector& vec) const
{
  return vec.x()*XAxis(BendAngle()) + vec.y()*YAxis()
    + vec.z()*Axis(BendAngle());
}

K18BeamlineField::K18BeamlineField()
  : K18BeamlineField("")
{
}

K18BeamlineField::K18BeamlineField(const G4String& scale_prefix)
  : m_field_map(configured_k18_field_map())
{
  m_magnet[0].type = kQuad;
  m_magnet[0].name = "K18Q10";
  m_magnet[0].scale = magnet_scale("K18Q10", scale_prefix);
  m_magnet[0].theta = 0.;
  m_magnet[0].quad = {
    kQ10B0, kQ10a0, kQ10z,
    -(d4_tan_half() + kDriftL2 + kQ11z + kDriftL1 + 0.5*kQ10z)
  };

  m_magnet[1].type = kQuad;
  m_magnet[1].name = "K18Q11";
  m_magnet[1].scale = magnet_scale("K18Q11", scale_prefix);
  m_magnet[1].theta = 0.;
  m_magnet[1].quad = {
    kQ11B0, kQ11a0, kQ11z,
    -(d4_tan_half() + kDriftL2 + 0.5*kQ11z)
  };

  m_magnet[2].type = kSect;
  m_magnet[2].name = "K18D4";
  m_magnet[2].scale = magnet_scale("K18D4", scale_prefix);
  m_magnet[2].theta = 0.;
  m_magnet[2].sect = {
    kD4B0,
    kRhoD4,
    kD4SpaceX,
    kD4HalfGap,
    G4ThreeVector(-d4_tan_half(), kRhoD4, 0.),
    kBendAngleD4,
    kAlphaD4
  };

  m_magnet[3].type = kQuad;
  m_magnet[3].name = "K18Q12";
  m_magnet[3].scale = magnet_scale("K18Q12", scale_prefix);
  m_magnet[3].theta = kBendAngleD4;
  m_magnet[3].quad = {
    kQ12B0, kQ12a0, kQ12z,
    d4_tan_half() + kDriftL3 + 0.5*kQ12z
  };

  m_magnet[4].type = kQuad;
  m_magnet[4].name = "K18Q13";
  m_magnet[4].scale = magnet_scale("K18Q13", scale_prefix);
  m_magnet[4].theta = kBendAngleD4;
  m_magnet[4].quad = {
    kQ13B0, kQ13a0, kQ13z,
    d4_tan_half() + kDriftL3 + kQ12z + kDriftL4 + 0.5*kQ13z
  };

  (void)kQ10x;
  (void)kQ10y;
  (void)kQ11x;
  (void)kQ11y;
  (void)kQ12x;
  (void)kQ12y;
  (void)kQ13x;
  (void)kQ13y;
  (void)kD4Length;
}

K18BeamlineField::~K18BeamlineField() = default;

G4bool
K18BeamlineField::Contains(const G4ThreeVector& pos, G4int i) const
{
  return 0 <= i && i < 5 && InMagnet(pos, i);
}

G4bool
K18BeamlineField::ContainsNominal(const G4ThreeVector& pos, G4int i) const
{
  if(i < 0 || i >= 5)
    return false;
  const auto& mag = m_magnet[i];
  if(mag.type == kQuad)
    return InMagnet(pos, i);

  const auto& sec = mag.sect;
  const G4double lvecx = pos.x() - sec.center.x();
  const G4double lvecy = pos.y() - sec.center.y();
  const G4double r = std::sqrt(lvecx*lvecx + lvecy*lvecy);
  if(r <= 0.)
    return false;
  const G4double ctheta = std::max(-1., std::min(1., lvecx/r));
  G4double theta = std::acos(ctheta)*kRad2Deg;
  if(lvecy < 0.)
    theta = 360. - theta;
  const G4double entrance_line =
    -1./std::tan(sec.alpha*kDeg2Rad)
    *(pos.x() + kRhoD4*std::tan(0.5*sec.bend_angle*kDeg2Rad));
  return r >= sec.rho - 0.5*sec.width
    && r <= sec.rho + 0.5*sec.width
    && theta >= 270.
    && theta <= 270. + sec.bend_angle
    && std::abs(pos.z()) <= sec.half_gap
    && pos.y() >= entrance_line;
}

void
K18BeamlineField::GetFieldValue(const G4double point[4],
                                G4double* bfield) const
{
  if(m_field_map){
    m_field_map->GetFieldValue(point, bfield);
    return;
  }

  const G4ThreeVector pos(point[0], point[1], point[2]);
  G4double local_b[3] = {0., 0., 0.};

  for(G4int i=0; i<5; ++i){
    if(!InMagnet(pos, i))
      continue;
    if(m_magnet[i].type == kQuad)
      CalcQuad(i, pos, local_b);
    else
      CalcSect(i, pos, local_b);
    break;
  }

  bfield[0] = local_b[0]/10.*tesla;
  bfield[1] = local_b[1]/10.*tesla;
  bfield[2] = local_b[2]/10.*tesla;
}

G4bool
K18BeamlineField::InMagnet(const G4ThreeVector& pos, G4int i) const
{
  const auto& mag = m_magnet[i];
  if(mag.type == kQuad){
    auto gpos = pos;
    gpos.rotateZ(-mag.theta*deg);
    const G4double lx = gpos.y();
    const G4double ly = gpos.z();
    const G4double lz = gpos.x();
    return (std::abs(lz - mag.quad.z) <= 0.5*mag.quad.l
            && std::sqrt(lx*lx + ly*ly) <= 2.*mag.quad.a0);
  }

  const auto& sec = mag.sect;
  const G4double lvecx = pos.x() - sec.center.x();
  const G4double lvecy = pos.y() - sec.center.y();
  const G4double r = std::sqrt(lvecx*lvecx + lvecy*lvecy);
  if(r <= 0.)
    return false;
  const G4double ctheta = std::max(-1., std::min(1., lvecx/r));
  G4double theta = std::acos(ctheta)*kRad2Deg;
  if(lvecy < 0.)
    theta = 360. - theta;

  const G4double entrance_line =
    -1./std::tan(sec.alpha*kDeg2Rad)*(pos.x() + kRhoD4*std::tan(0.5*sec.bend_angle*kDeg2Rad));

  if(r >= sec.rho - 0.5*sec.width
     && r <= sec.rho + 0.5*sec.width
     && theta >= 270.
     && theta <= 270. + sec.bend_angle
     && std::abs(pos.z()) <= sec.half_gap){
    return pos.y() >= entrance_line;
  }

  if(pos.y() >= entrance_line
     && pos.x() <= -kRhoD4*std::tan(0.5*sec.bend_angle*kDeg2Rad)
     && pos.y() <= 0.5*sec.width){
    return true;
  }

  return false;
}

void
K18BeamlineField::CalcQuad(G4int i, const G4ThreeVector& pos,
                           G4double* bfield) const
{
  const auto& mag = m_magnet[i];
  auto gpos = pos;
  gpos.rotateZ(-mag.theta*deg);
  const G4double lx = gpos.y();
  const G4double ly = gpos.z();
  const G4double g = mag.scale*mag.quad.b0/mag.quad.a0;
  const G4double local_bx = g*ly;
  const G4double local_by = g*lx;
  G4ThreeVector global_b(0., local_bx, local_by);
  global_b.rotateZ(mag.theta*deg);
  bfield[0] = global_b.x();
  bfield[1] = global_b.y();
  bfield[2] = global_b.z();
}

void
K18BeamlineField::CalcSect(G4int i, const G4ThreeVector&,
                           G4double* bfield) const
{
  bfield[0] = 0.;
  bfield[1] = 0.;
  bfield[2] = m_magnet[i].scale*m_magnet[i].sect.b0;
}

K18BeamlineFieldAdapter::K18BeamlineFieldAdapter(
  const G4ThreeVector& target_global,
  G4double target_l)
  : m_frame(target_global, target_l),
    m_field()
{
}

K18BeamlineFieldAdapter::~K18BeamlineFieldAdapter() = default;

void
K18BeamlineFieldAdapter::GetFieldValue(const G4double point[4],
                                       G4double* bfield) const
{
  const G4ThreeVector s2s_pos(point[0], point[1], point[2]);
  const G4ThreeVector internal_pos = m_frame.S2SToInternalPoint(s2s_pos);
  G4double internal_point[4] =
    {internal_pos.x(), internal_pos.y(), internal_pos.z(), point[3]};
  G4double internal_b[3] = {0., 0., 0.};
  m_field.GetFieldValue(internal_point, internal_b);
  const G4ThreeVector b_s2s =
    m_frame.InternalToS2SVector(G4ThreeVector(internal_b[0],
                                             internal_b[1],
                                             internal_b[2]));
  bfield[0] = b_s2s.x();
  bfield[1] = b_s2s.y();
  bfield[2] = b_s2s.z();
}

K18CombinedField::K18CombinedField(const G4String& s2s_field_map,
                                   const G4ThreeVector& target_global,
                                   G4double target_l)
  : m_s2s_field(new S2SField(s2s_field_map)),
    m_k18_field(target_global, target_l)
{
}

K18CombinedField::~K18CombinedField() = default;

void
K18CombinedField::GetFieldValue(const G4double point[4],
                                G4double* bfield) const
{
  G4double s2s_b[3] = {0., 0., 0.};
  G4double k18_b[3] = {0., 0., 0.};
  m_s2s_field->GetFieldValue(point, s2s_b);
  m_k18_field.GetFieldValue(point, k18_b);
  bfield[0] = s2s_b[0] + k18_b[0];
  bfield[1] = s2s_b[1] + k18_b[1];
  bfield[2] = s2s_b[2] + k18_b[2];
}
