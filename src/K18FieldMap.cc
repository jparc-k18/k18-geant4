// -*- C++ -*-

#include "K18FieldMap.hh"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>

#include <CLHEP/Units/SystemOfUnits.h>
#include <G4ios.hh>

#include "FuncName.hh"

namespace
{
G4bool
GridIndex(G4double value, G4double origin, G4double step, G4int count,
          G4int& index)
{
  const G4double grid = (value-origin)/step;
  const long long nearest = std::llround(grid);
  const G4double tolerance = std::max(1.e-8, 1.e-6*std::abs(step));
  if(nearest < 0 || nearest >= count ||
     std::abs(value-(origin+nearest*step)) > tolerance)
    return false;
  index = static_cast<G4int>(nearest);
  return true;
}

G4bool
AxisWeights(G4double value, G4double origin, G4double step, G4int count,
            G4int& i0, G4int& i1, G4double& w0, G4double& w1)
{
  const G4double last = origin+(count-1)*step;
  const G4double tolerance = std::max(1.e-10, 1.e-9*std::abs(step));
  if(value < origin-tolerance || value > last+tolerance)
    return false;

  const G4double grid = std::max(0., std::min(G4double(count-1),
                                               (value-origin)/step));
  if(grid >= count-1){
    i0 = i1 = count-1;
    w0 = 1.;
    w1 = 0.;
    return true;
  }

  i0 = static_cast<G4int>(std::floor(grid));
  i1 = i0+1;
  w1 = grid-i0;
  w0 = 1.-w1;
  return true;
}
}

//_____________________________________________________________________________
K18FieldMap::K18FieldMap(const G4String& file_name,
                         G4double value_measure, G4double value_calc)
  : m_is_ready(false),
    m_file_name(file_name),
    m_field(),
    m_nx(0), m_ny(0), m_nz(0),
    m_xmin(0.), m_ymin(0.), m_zmin(0.),
    m_dx(0.), m_dy(0.), m_dz(0.),
    m_scale(value_calc != 0. ? value_measure/value_calc
                             : std::numeric_limits<G4double>::quiet_NaN())
{
}

//_____________________________________________________________________________
std::size_t
K18FieldMap::Index(G4int ix, G4int iy, G4int iz) const
{
  return (static_cast<std::size_t>(ix)*m_ny+iy)*m_nz+iz;
}

//_____________________________________________________________________________
void
K18FieldMap::Clear()
{
  m_is_ready = false;
  m_field.clear();
}

//_____________________________________________________________________________
G4bool
K18FieldMap::Initialize()
{
  if(m_is_ready){
    G4cerr << FUNC_NAME << " already initialized" << G4endl;
    return false;
  }
  if(!std::isfinite(m_scale)){
    G4cerr << FUNC_NAME << " invalid K18FLDNMR/K18FLDCALC scale" << G4endl;
    return false;
  }

  std::ifstream ifs(m_file_name);
  if(!ifs.is_open()){
    G4cerr << FUNC_NAME << " file open fail: " << m_file_name << G4endl;
    return false;
  }

  Clear();
  if(!(ifs >> m_nx >> m_ny >> m_nz
       >> m_xmin >> m_ymin >> m_zmin >> m_dx >> m_dy >> m_dz)){
    G4cerr << FUNC_NAME << " invalid header: " << m_file_name << G4endl;
    return false;
  }
  if(m_nx < 2 || m_ny < 2 || m_nz < 2 ||
     !std::isfinite(m_dx) || !std::isfinite(m_dy) ||
     !std::isfinite(m_dz) || m_dx <= 0. || m_dy <= 0. || m_dz <= 0.){
    G4cerr << FUNC_NAME << " invalid grid dimensions or spacing" << G4endl;
    return false;
  }

  const std::size_t nx = static_cast<std::size_t>(m_nx);
  const std::size_t ny = static_cast<std::size_t>(m_ny);
  const std::size_t nz = static_cast<std::size_t>(m_nz);
  if(nx > std::numeric_limits<std::size_t>::max()/ny ||
     nx*ny > std::numeric_limits<std::size_t>::max()/nz){
    G4cerr << FUNC_NAME << " grid size overflow" << G4endl;
    return false;
  }
  const std::size_t expected = nx*ny*nz;
  m_field.assign(expected, G4ThreeVector());
  std::vector<std::uint8_t> seen(expected, 0);

  G4cout << FUNC_NAME << " file = " << m_file_name << G4endl
         << "  grid = " << m_nx << " x " << m_ny << " x " << m_nz
         << ", origin(cm) = (" << m_xmin << ", " << m_ymin
         << ", " << m_zmin << ")"
         << ", step(cm) = (" << m_dx << ", " << m_dy
         << ", " << m_dz << ")" << G4endl
         << "  uniform field scale K18FLDNMR/K18FLDCALC = "
         << m_scale << G4endl
         << " reading K18 fieldmap " << std::flush;

  G4double x, y, z, bx, by, bz;
  std::size_t records = 0;
  while(ifs >> x >> y >> z >> bx >> by >> bz){
    G4int ix, iy, iz;
    if(!GridIndex(x, m_xmin, m_dx, m_nx, ix) ||
       !GridIndex(y, m_ymin, m_dy, m_ny, iy) ||
       !GridIndex(z, m_zmin, m_dz, m_nz, iz)){
      G4cerr << G4endl << FUNC_NAME << " off-grid record at ("
             << x << ", " << y << ", " << z << ") cm" << G4endl;
      Clear();
      return false;
    }
    const std::size_t index = Index(ix, iy, iz);
    if(seen[index]){
      G4cerr << G4endl << FUNC_NAME << " duplicate grid record at ("
             << x << ", " << y << ", " << z << ") cm" << G4endl;
      Clear();
      return false;
    }
    seen[index] = 1;
    m_field[index].set(bx*m_scale, by*m_scale, bz*m_scale);
    ++records;
    if(records%1000000 == 0)
      G4cout << "." << std::flush;
  }

  if(!ifs.eof() || records != expected ||
     std::find(seen.begin(), seen.end(), std::uint8_t(0)) != seen.end()){
    G4cerr << G4endl << FUNC_NAME << " incomplete or malformed map: "
           << records << " records, expected " << expected << G4endl;
    Clear();
    return false;
  }

  G4cout << " done (" << records << " records)" << G4endl;
  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
G4bool
K18FieldMap::GetFieldValue(const G4double point[4], G4double* bfield) const
{
  bfield[0] = bfield[1] = bfield[2] = 0.;
  if(!m_is_ready)
    return false;

  const G4double position_cm[3] = {
    point[0]/CLHEP::cm, point[1]/CLHEP::cm, point[2]/CLHEP::cm
  };
  G4int index0[3], index1[3];
  G4double weight0[3], weight1[3];
  if(!AxisWeights(position_cm[0], m_xmin, m_dx, m_nx,
                  index0[0], index1[0], weight0[0], weight1[0]) ||
     !AxisWeights(position_cm[1], m_ymin, m_dy, m_ny,
                  index0[1], index1[1], weight0[1], weight1[1]) ||
     !AxisWeights(position_cm[2], m_zmin, m_dz, m_nz,
                  index0[2], index1[2], weight0[2], weight1[2]))
    return false;

  G4ThreeVector field;
  for(G4int cx=0; cx<2; ++cx){
    const G4int ix = cx ? index1[0] : index0[0];
    const G4double wx = cx ? weight1[0] : weight0[0];
    for(G4int cy=0; cy<2; ++cy){
      const G4int iy = cy ? index1[1] : index0[1];
      const G4double wy = cy ? weight1[1] : weight0[1];
      for(G4int cz=0; cz<2; ++cz){
        const G4int iz = cz ? index1[2] : index0[2];
        const G4double wz = cz ? weight1[2] : weight0[2];
        field += wx*wy*wz*m_field[Index(ix, iy, iz)];
      }
    }
  }
  bfield[0] = field.x()*CLHEP::tesla;
  bfield[1] = field.y()*CLHEP::tesla;
  bfield[2] = field.z()*CLHEP::tesla;
  return true;
}
