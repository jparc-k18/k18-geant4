// -*- C++ -*-

#include "DCGeomRecord.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>

#include <G4String.hh>

#include <TMath.h>

#include "FuncName.hh"
#include "PrintHelper.hh"

namespace
{
enum  EDefinition { kSks, kKurama };
const EDefinition GlobalCoordinate = kSks;
}

//_____________________________________________________________________________
DCGeomRecord::DCGeomRecord(G4int id, const G4String& name,
                           G4double x, G4double y, G4double z, G4double ta,
                           G4double ra1, G4double ra2, G4double length,
                           G4double resol,
                           G4double w0, G4double dd, G4double ofs)
  : m_id(id), m_name(name), m_pos(x,y,z), m_tilt_angle(ta),
    m_rot_angle1(ra1), m_rot_angle2(ra2),
    m_length(length), m_resolution(resol), m_w0(w0), m_dd(dd), m_offset(ofs)
{
  CalcVectors();
}

//_____________________________________________________________________________
DCGeomRecord::DCGeomRecord(G4int id, const G4String& name,
                           const G4ThreeVector pos, G4double ta,
                           G4double ra1, G4double ra2, G4double length,
                           G4double resol,
                           G4double w0, G4double dd, G4double ofs)
  : m_id(id), m_name(name), m_pos(pos),  m_tilt_angle(ta),
    m_rot_angle1(ra1), m_rot_angle2(ra2),
    m_length(length), m_resolution(resol), m_w0(w0), m_dd(dd), m_offset(ofs)
{
  CalcVectors();
}

//_____________________________________________________________________________
DCGeomRecord::~DCGeomRecord()
{
}

//_____________________________________________________________________________
void
DCGeomRecord::CalcVectors()
{
  G4double ct0 = TMath::Cos(m_tilt_angle*TMath::DegToRad());
  G4double st0 = TMath::Sin(m_tilt_angle*TMath::DegToRad());
  G4double ct1 = TMath::Cos(m_rot_angle1*TMath::DegToRad());
  G4double st1 = TMath::Sin(m_rot_angle1*TMath::DegToRad());
  G4double ct2 = TMath::Cos(m_rot_angle2*TMath::DegToRad());
  G4double st2 = TMath::Sin(m_rot_angle2*TMath::DegToRad());

  switch(GlobalCoordinate){
  case kSks: // SKS difinition
    m_dxds =  ct0*ct2-st0*ct1*st2;
    m_dxdt = -st0*ct2-ct0*ct1*st2;
    m_dxdu =  st1*st2;

    m_dyds =  ct0*st2+st0*ct1*ct2;
    m_dydt = -st0*st2+ct0*ct1*ct2;
    m_dydu = -st1*ct2;

    m_dzds =  st0*st1;
    m_dzdt =  ct0*st1;
    m_dzdu =  ct1;

    m_dsdx =  ct0*ct2-st0*ct1*st2;
    m_dsdy =  ct0*st2+st0*ct1*ct2;
    m_dsdz =  st0*st1;

    m_dtdx = -st0*ct2-ct0*ct1*st2;
    m_dtdy = -st0*st2+ct0*ct1*ct2;
    m_dtdz =  ct0*st1;

    m_dudx =  st1*st2;
    m_dudy = -st1*ct2;
    m_dudz =  ct1;
    break;
  case kKurama: // new difinition : RA2 is around Y-axis
    m_dxds =  ct0*ct2+st0*st1*st2;
    m_dxdt = -st0*ct2+ct0*st1*st2;
    m_dxdu =  ct1*st2;

    m_dyds =  st0*ct1;
    m_dydt =  ct0*ct1;
    m_dydu = -st1;

    m_dzds = -ct0*st2+st0*st1*ct2;
    m_dzdt =  st0*st2+ct0*st1*ct2;
    m_dzdu =  ct1*ct2;

    m_dsdx =  ct0*ct2+st0*st1*st2;
    m_dsdy =  st0*ct1;
    m_dsdz = -ct0*st2+st0*st1*ct2;

    m_dtdx = -st0*ct2+ct0*st1*st2;
    m_dtdy =  ct0*ct1;
    m_dtdz =  st0*st2+ct0*st1*ct2;

    m_dudx =  ct1*st2;
    m_dudy = -st1;
    m_dudz =  ct1*ct2;
    break;
  default:
    throw std::invalid_argument(G4String(FUNC_NAME +
                                         " Invalid global coordinate"));
  }
}

//_____________________________________________________________________________
G4ThreeVector
DCGeomRecord::NormalVector() const
{
  return G4ThreeVector(m_dxdu, m_dydu, m_dzdu);
}

//_____________________________________________________________________________
G4ThreeVector
DCGeomRecord::UnitVector() const
{
  return G4ThreeVector(m_dxds, m_dyds, m_dzds);
}

//_____________________________________________________________________________
G4double
DCGeomRecord::WirePos(G4double wire) const
{
  return m_dd*(wire - m_w0)+m_offset;
}

//_____________________________________________________________________________
G4int
DCGeomRecord::WireNumber(G4double pos) const
{
  G4double dw = ((pos-m_offset)/m_dd) + m_w0;
  G4int    iw = G4int(dw);
  if((dw-G4double(iw))>0.5)
    return iw+1;
  else
    return iw;
}

//_____________________________________________________________________________
void
DCGeomRecord::Print() const
{
  PrintHelper helper(3, std::ios::fixed, std::cout);
  std::cout << " id = "   << std::setw(3) << std::right << m_id << " "
	    << " name = " << std::setw(9) << std::left << m_name << " "
	    << " pos = (" << " " << std::right << std::fixed
	    << std::setw(8)  << m_pos.x() << ", "
	    << std::setw(8)  << m_pos.y() << ", "
	    << std::setw(10) << m_pos.z() << ") "
	    << " TA = "     << std::setw(7) << m_tilt_angle << " "
	    << " RA1 = "    << std::setw(6) << m_rot_angle1 << " "
	    << " RA2 = "    << std::setw(6) << m_rot_angle2 << " "
	    << " L = "      << std::setw(9) << m_length << " "
	    << " W0 = "     << std::setw(7) << m_w0 << " "
	    << " dWdX = "   << std::setw(6) << m_dd << " "
	    << " offset = " << std::setw(7) << m_offset << std::endl;
}
