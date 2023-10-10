// -*- C++ -*-

#include "DCGeomMan.hh"

#include <string>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>

#include <TString.h>

#include "DCGeomRecord.hh"
#include "FuncName.hh"

//_____________________________________________________________________________
DCGeomMan::DCGeomMan()
  : m_is_ready(false),
    m_file_name(),
    m_container(),
    m_detector_id_map(),
    m_global_z_map(),
    m_local_z_map()
{
}

//_____________________________________________________________________________
DCGeomMan::~DCGeomMan()
{
  Clear();
}

//_____________________________________________________________________________
G4int
DCGeomMan::CalcWireNumber(G4int lnum, G4double pos) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record){
    return record->WireNumber(pos);
  }
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
G4int
DCGeomMan::CalcWireNumber(const G4String& key, G4double wire) const
{
  return CalcWireNumber(GetDetectorId(key), wire);
}

//_____________________________________________________________________________
G4double
DCGeomMan::CalcWirePosition(G4int lnum, G4double wire) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record){
    return record->WirePos(wire);
  }
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
G4double
DCGeomMan::CalcWirePosition(const G4String& key, G4double wire) const
{
  return CalcWirePosition(GetDetectorId(key), wire);
}

//_____________________________________________________________________________
void
DCGeomMan::Clear()
{
  for(auto&& itr : m_container){
    delete itr.second;
  }
}

//_____________________________________________________________________________
G4double
DCGeomMan::GetLocalZ(G4int lnum) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record){
    return record->Length();
  }else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
G4double
DCGeomMan::GetLocalZ(const G4String& key) const
{
  return GetLocalZ(GetDetectorId(key));
}

//_____________________________________________________________________________
G4double
DCGeomMan::GetResolution(G4int lnum) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record) return record->Resolution();
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
G4double
DCGeomMan::GetResolution(const G4String& key) const
{
  return GetResolution(GetDetectorId(key));
}

//_____________________________________________________________________________
G4double
DCGeomMan::GetRotAngle1(G4int lnum) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record) return record->RotationAngle1();
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
G4double
DCGeomMan::GetRotAngle1(const G4String& key) const
{
  return GetRotAngle1(GetDetectorId(key));
}

//_____________________________________________________________________________
G4double
DCGeomMan::GetRotAngle2(G4int lnum) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record) return record->RotationAngle2();
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
G4double
DCGeomMan::GetRotAngle2(const G4String& key) const
{
  return GetRotAngle2(GetDetectorId(key));
}

//_____________________________________________________________________________
const G4ThreeVector&
DCGeomMan::GetGlobalPosition(G4int lnum) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record) return record->Pos();
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
const G4ThreeVector&
DCGeomMan::GetGlobalPosition(const G4String& key) const
{
  return GetGlobalPosition(GetDetectorId(key));
}

//_____________________________________________________________________________
G4double
DCGeomMan::GetTiltAngle(G4int lnum) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record) return record->TiltAngle();
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
G4double
DCGeomMan::GetTiltAngle(const G4String& key) const
{
  return GetTiltAngle(GetDetectorId(key));
}

//_____________________________________________________________________________
G4double
DCGeomMan::GetWirePitch(G4int lnum) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record) return record->WirePitch();
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
G4double
DCGeomMan::GetWirePitch(const G4String& key) const
{
  return GetWirePitch(GetDetectorId(key));
}

//_____________________________________________________________________________
G4bool
DCGeomMan::Initialize(const G4String& file_name)
{
  m_file_name = file_name;
  return Initialize();
}

//_____________________________________________________________________________
G4ThreeVector
DCGeomMan::NormalVector(G4int lnum) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record) return record->NormalVector();
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
G4ThreeVector
DCGeomMan::NormalVector(const G4String& key) const
{
  return NormalVector(GetDetectorId(key));
}

//_____________________________________________________________________________
G4ThreeVector
DCGeomMan::UnitVector(G4int lnum) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record) return record->UnitVector();
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
G4ThreeVector
DCGeomMan::UnitVector(const G4String& key) const
{
  return UnitVector(GetDetectorId(key));
}

//_____________________________________________________________________________
const DCGeomRecord*
DCGeomMan::GetRecord(G4int lnum) const
{
  DCGeomIterator itr = m_container.find(lnum);
  DCGeomIterator end = m_container.end();
  DCGeomRecord* record = 0;
  if(itr!=end)
    record = itr->second;
  if(!record){
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
  return record;
}

//_____________________________________________________________________________
const DCGeomRecord*
DCGeomMan::GetRecord(const G4String& key) const
{
  return GetRecord(GetDetectorId(key));
}

//_____________________________________________________________________________
G4bool
DCGeomMan::Initialize()
{
  if(m_is_ready){
    G4cerr << "#W " << FUNC_NAME
		<< " already initialied" << G4endl;
    return false;
  }

  std::ifstream ifs(m_file_name);
  if(!ifs.is_open()){
    G4cerr << FUNC_NAME << " file open fail : " << m_file_name << G4endl;
    return false;
  }

  Clear();

  TString line;
  while(ifs.good() && line.ReadLine(ifs)){
    if(line[0]=='#') continue;
    std::istringstream iss(line.Data());
    G4int id; G4String name;
    G4double gx, gy, gz, ta, ra1, ra2, l, res, w0, dd, ofs;
    if(iss >> id >> name >> gx >> gy >> gz >> ta >> ra1 >> ra2
	>> l >> res >> w0 >> dd >> ofs){
      DCGeomRecord* record =
	new DCGeomRecord(id, name, gx, gy, gz, ta, ra1, ra2,
			  l, res, w0, dd, ofs);
      if(m_container[id]){
	G4cerr << "#W " << FUNC_NAME << " "
		    << "duplicated key is deleted : " << id << G4endl;
	m_container[id]->Print();
	delete m_container[id];
      }
      m_container[id] = record;

      m_detector_id_map[name] = id;
      m_global_z_map[name]    = gz;
      m_local_z_map[name]     = l;
    }else{
      G4cerr << FUNC_NAME << " invalid format : " << line << G4endl;
    }
  }

  m_is_ready = true;
  return m_is_ready;
}

//_____________________________________________________________________________
std::vector<G4int>
DCGeomMan::GetDetectorIDList() const
{
  std::vector<G4int> vlist;
  vlist.reserve(m_container.size());
  DCGeomIterator itr, end=m_container.end();
  for(itr=m_container.begin(); itr!=end; ++itr){
    vlist.push_back(itr->first);
  }

  return vlist;
}

//_____________________________________________________________________________
G4ThreeVector
DCGeomMan::Local2GlobalPos(G4int lnum, const G4ThreeVector& in) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(!record){
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }

  G4double x = record->dxds()*in.x() + record->dxdt()*in.y()
    + record->dxdu()*in.z() + record->Pos().x();
  G4double y = record->dyds()*in.x() + record->dydt()*in.y()
    + record->dydu()*in.z() + record->Pos().y();
  G4double z = record->dzds()*in.x() + record->dzdt()*in.y()
    + record->dzdu()*in.z() + record->Pos().z();

  return G4ThreeVector(x, y, z);
}

//_____________________________________________________________________________
G4ThreeVector
DCGeomMan::Local2GlobalPos(const G4String& key, const G4ThreeVector& in) const
{
  return Local2GlobalPos(GetDetectorId(key), in);
}

//_____________________________________________________________________________
G4ThreeVector
DCGeomMan::Global2LocalPos(G4int lnum, const G4ThreeVector& in) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(!record){
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }

  G4double x
    = record->dsdx()*(in.x()-record->Pos().x())
    + record->dsdy()*(in.y()-record->Pos().y())
    + record->dsdz()*(in.z()-record->Pos().z());
  G4double y
    = record->dtdx()*(in.x()-record->Pos().x())
    + record->dtdy()*(in.y()-record->Pos().y())
    + record->dtdz()*(in.z()-record->Pos().z());
  G4double z
    = record->dudx()*(in.x()-record->Pos().x())
    + record->dudy()*(in.y()-record->Pos().y())
    + record->dudz()*(in.z()-record->Pos().z());

  return G4ThreeVector(x, y, z);
}

//_____________________________________________________________________________
G4ThreeVector
DCGeomMan::Global2LocalPos(const G4String& key, const G4ThreeVector& in) const
{
  return Global2LocalPos(GetDetectorId(key), in);
}

//_____________________________________________________________________________
G4ThreeVector
DCGeomMan::Local2GlobalDir(G4int lnum, const G4ThreeVector& in) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(!record){
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }

  G4double x = record->dxds()*in.x() + record->dxdt()*in.y()
    + record->dxdu()*in.z();
  G4double y = record->dyds()*in.x() + record->dydt()*in.y()
    + record->dydu()*in.z();
  G4double z = record->dzds()*in.x() + record->dzdt()*in.y()
    + record->dzdu()*in.z();

  return G4ThreeVector(x, y, z);
}

//_____________________________________________________________________________
G4ThreeVector
DCGeomMan::Local2GlobalDir(const G4String& key, const G4ThreeVector& in) const
{
  return Local2GlobalDir(GetDetectorId(key), in);
}

//_____________________________________________________________________________
G4ThreeVector
DCGeomMan::Global2LocalDir(G4int lnum, const G4ThreeVector& in) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(!record){
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }

  G4double x = record->dsdx()*in.x() + record->dsdy()*in.y()
    + record->dsdz()*in.z();
  G4double y = record->dtdx()*in.x() + record->dtdy()*in.y()
    + record->dtdz()*in.z();
  G4double z = record->dudx()*in.x() + record->dudy()*in.y()
    + record->dudz()*in.z();

  return G4ThreeVector(x, y, z);
}

//_____________________________________________________________________________
G4ThreeVector
DCGeomMan::Global2LocalDir(const G4String& key, const G4ThreeVector& in) const
{
  return Global2LocalDir(GetDetectorId(key), in);
}

//_____________________________________________________________________________
G4int
DCGeomMan::GetDetectorId(const G4String& key) const
{
  DCGeomIterator itr, end = m_container.end();
  for(itr=m_container.begin(); itr!=end; ++itr){
    if (itr->second->Name() == key)
      return itr->second->Id();
  }
  G4String e(FUNC_NAME + Form(" No record : %s", key.data()));
  throw std::invalid_argument(e);
}

//_____________________________________________________________________________
void
DCGeomMan::SetFileName(const G4String& file_name)
{
  m_file_name = file_name;
}

//_____________________________________________________________________________
void
DCGeomMan::SetResolution(G4int lnum, G4double res)
{
  DCGeomRecord* record = m_container[lnum];
  if(record) {
    record->SetResolution(res);
  }
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
void
DCGeomMan::SetResolution(const G4String& key, G4double res)
{
  SetResolution(GetDetectorId(key), res);
}
