// -*- C++ -*-

#include "BeamMan.hh"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>

#include <CLHEP/Units/SystemOfUnits.h>
#include <G4ThreeVector.hh>
#include <Randomize.hh>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "FuncName.hh"
#include "PrintHelper.hh"

//_____________________________________________________________________________
void
BeamInfo::FF(G4double ffplus)
{
  VO(1200*CLHEP::mm + ffplus);
}

//_____________________________________________________________________________
void
BeamInfo::VO(G4double voplus)
{
  if(mom.z() == 0.) return;
  pos.setX(x0 + u0*voplus);
  pos.setY(y0 + v0*voplus);
  pos.setZ(z0 + voplus);
}

//_____________________________________________________________________________
void
BeamInfo::Print() const
{
  PrintHelper helper(4, std::ios::fixed, G4cout);
  G4cout << "   "
         << "pos=" << pos << ",  "
         << "mom=" << mom << G4endl;
}

//_____________________________________________________________________________
BeamMan::BeamMan()
  : m_is_ready(false),
    m_file_name(),
    m_file(),
    m_param_array(),
    m_n_param()
{
}

//_____________________________________________________________________________
BeamMan::~BeamMan()
{
}

//_____________________________________________________________________________
G4bool
BeamMan::Initialize()
{
  const auto& confMan = ConfMan::GetInstance();
  // const auto& geomMan = DCGeomMan::GetInstance();
  const G4double p0 = confMan.Get<G4double>("PK18")*CLHEP::GeV;

  if(m_file_name.empty())
    return true;

  m_file = TFile::Open(m_file_name);
  if(!m_file || !m_file->IsOpen())
    return false;

  m_param_array.clear();

#if 1
  TTreeReader reader("bcout", m_file);
  TTreeReaderValue<Int_t> ntrack(reader, "ntrack");
  TTreeReaderArray<Double_t> chisqr(reader, "chisqr");
  TTreeReaderArray<Double_t> x0(reader, "x0");
  TTreeReaderArray<Double_t> y0(reader, "y0");
  TTreeReaderArray<Double_t> u0(reader, "u0");
  TTreeReaderArray<Double_t> v0(reader, "v0");
#else
  TTreeReader reader("k18track", m_file);
  TTreeReaderValue<Int_t> ntrack(reader, "ntK18");
  TTreeReaderArray<Double_t> chisqr(reader, "chisqrK18");
  TTreeReaderArray<Double_t> x0(reader, "xout");
  TTreeReaderArray<Double_t> y0(reader, "yout");
  TTreeReaderArray<Double_t> u0(reader, "uout");
  TTreeReaderArray<Double_t> v0(reader, "vout");
#endif

  while(reader.Next()){
    if(*ntrack != 1) continue;
    for(Int_t i=0, n=*ntrack; i<n; ++i){
      if(chisqr.At(i) > 5.) continue;
      BeamInfo beam;
      Double_t z0 = (-5326.5-1200.0)*CLHEP::mm; // VO
      beam.x0 = x0.At(i);
      beam.y0 = y0.At(i);
      beam.u0 = u0.At(i);
      beam.v0 = v0.At(i);
      beam.z0 = z0;
      beam.pos.set(x0.At(i), y0.At(i), z0);
      beam.mom.set(u0.At(i), v0.At(i), 1.);
      beam.mom.setMag(p0);
      m_param_array.push_back(beam);
    }
  }

  m_file->Close();
  m_n_param = m_param_array.size();
  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
G4bool
BeamMan::Initialize(const G4String& filename)
{
  m_file_name = filename;
  return Initialize();
}

//_____________________________________________________________________________
const BeamInfo&
BeamMan::Get() const
{
  if(m_is_ready && m_n_param > 0){
    return m_param_array.at(G4RandFlat::shootInt(m_n_param));
  }else{
    static BeamInfo nullinfo;
    return nullinfo;
  }
}

//_____________________________________________________________________________
void
BeamMan::Print() const
{
  PrintHelper helper(4, std::ios::fixed, G4cout);
  G4cout << FUNC_NAME << G4endl;
  for(const auto& b : m_param_array){
    G4cout << "   "
	   << "pos=" << b.pos << ",  "
	   << "mom=" << b.mom << G4endl;
  }
  G4cout << "   nparam = " << m_param_array.size() << G4endl;
}
