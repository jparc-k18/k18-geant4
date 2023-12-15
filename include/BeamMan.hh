// -*- C++ -*-

#ifndef BEAM_MAN_HH
#define BEAM_MAN_HH

#include <string>
#include <map>
#include <vector>

#include <globals.hh>
#include <G4ThreeVector.hh>

#include <Rtypes.h>

class TFile;

//_____________________________________________________________________________
struct BeamInfo
{
  G4double x0;
  G4double y0;
  G4double u0;
  G4double v0;
  G4double z0; // VO
  G4ThreeVector pos;
  G4ThreeVector mom;
  void FF(G4double ffplus=0);
  void VO(G4double voplus=0);
  void Print() const;
};

//_____________________________________________________________________________
class BeamMan
{
public:
  static G4String ClassName();
  static BeamMan& GetInstance();
  ~BeamMan();

private:
  BeamMan();
  BeamMan(const BeamMan& );
  BeamMan& operator =(const BeamMan&);

private:
  using ParamArray = std::vector<BeamInfo>;
  G4bool     m_is_ready;
  G4String   m_file_name;
  TFile*     m_file;
  ParamArray m_param_array;
  G4int      m_n_param;

public:
  const BeamInfo& Get() const;
  G4bool          Initialize();
  G4bool          Initialize(const G4String& filename);
  G4bool          IsReady() const { return m_is_ready; }
  void            Print() const;
};

//_____________________________________________________________________________
inline G4String
BeamMan::ClassName()
{
  static G4String s_name("BeamMan");
  return s_name;
}

//_____________________________________________________________________________
inline BeamMan&
BeamMan::GetInstance()
{
  static BeamMan s_instance;
  return s_instance;
}

#endif
