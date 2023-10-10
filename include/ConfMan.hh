// -*- C++ -*-

#ifndef ConfMan_h
#define ConfMan_h 1

#include <iomanip>
#include <map>
#include <string>

#include <G4String.hh>

#include <TMath.h>

//class EvDisp;

//_____________________________________________________________________________
class ConfMan
{
public:
  static G4String ClassName();
  static ConfMan& GetInstance();
  ~ConfMan();
private:
  ConfMan();
  ConfMan(const ConfMan&);
  ConfMan& operator =(const ConfMan&);

private:
  using StrList = std::map<G4String, G4String>;
  using DoubleList = std::map<G4String, G4double>;
  using IntList = std::map<G4String, G4int>;
  using BoolList = std::map<G4String, G4bool>;
  std::string ConfFileName_;
  G4String   m_conf_key;
  G4String   m_conf_dir;
  G4String   m_conf_buf;
  StrList    m_file;
  StrList    m_string;
  DoubleList m_double;
  IntList    m_int;
  BoolList   m_bool;


  std::string DCGeomFileName_;
  std::string BfieldMap_;
  std::string BfieldMap_link; // symbolic link for Bfield
  std::string oROOTFile;
  //std::string LogFileName;
  //G4String oROOTFile;
  //EvDisp *evDisp_;
  //bool FlagEvDisp_;
  //double K18Momentum_;
  int GeomFlag_;
  Int_t EMFlag_;
  Int_t DecayFlag_;
  Int_t HadronFlag_;
  Int_t GenPID;
  Int_t TargetID;
  Int_t generator;
  Int_t link_len; // symbolic link for Bfield
  double tposz;
  double tthickness;
  double momcent, mombite;
  double thetamax;
  double tof_overlap;
  double tof_distance;
  double mag_scale, mag_scale_Q1, mag_scale_Q2;
  double target_pos;
  double beamx, beamy, beamz;

public:
  G4bool Initialize(const G4String& file_name);
  const G4String& ConfBuf() const { return m_conf_buf; }
  template <typename T>
  static const T& Get(const G4String& key);
  std::string GetFieldMap( void ) const { return BfieldMap_; };
  //double K18Momentum( void ) const { return K18Momentum_;};
  int GeomFlag( void )    const  { return GeomFlag_;    };
  Int_t GetEMFlag(void)     const  { return EMFlag_;    };
  Int_t GetDecayFlag(void)  const  { return DecayFlag_; };
  Int_t GetHadronFlag(void)  const { return HadronFlag_;};
  G4String GetOutROOTFile() const  { return oROOTFile;  };
  double GetGenMomCent() const     { return momcent;    };
  double GetGenMomBite() const     { return mombite;    };
  double GetGenTheta()   const     { return thetamax;   };
  Int_t  GetGenPID()     const     { return GenPID;     };
  double GetTOFDistance()const     { return tof_distance;};
  double GetTOFOverlap() const     { return tof_overlap; };
  double GetMagScale()   const     { return mag_scale;   };
  double GetMagScaleQ1()   const   { return mag_scale_Q1;   };
  double GetMagScaleQ2()   const   { return mag_scale_Q2;   };
  double GetTargetPositionZ() const { return tposz;       };
  Int_t  GetTargetID() const       { return TargetID;    };
  Int_t  GetGenerator() const      { return generator;    };
  double GetTThickness() const     { return tthickness;  };
  double GetBeamX() const          { return beamx;  };
  double GetBeamY() const          { return beamy;  };
  double GetBeamZ() const          { return beamz;  };
  void ShowParam();
  void OutputLog();
  void SetTargetPositionZ(double tpos_z) { target_pos = tpos_z; };
  double GetTargetPosizionZ(){ return target_pos; };

  //bool GetEvDispFlag( void ) const { return FlagEvDisp_; }
  //EvDisp *GetEvDisp( void ) { return evDisp_; }
  //bool InitializeEvDisp( void );


private:
  G4String FilePath(const G4String& src) const;
  G4bool   InitializeParameterFiles();
  template <typename T>
  G4bool    InitializeParameter();
  template <typename T>
  G4bool    InitializeParameter(const G4String& key);
  template <typename T>
  G4bool    InitializeParameter(const G4String& key1,
                                const G4String& key2);
  G4bool   ShowResult(G4bool status, const G4String& name) const;
};

//_____________________________________________________________________________
inline G4String
ConfMan::ClassName()
{
  static const G4String s_name("ConfMan");
  return s_name;
}

//_____________________________________________________________________________
inline ConfMan&
ConfMan::GetInstance()
{
  static ConfMan s_instance;
  return s_instance;
}

//_____________________________________________________________________________
template <>
inline const G4String&
ConfMan::Get<G4String>(const G4String& key)
{
  return GetInstance().m_string[key];
}

//_____________________________________________________________________________
template <>
inline const G4double&
ConfMan::Get<G4double>(const G4String& key)
{
  return GetInstance().m_double[key];
}

//_____________________________________________________________________________
template <>
inline const G4int&
ConfMan::Get<G4int>(const G4String& key)
{
  return GetInstance().m_int[key];
}

//_____________________________________________________________________________
template <>
inline const G4bool&
ConfMan::Get<G4bool>(const G4String& key)
{
  return GetInstance().m_bool[key];
}

//_____________________________________________________________________________
inline G4bool
ConfMan::ShowResult(G4bool status, const G4String& name) const
{
  if(status)
    std::cout << std::setw(20) << std::left
	      << " ["+name+"]"
	      << "-> Initialized" << std::endl;
  else
    std::cout << std::setw(20) << std::left
	      << " ["+name+"]"
	      << "-> Failed" << std::endl;
  return status;
}

//_____________________________________________________________________________
template <typename T>
inline G4bool
ConfMan::InitializeParameter()
{
  return
    ShowResult(T::GetInstance().Initialize(),
               T::GetInstance().ClassName());
}

//_____________________________________________________________________________
template <typename T>
inline G4bool
ConfMan::InitializeParameter(const G4String& key)
{
  return
    ShowResult(T::GetInstance().Initialize(m_file[key]),
               T::GetInstance().ClassName());
}

//_____________________________________________________________________________
template <typename T>
inline G4bool
ConfMan::InitializeParameter(const G4String& key1,
                             const G4String& key2)
{
  return
    ShowResult(T::GetInstance().Initialize(m_file[key1],
                                           m_file[key2]),
               T::GetInstance().ClassName());
}

#endif
