/*
  "ConfMan.hh"
  
  2007/4  K.Shirotori
  Modified by Toshiyuki Gogami , 21Nov2014
*/

#ifndef ConfMan_h
#define ConfMan_h 1

#include <string>
#include "G4String.hh"
#include "TMath.h"

class DCGeomMan;
//class EvDisp;

class ConfMan
{
public:
  explicit ConfMan( const std::string & filename );
  explicit ConfMan( const std::string & filename, const std::string & filename2 );
  ~ConfMan();
private:
  ConfMan( const ConfMan & );
  ConfMan & operator = ( const ConfMan & );

private:
  std::string ConfFileName_;
  static ConfMan *confManager_;

private:
  std::string DCGeomFileName_;
  DCGeomMan *DCGeomManager_;
  std::string BfieldMap_;
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
  bool Initialize( void );
  static ConfMan *GetConfManager( void ) { return confManager_; };
  std::string GetFieldMap( void ) { return BfieldMap_; };
  DCGeomMan *GetDCGeomManager( void ) { return DCGeomManager_; };
  //double K18Momentum( void ) const { return K18Momentum_;};
  int GeomFlag( void )    const  { return GeomFlag_;    };
  Int_t GetEMFlag(void)     const  { return EMFlag_;    };
  Int_t GetDecayFlag(void)  const  { return DecayFlag_; };
  Int_t GetHadronFlag(void)  const { return HadronFlag_;};
  G4String GetOutROOTFile()        { return oROOTFile;  };
  double GetGenMomCent() const     { return momcent;    };
  double GetGenMomBite() const     { return mombite;    };
  double GetGenTheta()   const     { return thetamax;   };
  Int_t  GetGenPID()     const     { return GenPID;     };
  double GetTOFDistance()const     { return tof_distance;};
  double GetTOFOverlap() const     { return tof_overlap; };
  double GetMagScale()   const     { return mag_scale;   };
  double GetMagScaleQ1()   const     { return mag_scale_Q1;   };
  double GetMagScaleQ2()   const     { return mag_scale_Q2;   };
  double GetTargetPositionZ()      { return tposz;       };
  Int_t  GetTargetID()             { return TargetID;    };
  Int_t  GetGenerator()             { return generator;    };
  double GetTThickness()           { return tthickness;  };
  double GetBeamX()                { return beamx;  };
  double GetBeamY()                { return beamy;  };
  double GetBeamZ()                { return beamz;  };
  void ShowParam();
  void OutputLog();
  void SetTargetPositionZ(double tpos_z) { target_pos = tpos_z; };
  double GetTargetPosizionZ(){ return target_pos; };

  //bool GetEvDispFlag( void ) const { return FlagEvDisp_; }
  //EvDisp *GetEvDisp( void ) { return evDisp_; }
  //bool InitializeEvDisp( void );


private:
  bool InitializeParameterFiles( void );
  bool EndAnalysis( void );

};
 
#endif
