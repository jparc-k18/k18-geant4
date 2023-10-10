// -*- C++ -*-

#include "ConfMan.hh"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <libgen.h>
#include <unistd.h>

#include <TString.h>

#include "DCGeomMan.hh"
#include "FuncName.hh"
#include "HistMan.hh"

ConfMan::ConfMan()
  : ConfFileName_(),
    m_conf_key("CONF"),
    m_conf_dir(),
    m_conf_buf(),
    // K18Momentum_(1.8),
    oROOTFile("0"),
    GeomFlag_(0),
    EMFlag_(0),
    DecayFlag_(0),
    generator(0),
    momcent(0.0),
    mombite(0.0),
    tof_overlap(0.0),
    tof_distance(-20.0),
    mag_scale(1.00),
    mag_scale_Q1(1.00),
    mag_scale_Q2(1.00)
{
}

ConfMan::~ConfMan()
{
  ShowParam();
  OutputLog();
}

const int BufSize = 144;

G4bool ConfMan::Initialize(const G4String& file_name)
{
  static const std::string funcname = "[ConfMan::Initialize]";

  ConfFileName_ = file_name;
  m_file[m_conf_key] = file_name;
  m_conf_dir = ::dirname(const_cast<char*>(m_file[m_conf_key].data()));
  m_conf_buf.clear();
  m_conf_buf += "\n";

  char buf[BufSize], buf1[BufSize], buf2[BufSize+1];
  int    intval;

  std::ifstream ifs(ConfFileName_);
  if(!ifs.is_open()){
    G4cerr << FUNC_NAME << " file open fail : " << ConfFileName_ << G4endl;
    return false;
  }

  bool mag_Q1 =false;
  bool mag_Q2 =false;
  // ~~~~~~~~~ Read configuration file ~~~~~~~~~~~~~~~~~~~~~~

  G4String line;
  while(ifs.good() && std::getline(ifs, line)){
    if(line.empty() || line[0]=='#') continue;
    m_conf_buf += line + "\n";
    std::istringstream iss(line);
    G4String key, val;
    iss >> key >> val;
    if(key.empty() || val.empty())
      continue;
    if(key.back() == ':')
      key.pop_back();
    G4cout << " key = "   << std::setw(20) << std::left << key
	   << " value = " << std::setw(30) << std::left << val
	   << G4endl;

    m_file[key]   = FilePath(val);
    m_string[key] = val;
    m_double[key] = std::strtod(val, nullptr);
    m_int[key]    = std::strtol(val, nullptr, 10);
    m_bool[key]   = static_cast<G4bool>(std::strtol(val, nullptr, 10));

    if(key == "DCGEO"){
      DCGeomFileName_ = val;
    }
    else if(key == "BMAP"){
      BfieldMap_ = val;
    }
    else if(key == "GenMomCent[GeV/c]"){
      momcent = std::strtod(val, nullptr);
    }
    else if(key == "GenMomBite[GeV/c]"){
      mombite = std::strtod(val, nullptr);
    }
    else if(key == "GenTheta[deg]"){
      thetamax = std::strtod(val, nullptr);
    }
    else if(key == "TPOSZ[mm]"){
      tposz = std::strtod(val, nullptr);
    }
    else if(key == "TARGET"){
      TargetID = std::strtod(val, nullptr);
    }
    else if(key == "TTHICKNESS[g/cm2]"){
      tthickness = std::strtod(val, nullptr);
    }
    else if(key == "GenPID"){
      GenPID = std::strtod(val, nullptr);
    }
    // else if( sscanf(buf,"BeamWidth[mm]: %lf %lf %lf",
    //                 &beamx, &beamy, &beamz  )==3 ){}
    //else if( sscanf(buf,"PK18: %lf", &val )==1 ){
    //	K18Momentum_=val;
    //}
    else if(key == "EM"){
      EMFlag_ = std::strtod(val, nullptr);
    }
    // else if( sscanf(buf,"DECAY: %d", &intval )==1 ){
    //   DecayFlag_=intval;
    // }
    // else if( sscanf(buf,"HADRON: %d", &intval )==1 ){
    //   HadronFlag_=intval;
    // }
    // else if( sscanf(buf,"ROOTFile: %s",buf1)==1 ){
    //   if(oROOTFile[0]=='0')
    //     oROOTFile = buf1;
    // }
    // else if( sscanf(buf,"TOF_OVERLAP[mm]: %lf", &val )==1 ){
    //   tof_overlap=val;
    // }
    // else if( sscanf(buf,"TOF_DISTANCE[mm]: %lf", &val )==1 ){
    //   tof_distance=val;
    // }
    else if(key == "Mag_Scale"){
      mag_scale = std::strtod(val, nullptr);
    }
    else if(key == "Mag_Scale_Q1"){
      mag_scale_Q1 = std::strtod(val, nullptr);
      mag_Q1 = true;
    }
    else if(key == "Mag_Scale_Q2"){
      mag_scale_Q2 = std::strtod(val, nullptr);
      mag_Q2 = true;
    }
    else if(key == "Generator"){
      generator = std::strtod(val, nullptr);
    }

    // // Event display
    // //else if( sscanf(buf,"EVDISP: %d",&intval)==1 )
    // //if(intval==1) FlagEvDisp_=true;
    // //else         FlagEvDisp_=false;
  }

  if(!mag_Q1)
    mag_scale_Q1 = mag_scale;
  if(!mag_Q2)
    mag_scale_Q2 = mag_scale;

  /*
    std::cout << "----------"  << ConfFileName_ << "--------" << std::endl;
    std::cout << "**********Geometry**********" << std::endl;
    std::cout << "DC Geom. Param.:  "  << DCGeomFileName_   << std::endl;
    std::cout << "**********BfieldMap*********" << std::endl;
    std::cout << "B field Map :     "  << BfieldMap_  << std::endl;
    std::cout << "**********Primary Action**********" << std::endl;
    std::cout << "Momentum [GeV/c]:    "  << momcent << " +/- "
    << mombite << std::endl;
    std::cout << "**********Physics Process**********" << std::endl;
    std::cout << "EM Process:    "  << EMFlag_              << std::endl;
    std::cout << "Decay Process:    "  << DecayFlag_              << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
  */

  return InitializeParameterFiles();
}

bool ConfMan::InitializeParameterFiles( void )
{
  return (true
          && InitializeParameter<DCGeomMan>("DCGEO")
          // && InitializeParameter<BeamMan>("BEAM")
          // && InitializeParameter<DetSizeMan>("DSIZE")
          && InitializeParameter<HistMan>("HIST")
          // && InitializeParameter<JamMan>("JAM")
          // && InitializeParameter<IncMan>("INC")
          );
}

void ConfMan::ShowParam()
{
  G4cout << G4endl;
  G4cout << " ------------ Used Parameters ------------ "  << G4endl;
  G4cout << "ROOT file:   " << oROOTFile    << G4endl;
  G4cout << "DCGeometry:  " << DCGeomFileName_ << G4endl;
  G4cout << "FieldMap:    " << BfieldMap_  << BfieldMap_link  << G4endl;
  G4cout << "Mag Scale:   " << mag_scale   << G4endl;
  G4cout << "Mag Scale Q1:   " << mag_scale_Q1   << G4endl;
  G4cout << "Mag Scale Q2:   " << mag_scale_Q2   << G4endl;
  //G4cout << "Momentum:    " << K18Momentum_ << G4endl;
  G4cout << "TargetID:    " << TargetID << G4endl;
  G4cout << "TargetThickness: " << tthickness  << " g/cm^{2}"<< G4endl;
  G4cout << "TargetPosZ:  " << tposz  << " mm"<< G4endl;
  G4cout << "Momentum:    " << momcent << " +/- "
	 << mombite << " GeV/c" << G4endl;
  G4cout << "Theta:       0 - " << thetamax << " deg " <<G4endl;
  G4cout << "GenPID:      " << GenPID
	 << " (1:K+ 2:K- 3:pi+ 4:pi- 5:p 6:e- 7:mu- 8:xi-)" << G4endl;
  G4cout << "BeamWidth:   " << beamx << "(sigma), " << beamy << "(sigma), " << beamz
	 << "(uniform) mm" << G4endl;
  G4cout << "EMFlag:      " << EMFlag_      << G4endl;
  G4cout << "DecayFlag:   " << DecayFlag_   << G4endl;
  G4cout << "HadronFlag:  " << HadronFlag_  << G4endl;
  G4cout << "Generator:  " << generator  << G4endl;
  G4cout << "                 ╭( ･ㅂ･)و ̑̑" << G4endl;
}

void ConfMan::OutputLog(){

  std::string LogFoot("_Log");
  std::string LogFileName = oROOTFile+LogFoot;

  auto ofs = new std::ofstream(LogFileName);
  //*ofs << G4endl;
  *ofs << "  /// S-2S Geant4 simulation ///" << G4endl;
  *ofs << "  /// Used Parameters -->    ///"  << G4endl;
  *ofs << "ROOT file:   " << oROOTFile    << G4endl;
  *ofs << "DCGeometry:  " << DCGeomFileName_ << G4endl;
  *ofs << "FieldMap:    " << BfieldMap_  << BfieldMap_link  << G4endl;
  *ofs << "Mag Scale:   " << mag_scale   << G4endl;
  *ofs << "Mag Scale Q1:   " << mag_scale_Q1   << G4endl;
  *ofs << "Mag Scale Q2:   " << mag_scale_Q2   << G4endl;
  *ofs << "TargetID:    " << TargetID << G4endl;
  *ofs << "TargetThickness: " << tthickness << " g/cm^{2}" << G4endl;
  *ofs << "TargetPosZ:  " << tposz << " mm" << G4endl;
  //G4cout << "Momentum:    " << K18Momentum_ << G4endl;
  *ofs << "Momentum:    " << momcent << " +/- "
       << mombite << " GeV/c" << G4endl;
  *ofs << "Theta:       0 - " << thetamax << " deg " <<G4endl;
  *ofs << "GenPID:      " << GenPID
       << " (1:K+ 2:K- 3:pi+ 4:pi- 5:p 6:e- 7:mu- 8:xi-)" << G4endl;
  *ofs << "BeamWidth:   " << beamx << "(sigma), " << beamy << "(sigma), "
       << beamz << "(uniform) mm" << G4endl;
  *ofs << "EMFlag:      " << EMFlag_      << G4endl;
  *ofs << "DecayFlag:   " << DecayFlag_   << G4endl;
  *ofs << "HadronFlag:  " << HadronFlag_  << G4endl;
  *ofs << "Generator:  " << generator  << G4endl;
  //*ofs << "                 ╭( ･ㅂ･)و ̑̑" << G4endl;
  ofs->close();
}

/*
  bool ConfMan::InitializeEvDisp( void )
  {
  static const std::string funcname = "[ConfMan::InitializeEvDisp]";
  evDisp_ = & EvDisp::GetInstance();
  evDisp_->Initialize();
  }
*/

//_____________________________________________________________________________
G4String
ConfMan::FilePath(const G4String& src) const
{
  std::ifstream tmp(src);
  if (tmp.good())
    return src;
  else
    return m_conf_dir + "/" + src;
}
