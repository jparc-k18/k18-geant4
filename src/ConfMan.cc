/*
  "ConfMan.cc"
  
  2007/4  K.Shirotori
  Modified by Toshi Gogami , 21Nov2014
*/

#include "ConfMan.hh"
#include "DCGeomMan.hh"
//#include "EvDisp.hh"
#include "TString.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
using namespace std;

ConfMan * ConfMan::confManager_ = 0;

ConfMan::ConfMan( const std::string & filename )
  : ConfFileName_(filename),DCGeomManager_(0),// K18Momentum_(1.8),
    GeomFlag_(0), EMFlag_(0), DecayFlag_(0),
    momcent(0.0), mombite(0.0),
    tof_overlap(0.0),tof_distance(-20.0),
    mag_scale(1.00), mag_scale_Q1(1.00), mag_scale_Q2(1.00),
    generator(0), oROOTFile("0")
{
  static const std::string funcname = "[ConfMan::ConfMan]";
  if( confManager_ ){
    std::cerr << funcname << ": constructing twice" << std::endl;
    exit(-1);
  }
  confManager_ = this;
}

ConfMan::ConfMan( const std::string & filename, const std::string & filename2)
  : ConfFileName_(filename),DCGeomManager_(0),// K18Momentum_(1.8),
    GeomFlag_(0), EMFlag_(0), DecayFlag_(0),
    momcent(0.0), mombite(0.0),
    tof_overlap(0.0),tof_distance(-20.0),
    mag_scale(1.00), mag_scale_Q1(1.00), mag_scale_Q2(1.00),
    generator(0), oROOTFile(filename2)
{
  static const std::string funcname = "[ConfMan::ConfMan]";
  if( confManager_ ){
    std::cerr << funcname << ": constructing twice" << std::endl;
    exit(-1);
  }
  confManager_ = this;
}

ConfMan::~ConfMan()
{
  EndAnalysis();
  confManager_=0;
}

bool ConfMan::EndAnalysis( void )
{

  if(DCGeomManager_){
    delete DCGeomManager_; DCGeomManager_=0;
  }

  return true;
}


const int BufSize = 144;

bool ConfMan::Initialize( void )
{
  static const std::string funcname = "[ConfMan::Initialize]";

  FILE *fp;
  char buf[BufSize], buf1[BufSize], buf2[BufSize+1];
  int    intval;
  double val;

  if((fp=fopen(ConfFileName_.c_str(),"r"))==0){
    std::cerr << funcname << ": file open fail" << std::endl;
    exit(-1);
  }
  bool mag_Q1 =false;
  bool mag_Q2 =false;
  // ~~~~~~~~~ Read configuration file ~~~~~~~~~~~~~~~~~~~~~~
  while( fgets(buf,BufSize,fp)!=0 ){
    if( buf[0]!='#' ){
      //Geometry
      if( sscanf(buf,"DCGEO: %s",buf1)==1 ){
	DCGeomFileName_=buf1;
      }
      //Primary generator
	
      else if( sscanf(buf,"BMAP: %s",buf1)==1 ){
	BfieldMap_=buf1;
        link_len = readlink(buf1, buf2, BufSize); // read symboric link
        if (link_len > 0) BfieldMap_link= Form(" -> %s", buf2);
	else BfieldMap_link = "";
      }
      else if( sscanf(buf,"GenMomCent[GeV/c]: %lf", &val )==1 ){
	momcent=val;
      }
      else if( sscanf(buf,"GenMomBite[GeV/c]: %lf", &val )==1 ){
	mombite=val;
      }
      else if( sscanf(buf,"GenTheta[deg]: %lf", &val )==1 ){
	thetamax=val;
      }
      else if( sscanf(buf,"TPOSZ[mm]: %lf", &val )==1 ){
	tposz = val;
      }
      else if( sscanf(buf,"TARGET: %d", &intval )==1 ){
	TargetID = intval;
      }
      else if( sscanf(buf,"TTHICKNESS[g/cm2]: %lf", &val )==1 ){
	tthickness = val;
      }
      else if( sscanf(buf,"GenPID: %d", &intval )==1 ){
	GenPID=intval;
      }
      else if( sscanf(buf,"BeamWidth[mm]: %lf %lf %lf", 
		      &beamx, &beamy, &beamz  )==3 ){}
      //else if( sscanf(buf,"PK18: %lf", &val )==1 ){
      //	K18Momentum_=val;
      //}
      else if( sscanf(buf,"EM: %d", &intval )==1 ){
	EMFlag_=intval;
      }
      else if( sscanf(buf,"DECAY: %d", &intval )==1 ){
	DecayFlag_=intval;
      }
      else if( sscanf(buf,"HADRON: %d", &intval )==1 ){
	HadronFlag_=intval;
      }
      else if( sscanf(buf,"ROOTFile: %s",buf1)==1 ){
	if(oROOTFile[0]=='0')
	oROOTFile = buf1;
      }
      else if( sscanf(buf,"TOF_OVERLAP[mm]: %lf", &val )==1 ){
	tof_overlap=val;
      }
      else if( sscanf(buf,"TOF_DISTANCE[mm]: %lf", &val )==1 ){
	tof_distance=val;
      }
      else if( sscanf(buf,"Mag_Scale: %lf", &val )==1 ){
	mag_scale=val;
      }
      else if( sscanf(buf,"Mag_Scale_Q1: %lf", &val )==1 ){
	mag_scale_Q1=val;
	mag_Q1 = true;
      }
      else if( sscanf(buf,"Mag_Scale_Q2: %lf", &val )==1 ){
	mag_scale_Q2=val;
	mag_Q2 = true;
      }
      else if( sscanf(buf,"Generator: %lf", &val )==1 ){
	generator=val;
      }
      
      // Event display
      //else if( sscanf(buf,"EVDISP: %d",&intval)==1 )
      //if(intval==1) FlagEvDisp_=true;
      //else         FlagEvDisp_=false;

    } /* if( buf[0]... ) */
  } /* while(...) */

  if(!mag_Q1)
    mag_scale_Q1 = mag_scale;
  if(!mag_Q2)
    mag_scale_Q2 = mag_scale;
  
  fclose(fp);

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
  InitializeParameterFiles();

  return true;
}

bool ConfMan::InitializeParameterFiles( void )
{
  DCGeomManager_ = & DCGeomMan::GetInstance();
  
  if( DCGeomFileName_!="" ){
    DCGeomManager_->Initialize(DCGeomFileName_);
  }
  else {
    DCGeomManager_->Initialize();
  }
  
  return true;
}

void ConfMan::ShowParam(){
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
	 << "(uniform) mm" << endl;
  G4cout << "EMFlag:      " << EMFlag_      << G4endl;
  G4cout << "DecayFlag:   " << DecayFlag_   << G4endl;
  G4cout << "HadronFlag:  " << HadronFlag_  << G4endl;
  G4cout << "Generator:  " << generator  << G4endl;
  G4cout << "                 ╭( ･ㅂ･)و ̑̑" << G4endl;
}

void ConfMan::OutputLog(){
  
  std::string LogFoot("_Log");
  std::string LogFileName = oROOTFile+LogFoot;
  
  ofstream* ofs = new ofstream(LogFileName);
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
       << beamz << "(uniform) mm" << endl;
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
