/*
  DCConfMan.cc
*/

#include "DCConfMan.hh"
#include "DCGeomRecord.hh"

#include <string>
#include <stdexcept>
#include <cstdio>
#include <cstring>

const int MaxChar = 200;

DCConfMan *DCConfMan::confMan_=0;

DCConfMan::DCConfMan()
  : TOFid_(61)
{}

DCConfMan::~DCConfMan()
{}


DCConfMan & DCConfMan::GetInstance( void )
{
  if( !confMan_ ){
    confMan_ = new DCConfMan();
  }
  return *confMan_;
}

double DCConfMan::GetLocalZ( int lnum ) const
{
  static const std::string funcname = "[DCConfMan::GetLocalZ(int)]"; 
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->length_;
  else{
    std::cerr << funcname << ": No record. Layer yade #=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

double DCConfMan::GetResolution( int lnum ) const
{
  static const std::string funcname = "[DCConfMan::GetResolution(int)]"; 
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->resol_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

double DCConfMan::GetTiltAngle( int lnum ) const
{
  static const std::string funcname = "[DCConfMan::GetTiltAngle(int)]";
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->tiltAngle_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

double DCConfMan::GetRotAngle1( int lnum ) const
{
  static const std::string funcname = "[DCConfMan::GetRotAngle1(int)]";
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->rotAngle1_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

double DCConfMan::GetRotAngle2( int lnum ) const
{
  static const std::string funcname = "[DCConfMan::GetRotAngle2(int)]";
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->rotAngle2_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

const ThreeVector & DCConfMan::GetGlobalPosition( int lnum ) const
{
  static const std::string funcname = "[DCConfMan::GetGlobalPosition(int)]";
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->pos_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

ThreeVector DCConfMan::NormalVector( int lnum ) const
{
  static const std::string funcname = "[DCConfMan::NormalVector(int)]";
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->NormalVector();
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

ThreeVector DCConfMan::UnitVector( int lnum ) const
{
  static const std::string funcname = "[DCConfMan::UnitVector(int)]";
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo->UnitVector();
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

const DCGeomRecord * DCConfMan::GetRecord( int lnum ) const
{
  static const std::string funcname = "[DCConfMan::GetRecord(int)]";
  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( pGeo ) return pGeo;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}
 
  

void DCConfMan::clearElements( void )
{
  //  for_each( geomRecord_.begin(), geomRecord_.end(), DeleteObject() );
  std::map <int, DCGeomRecord *>::iterator itr;
  for( itr=geomRecord_.begin(); itr!=geomRecord_.end(); ++itr ){
    delete itr->second;
  }
  geomRecord_.clear();
  TOFid_=61;
}

bool DCConfMan::Initialize( void )
{
  static const std::string funcname = "[DCConfMan::Initialize]";
  char str[MaxChar];
  char cname[MaxChar];
  int id;
  double xs, ys, zs, ta, ra1, ra2, l, res;

  FILE *fp;

  if( ( fp = fopen( filename_.c_str(), "r" ) ) == 0 ){
    throw std::invalid_argument(funcname+": file open fail");
  }

  clearElements();

  while( fgets( str, MaxChar, fp ) != 0 ){
    if( str[0]!='#' ){
      if( sscanf( str, "%d %s %lf %lf %lf %lf %lf %lf %lf %lf",
		  &id, cname, &xs, &ys, &zs, &ta, &ra1, &ra2, &l, &res )
	  == 10 ){
	DCGeomRecord *pRec = 
	  new DCGeomRecord( id, cname, xs, ys, zs, ta, ra1, ra2, l, res );
	DCGeomRecord *pOld = geomRecord_[id];
	geomRecord_[id] = pRec;

	if(strcmp(cname,"TOF")==0) {
	  TOFid_=id;
	}

	if( pOld ){
	  std::cerr << funcname << ": duplicated id number. "
		    << " following record is deleted." << std::endl;
	  std::cerr << "Id=" << pOld->id_ << " " << pOld->pos_
		    << " ) ... " << std::endl;
	  delete pOld;
	}
      }
      else {
	std::string strtemp=str;
	std::cerr << funcname << ": Invalid format " << strtemp << std::endl;
      }
    }
  }

  fclose(fp);

  return true;
}


std::vector <int> DCConfMan::GetDetectorIDList( void ) const
{
  std::vector<int> vlist;
  vlist.reserve(geomRecord_.size());
  std::map <int, DCGeomRecord *>::const_iterator 
    itr=geomRecord_.begin(), end=geomRecord_.end();

  for(; itr!=end; ++itr ){
    vlist.push_back( itr->first );
  }

  return vlist;
}

ThreeVector DCConfMan::Local2GlobalPos( int lnum, 
					const ThreeVector &in ) const
{
  static const std::string funcname = 
    "[DCConfMan::Local2GlobalPos(ThreeVecor &)]";

  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( !pGeo ) 
    throw std::out_of_range(funcname+": No record" );

  double x = pGeo->dxds_*in.x() + pGeo->dxdt_*in.y()
    + pGeo->dxdu_*in.z() + pGeo->pos_.x();
  double y = pGeo->dyds_*in.x() + pGeo->dydt_*in.y()
    + pGeo->dydu_*in.z() + pGeo->pos_.y();
  double z = pGeo->dzds_*in.x() + pGeo->dzdt_*in.y()
    + pGeo->dzdu_*in.z() + pGeo->pos_.z();

  return ThreeVector( x, y, z );
}

ThreeVector DCConfMan::Global2LocalPos( int lnum,
					const ThreeVector &in ) const
{
  static const std::string funcname = 
    "[DCConfMan::Global2LocalPos(ThreeVecor &)]";

  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( !pGeo ){
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }

  double x 
    = pGeo->dsdx_*(in.x()-pGeo->pos_.x())
    + pGeo->dsdy_*(in.y()-pGeo->pos_.y())
    + pGeo->dsdz_*(in.z()-pGeo->pos_.z());
  double y 
    = pGeo->dtdx_*(in.x()-pGeo->pos_.x())
    + pGeo->dtdy_*(in.y()-pGeo->pos_.y())
    + pGeo->dtdz_*(in.z()-pGeo->pos_.z());
  double z 
    = pGeo->dudx_*(in.x()-pGeo->pos_.x())
    + pGeo->dudy_*(in.y()-pGeo->pos_.y())
    + pGeo->dudz_*(in.z()-pGeo->pos_.z());

  return ThreeVector( x, y, z );
}

ThreeVector DCConfMan::Local2GlobalDir( int lnum,
					const ThreeVector &in ) const
{
  static const std::string funcname = 
    "[DCConfMan::Local2GlobalDir(ThreeVecor &)]";

  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( !pGeo ){
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }

  double x = pGeo->dxds_*in.x() + pGeo->dxdt_*in.y()
    + pGeo->dxdu_*in.z();
  double y = pGeo->dyds_*in.x() + pGeo->dydt_*in.y()
    + pGeo->dydu_*in.z();
  double z = pGeo->dzds_*in.x() + pGeo->dzdt_*in.y()
    + pGeo->dzdu_*in.z();

  return ThreeVector( x, y, z );
}

ThreeVector DCConfMan::Global2LocalDir( int lnum,
					const ThreeVector &in ) const
{
  static const std::string funcname = 
    "[DCConfMan::Global2LocalDir(ThreeVecor &)]";

  DCGeomRecord *pGeo = geomRecord_[lnum];
  if( !pGeo ){
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }

  double x = pGeo->dsdx_*in.x() + pGeo->dsdy_*in.y()+ pGeo->dsdz_*in.z();
  double y = pGeo->dtdx_*in.x() + pGeo->dtdy_*in.y()+ pGeo->dtdz_*in.z();
  double z = pGeo->dudx_*in.x() + pGeo->dudy_*in.y()+ pGeo->dudz_*in.z();

  return ThreeVector( x, y, z );
}
