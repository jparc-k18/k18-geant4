/*
  DCConfMan.hh
*/

#ifndef DCConfMan_h
#define DCConfMan_h 1

#include "ThreeVector.hh"
#include <string>
#include <vector>
#include <map>

class DCGeomRecord;

class DCConfMan 
{
private:
  DCConfMan();
public:
  ~DCConfMan();

public:
  void SetFileName( const char *filename ) { filename_=filename; }
  void SetFileName( const std::string &filename ) { filename_=filename; }

  bool Initialize( void );
  bool Initialize( const char *filename )
  { filename_=filename; Initialize(); }
  bool Initialize( const std::string &filename )
  { filename_=filename; Initialize(); }

  static DCConfMan & GetInstance( void );
  double GetLocalZ( int lnum ) const;
  double GetResolution( int lnum ) const;
  double GetTiltAngle( int lnum ) const;
  double GetRotAngle1( int lnum ) const;
  double GetRotAngle2( int lnum ) const;
  const ThreeVector & GetGlobalPosition( int lnum ) const;
  ThreeVector NormalVector( int lnum ) const;
  ThreeVector UnitVector( int lnum ) const;
  const DCGeomRecord *GetRecord( int lnum ) const;

  ThreeVector Local2GlobalPos( int lnum, const ThreeVector &in ) const;
  ThreeVector Global2LocalPos( int lnum, const ThreeVector &in ) const;
  ThreeVector Local2GlobalDir( int lnum, const ThreeVector &in ) const;
  ThreeVector Global2LocalDir( int lnum, const ThreeVector &in ) const;


  std::vector<int> GetDetectorIDList( void ) const;
  int GetTofId( void ) const { return TOFid_; }
  
private:
  static DCConfMan *confMan_;
  std::string filename_;
  mutable std::map <int, DCGeomRecord *> geomRecord_;
  int TOFid_;

private:
  void clearElements( void );
};


#endif
