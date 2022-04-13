/*
  DCGeomMan.hh
*/

#ifndef DCGeomMan_h
#define DCGeomMan_h 1

#include <string>
#include <vector>
#include <map>

#include "G4ThreeVector.hh"

class DCGeomRecord;

class DCGeomMan 
{
private:
  DCGeomMan();
public:
  ~DCGeomMan();

public:
  void SetFileName( const char *filename ) { filename_=filename; }
  void SetFileName( const std::string &filename ) { filename_=filename; }

  bool Initialize( void );
  bool Initialize( const char *filename )
  { filename_=filename; bool ret = Initialize(); return ret; }
  bool Initialize( const std::string &filename )
  { filename_=filename; bool ret = Initialize(); return ret; }

  static DCGeomMan & GetInstance( void );
  double GetLocalZ( int lnum ) const;
  double GetResolution( int lnum ) const;
  double GetTiltAngle( int lnum ) const;
  double GetRotAngle1( int lnum ) const;
  double GetRotAngle2( int lnum ) const;
  const G4ThreeVector & GetGlobalPosition( int lnum ) const;
  G4ThreeVector NormalVector( int lnum ) const;
  G4ThreeVector UnitVector( int lnum ) const;
  const DCGeomRecord *GetRecord( int lnum ) const;

  G4ThreeVector Local2GlobalPos( int lnum, const G4ThreeVector &in ) const;
  G4ThreeVector Global2LocalPos( int lnum, const G4ThreeVector &in ) const;
  G4ThreeVector Local2GlobalDir( int lnum, const G4ThreeVector &in ) const;
  G4ThreeVector Global2LocalDir( int lnum, const G4ThreeVector &in ) const;

  double calcWirePosition( int lnum, int wire ) const;
  int calcWireNumber( int lnum, double position ) const;

  std::vector<int> GetDetectorIDList( void ) const;
  int GetDetectorId( const std::string &detName ) const;
  int GetTofId( void ) const { return TOFid_; }
  int GetSMFId( void ) const { return SMFid_; }
  
private:
  static DCGeomMan *geomMan_;
  std::string filename_;
  mutable std::map <int, DCGeomRecord *> geomRecord_;
  int TOFid_, SMFid_;

private:
  void clearElements( void );
};


#endif
