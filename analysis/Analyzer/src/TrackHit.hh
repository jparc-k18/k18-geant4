/*
  TrackHit.hh
*/

#ifndef TrackHit_h
#define TrackHit_h 1

#include "DCHit.hh" 
#include "ThreeVector.hh"

class TrackHit
{
public:
  explicit TrackHit( DCHit *hit );
  ~TrackHit();
private:
  TrackHit( const TrackHit & );
  TrackHit & operator= ( const TrackHit & );

public:
  DCHit *GetDCHit( void ) { return dchitp_; }
  void SetCalGPos( const ThreeVector &pos ) { calGPos_=pos; }
  void SetCalLPos( double pos ) { calLPos_=pos; }

  int GetLayer( void ) const { return dchitp_->GetLayer(); }
  double GetLocalHitPos( void ) const 
  { return dchitp_->GetLocalHitPos(); }
  const ThreeVector & GetCalGPos( void ) const 
  { return calGPos_; }
  double GetCalLPos( void ) const { return calLPos_; }
  double GetResidual( void ) const 
  { return dchitp_->GetLocalHitPos()-calLPos_; }

private:
  DCHit *dchitp_;
  ThreeVector calGPos_;
  double calLPos_;

};

#endif

