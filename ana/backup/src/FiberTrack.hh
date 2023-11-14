/*
  FiberTrack.hh
*/

#ifndef FiberTrack_h
#define FiberTrack_h 1

#include <vector>

class FiberHit;

class FiberTrack
{
public:
  explicit FiberTrack();
  ~FiberTrack();
private:
  FiberTrack( const FiberTrack & );
  FiberTrack & operator = ( const FiberTrack & );

private:
  std::vector <FiberHit *> hitArray;

public:
  void AddHit( FiberHit *hitp ) { hitArray.push_back( hitp ); }
  bool DoFit( void );
  std::size_t GetNHit( void ) const { return hitArray.size(); }
  FiberHit * GetHit( std::size_t nth ) const;
  FiberHit * GetHitOfLayerNumber( int lnum ) const;

  double GetX0( void ) const { return x0_; }
  double GetY0( void ) const { return y0_; }
  double GetU0( void ) const { return u0_; }
  double GetV0( void ) const { return v0_; }
  double GetChiSquare( void ) const { return chisqr_; }
  double GetX( double z ) const { return x0_+u0_*z; } 
  double GetY( double z ) const { return y0_+v0_*z; } 
  bool GetStatus( void ) const { return status_; } 
private:
  bool status_;
  double x0_, y0_, u0_, v0_;
  double chisqr_;

};


#endif
