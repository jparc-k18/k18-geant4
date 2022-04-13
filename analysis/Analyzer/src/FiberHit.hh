/*
  FiberHit.hh
*/

#ifndef FiberHit_h 
#define FiberHit_h

class FiberHit
{
public:
  explicit FiberHit();
  FiberHit( int layer, int segment, double edep, double time, int pid );
  ~FiberHit();

private:
  FiberHit( const FiberHit & );
  FiberHit & operator = ( const FiberHit & );

private:
  int layer_, segment_;
  double x_, y_, z_, r_;
  double edep_, time_;
  int    pid_;
public:
  void SetLayer( int layer ) { layer_=layer; }
  void SetSegment( int segment ) { segment_=segment; }
  void SetEdep( double de ) { edep_=de; }
  void SetTime( double time ) { time_=time; }
  void SetParticleId( int id ) { pid_=id; }

  int GetLayer( void ) const { return layer_; }
  int GetSegment( void ) const { return segment_; }
  double GetEdep( void ) const { return edep_; }
  double GetTime( void ) const { return time_; }
  int GetParticleId( void ) const { return pid_; }
};

#endif
