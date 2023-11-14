/*
  DCHit.hh
*/

#ifndef DCHit_h 
#define DCHit_h

class DCHit
{
public:
  explicit DCHit();
  DCHit( int layer, int wire, int tdc ); 
  ~DCHit();

private:
  DCHit( const DCHit & );
  DCHit & operator = ( const DCHit & );

private:
  int layer_, wire_, tdc_;
  double dt_, dl_;
  double angle_, xl_;
  double xcal_, ycal_; 
  bool belongTrack_;

public:
  void SetLayer( int layer ) { layer_=layer; }
  void SetWire( int wire ) { wire_=wire; }
  void SetTdcVal( int tdc ) { tdc_=tdc; }
  void SetDriftTime( double dt ) { dt_=dt; }
  void SetDriftLength( double dl ) { dl_=dl; }
  void SetTiltAngle( double angleDegree ) { angle_=angleDegree; }
  void SetLocalHitPos( double xl ) { xl_=xl; }
  void SetCalPosition( double x, double y ) { xcal_=x; ycal_=y; }

  int GetLayer( void ) const { return layer_; }
  int GetWire( void ) const { return wire_; }
  int GetTdcVal( void ) const { return tdc_; }
  double GetDriftTime( void ) const { return dt_; }
  double GetDriftLength( void ) const { return dl_; }
  double GetTiltAngle( void ) const { return angle_; }
  double GetLocalHitPos( void ) const { return xl_; }
  double GetXcal( void ) const { return xcal_; }
  double GetYcal( void ) const { return ycal_; }
  double GetResidual( void ) const;

  void setFlags( void ) { belongTrack_=true; }
  void clearFlags( void ) { belongTrack_=false; }
  bool showFlags( void ) const { return belongTrack_; }


};

#endif
