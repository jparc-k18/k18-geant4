/*
  DCLocalTrack.hh
*/

#ifndef DCLocalTrack_h
#define DCLocalTrack_h 1

#include <vector>
#include <functional>


class DCHit;

class DCLocalTrack
{
public:
  explicit DCLocalTrack();
  ~DCLocalTrack();
private:
  DCLocalTrack( const DCLocalTrack & );
  DCLocalTrack & operator = ( const DCLocalTrack & );

private:
  std::vector <DCHit *> hitArray;

public:
  void AddHit( DCHit *hitp ) { hitArray.push_back( hitp ); }
  bool DoFit( void );
  std::size_t GetNHit( void ) const { return hitArray.size(); }
  DCHit * GetHit( std::size_t nth ) const;
  DCHit * GetHitOfLayerNumber( int lnum ) const;

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

struct DCLTrackComp 
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
{
  bool operator()( const DCLocalTrack * const p1, 
		   const DCLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    if(n1>n2+1) return true;
    else if(n2>n1+1) return false;
    else
      return (p1->GetChiSquare())<(p2->GetChiSquare());
  }

};

struct DCLTrackCompSdcOut 
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
{
  bool operator()( const DCLocalTrack * const p1, 
		   const DCLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    double chi1=p1->GetChiSquare(), chi2=p2->GetChiSquare();
    double chi_th=3.;

    if (chi1>chi_th || chi2>chi_th) {
      if(n1>n2+1) return true;
      else if(n2>n1+1) return false;
      else
	return (p1->GetChiSquare())<(p2->GetChiSquare());
    }
    else {
      if(n1>n2) return true;
      else if(n2>n1) return false;
      else
	return (p1->GetChiSquare())<(p2->GetChiSquare());
    }
  }
};

#endif
