/*
  K18Track.hh
*/

#ifndef K18Track_h
#define K18Track_h 1

#include "ThreeVector.hh"

#include <vector>
#include <functional>

class DCHit;
class DCLocalTrack;
class TrackHit;
class DCAnalyzer;

class K18Track
{
public: 
  K18Track( DCLocalTrack *tin, DCLocalTrack *tout, double P0, DCHit *SF1, DCHit* SF2);
  ~K18Track(); 

private:
  K18Track( const K18Track & );
  K18Track & operator = ( const K18Track & );

private:
  DCLocalTrack *TrIn_, *TrOut_;
  DCHit *sf1_, *sf2_;
  double P0_;
  bool Status_;
  double chisqr_;
  double Xi_, Yi_, Ui_, Vi_;
  double Xo_, Yo_, Uo_, Vo_;
  double Delta_; 

  bool StatusD2U_;
  double P0D2U_;
  double DeltaD2U_; 

  bool gfastatus_;

  std::vector <TrackHit *> hitContIn, hitContOut;
  
public:
  bool doFit( void );
  bool doFitD2U( void );
  bool CalcMomentumD2U( void );
  bool Status( void ) const { return Status_; }
  bool StatusD2U( void ) const { return StatusD2U_; }

  DCLocalTrack *TrackIn( void )  { return TrIn_; }
  DCLocalTrack *TrackOut( void ) { return TrOut_; }

  double Xin( void ) const { return Xi_; }
  double Yin( void ) const { return Yi_; }
  double Uin( void ) const { return Ui_; }
  double Vin( void ) const { return Vi_; }
  double Xout( void ) const { return Xo_; }
  double Yout( void ) const { return Yo_; }
  double Uout( void ) const { return Uo_; }
  double Vout( void ) const { return Vo_; }
  double Delta( void ) const { return Delta_; }
  double P( void ) const { return P0_*(1.+Delta_); }
  double chisquare() const { return chisqr_; }

  double DeltaD2U( void ) const { return DeltaD2U_; }
  double PD2U( void ) const { return P0_*(1.+DeltaD2U_); }

  ThreeVector BeamMomentum( void ) const;
  double Xtgt( void ) const;
  double Ytgt( void ) const;
  double Utgt( void ) const;
  double Vtgt( void ) const;

  int GetNHitsIn( void )  const { return hitContIn.size(); }
  int GetNHitsOut( void ) const { return hitContOut.size(); }
  int GetNHitsTotal( void ) const
  { return hitContIn.size()+hitContOut.size(); }

  TrackHit * GetK18HitIn( int i );
  TrackHit * GetK18HitOut( int i );
  TrackHit * GetK18HitTotal( int i );
  TrackHit * GetK18HitByPlaneId( int PlId );

  bool GoodForAnalysis( void ) const { return gfastatus_; }
  bool GoodForAnalysis( bool status )
  { bool ret=gfastatus_; gfastatus_=status; return ret; } 

  //bool ReCalc( bool ApplyRecursively=false );

private:
  void deleteHits( void );
  void addHits( void );

  // for DS
// private:
//   unsigned int DSKey_;
//   struct DSFormat {
//     unsigned int header_;
//     unsigned int pkeyIn_, pkeyOut_;
//     float xi_, yi_, ui_, vi_;
//     float xo_, yo_, uo_, vo_;
//     float p0_, delta_, chisqr_;
//     short status_, gfastatus_;
//   };
// public:
//   explicit K18Track( unsigned int *bufp, DCAnalyzer *DCana )
//   { DSRestore( bufp, DCana ); }
//   std::size_t DSSize( void ) const;
//   std::size_t DSSave( unsigned int *bufp ) const;
//   bool DSRestore( unsigned int *bufp, DCAnalyzer *DCana );
//   void DSSetSaveKey( unsigned int key ) { DSKey_ = key; }
//   unsigned int DSGetSaveKey( void ) const { return DSKey_; } 

};

struct K18TrackComp
  : public std::binary_function <K18Track *, K18Track *,bool>
{
  bool operator()( const K18Track * const p1,
		   const K18Track * const p2 ) const
  {
    return (p1->chisquare())<(p2->chisquare());
  }

};

#endif
