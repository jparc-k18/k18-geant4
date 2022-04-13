/*
  DCAnalyzer.hh
*/

#ifndef DCAnalyzer_h 
#define DCAnalyzer_h 1

#include "DetectorId.hh"
#include "ThreeVector.hh"
#include <vector>

class DCHit;
class DCLocalTrack;
class SksTrack;
class K18Track;
class RawData;
//class SimuData;

typedef std::vector <DCHit *> DCHitContainer;


class DCAnalyzer
{
public:
  DCAnalyzer();
  ~DCAnalyzer();
private:
  DCAnalyzer( const DCAnalyzer & );
  DCAnalyzer & operator = ( const DCAnalyzer & );

private:
  // DCHit
  DCHitContainer BcInHC[NumOfLayersBcIn+1], BcOutHC[NumOfLayersBcOut+1];
  DCHitContainer SdcInHC[NumOfLayersSdcIn+1], SdcOutHC[NumOfLayersSdcOut+1];

  // DCLocalTrack
  std::vector <DCLocalTrack *> TrackBcInCol, TrackBcOutCol;
  std::vector <DCLocalTrack *> TrackSdcInCol, TrackSdcOutCol;

  // SksTrack
  std::vector <SksTrack *> SksTrackCol;

  // K18Track
  std::vector <K18Track *> K18TrackCol;

public:
  ///bool DecodeRawHits( RawData *rawData );
  /* New Added by K.Miwa */
 //  bool DecodeSimuHits( SimuData *simuData );
  bool SetDCHits( int lnum, double lpos );

  inline const DCHitContainer & GetBcInHC( int layer ) const;
  inline const DCHitContainer & GetBcOutHC( int layer ) const;
  inline const DCHitContainer & GetSdcInHC( int layer ) const;
  inline const DCHitContainer & GetSdcOutHC( int layer ) const;

  bool TrackSearchBcIn( void );
  bool TrackSearchBcOut( void );
  bool TrackSearchSdcIn( void );
  bool TrackSearchSdcOut( void );

  int GetNtracksBcIn( void ) const  { return TrackBcInCol.size(); }
  int GetNtracksBcOut( void ) const { return TrackBcOutCol.size(); }
  int GetNtracksSdcIn( void ) const  { return TrackSdcInCol.size(); }
  int GetNtracksSdcOut( void ) const { return TrackSdcOutCol.size(); }

  inline DCLocalTrack * GetTrackBcIn( int i ) const;
  inline DCLocalTrack * GetTrackBcOut( int i ) const;
  inline DCLocalTrack * GetTrackSdcIn( int i ) const;
  inline DCLocalTrack * GetTrackSdcOut( int i ) const;

  bool TrackSearchSks( void );
  bool TrackSearchK18( void );

  int GetNTracksSks( void ) const { return SksTrackCol.size(); }
  int GetNTracksK18( void ) const { return K18TrackCol.size(); }

  inline SksTrack * GetSksTrack( int i ) const;
  inline K18Track * GetK18Track( int i ) const;

private:
  void clearDCHits( void );
  void clearTracksBcIn( void );
  void clearTracksBcOut( void );
  void clearTracksSdcIn( void );
  void clearTracksSdcOut( void );
  void clearSksTracks( void );
  void clearK18Tracks( void );

public:
  void resetTracksBcIn( void ) { clearTracksBcIn(); }
  void resetTracksBcOut( void ) { clearTracksBcOut(); }
  void resetTracksSdcIn( void ) { clearTracksSdcIn(); }
  void resetTracksSdcOut( void ) { clearTracksSdcOut(); }


};


inline const DCHitContainer & DCAnalyzer::GetBcInHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersBcIn ) layer=0;
  return BcInHC[layer];
}

inline const DCHitContainer & DCAnalyzer::GetBcOutHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersBcOut ) layer=0;
  return BcOutHC[layer];
}

inline const DCHitContainer & DCAnalyzer::GetSdcInHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersSdcIn ) layer=0;
  return SdcInHC[layer];
}

inline const DCHitContainer & DCAnalyzer::GetSdcOutHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersSdcOut ) layer=0;
  return SdcOutHC[layer];
}

inline DCLocalTrack * DCAnalyzer::GetTrackBcIn( int i ) const
{
  if( i>=0 && i<TrackBcInCol.size() )
    return TrackBcInCol[i];
  else
    return 0;
}

inline DCLocalTrack * DCAnalyzer::GetTrackBcOut( int i ) const
{
  if( i>=0 && i<TrackBcOutCol.size() )
    return TrackBcOutCol[i];
  else
    return 0;
}

inline DCLocalTrack * DCAnalyzer::GetTrackSdcIn( int i ) const
{
  if( i>=0 && i<TrackSdcInCol.size() )
    return TrackSdcInCol[i];
  else
    return 0;
}

inline DCLocalTrack * DCAnalyzer::GetTrackSdcOut( int i ) const
{
  if( i>=0 && i<TrackSdcOutCol.size() )
    return TrackSdcOutCol[i];
  else
    return 0;
}

inline SksTrack * DCAnalyzer::GetSksTrack( int i ) const
{
  if( i>=0 && i<SksTrackCol.size() )
    return SksTrackCol[i];
  else
    return 0;
}

inline K18Track * DCAnalyzer::GetK18Track( int i ) const
{
  if( i>=0 && i<K18TrackCol.size() )
    return K18TrackCol[i];
  else
    return 0;
}


#endif 
