/*
  DCTrackSearch.hh
*/

#ifndef DCTrackSearch_h
#define DCTrackSearch_h 1

#include "DCAnalyzer.hh"

#include <vector>
#include <functional>
#include "ThreeVector.hh"

struct DCPairPlaneInfo;
class DCLocalTrack;

//For BC3&4 SDC1&2
int LocalTrackSearch( const DCHitContainer * HC,  
		      const DCPairPlaneInfo * PpInfo,
		      int npp, std::vector <DCLocalTrack *> & TrackCont,
		      int MinNumOfHits=6 );
#if 0
//For SDC3&4
int SdcOutLocalTrackSearch( const DCHitContainer * HC,
			    std::vector <DCLocalTrack *> &TrackCont );
#endif

#if 0
bool MakePairPlaneHitCluster( const DCHitContainer & HC1,
			      const DCHitContainer & HC2,
			      double CellSize,
			      std::vector <DCPairHitCluster *> & Cont );
#endif

bool MakeUnPairPlaneHitCluster( const DCHitContainer & HC,
				std::vector <DCHit *> & Cont );

DCLocalTrack *MakeTrack( const std::vector < std::vector <DCHit *> > &CandCont,
			 const int *combination );

std::vector< std::vector<int> > makeindex( int ndim, const int *index1 ); 
std::vector< std::vector<int> > makeindex_SdcOut( int ndim_org, int minimumHit, int ndim, const int *index1 ); 


#endif
