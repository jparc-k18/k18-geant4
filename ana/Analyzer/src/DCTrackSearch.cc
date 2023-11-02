/*
  DCTrackSearch.cc
*/

#include "DCTrackSearch.hh"
#include "DCParameters.hh"
#include "DCHit.hh"
#include "DCLocalTrack.hh"
#include "TemplateLib.hh"
#include "DetectorId.hh"
#include "DCConfMan.hh"


#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);
const double MaxChisquare = 10000.;
const double MaxNumberOfClusters = 10;

//For BC3&4, SDC1&2
int LocalTrackSearch( const DCHitContainer * HC,
		      const DCPairPlaneInfo * PpInfo,
		      int npp, std::vector <DCLocalTrack *> &TrackCont,
		      int MinNumOfHits )
{
  static const std::string funcname = "[LocalTrackSearch]";

  std::vector < std::vector <DCHit *> > CandCont;
  CandCont.resize(npp);

  for( int i=0; i<npp; ++i ){
    bool ppFlag=PpInfo[i].flag;
    int layer1=PpInfo[i].id1, layer2=PpInfo[i].id2;
    if(ppFlag) {
      std::cerr << funcname << " Invalid ppFlag, layer1 = "
		<< layer1 << ", layer2 = " << layer2 << std::endl;
    }
    else{
      MakeUnPairPlaneHitCluster( HC[layer1], CandCont[i] );
      //if(layer1>12) std::cout<<"layer1="<<layer1<<std::endl;
    }
  }

  std::vector <int> nCombi(npp);
  for( int i=0; i<npp; ++i ){ 
    nCombi[i]=(CandCont[i]).size();

    // If #Cluster>MaxNumberOfCluster,  error return

    if(nCombi[i]>MaxNumberOfClusters){
      for( int i=0; i<npp; ++i )
	for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
      return 0;
    } 
  }

#if 0
  std::cout << funcname << ": #Hits of each group" << std::endl;
  for( int i=0; i<npp; ++i ) std::cout << std::setw(4) << nCombi[i];
  std::cout << std::endl;
  for( int i=0; i<npp; ++i ){
    int n=CandCont[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((DCLTrackHit *)CandCont[i][j]->GetHit(0))->GetWire() << " ";
      //      std::cout << CandCont[i][j] << " ";
    }
    std::cout << std::endl;
  }
#endif

  std::vector < std::vector <int> > 
    CombiIndex = makeindex( npp, &nCombi[0] );
  int nnCombi=CombiIndex.size();

#if 0
  std::cout << " ===> " << nnCombi << " combinations will be checked.." 
	    << std::endl;
#endif

  for( int i=0; i<nnCombi; ++i ){
    DCLocalTrack *track = MakeTrack( CandCont, &((CombiIndex[i])[0]));
    if( !track ) continue;
    if( track->GetNHit()>=MinNumOfHits && track->DoFit() &&
	track->GetChiSquare()<MaxChisquare ){
      TrackCont.push_back(track);
    }    
    else
      delete track;
  }

  // Clear Flags
  int nbefore=TrackCont.size();
  for( int i=0; i<nbefore; ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
  }

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": Before Sorting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      DCLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  partial_sort( TrackCont.begin(), TrackCont.end(), 
		TrackCont.end(), DCLTrackComp() );

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Sorting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      DCLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  // Delete Duplicated Tracks

  for( int i=0; i<int(TrackCont.size()); ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();

    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
      DCLocalTrack *tp2=TrackCont[i2];
      int nh2=tp2->GetNHit(), flag=0;
      for( int j=0; j<nh2; ++j )
	if( tp2->GetHit(j)->showFlags() ) ++flag;
      if(flag){
	delete tp2;
	TrackCont.erase(TrackCont.begin()+i2);
      }
    }      
  }

  {
    int nn=TrackCont.size();
    for(int i=0; i<nn; ++i ){
      DCLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	int lnum = tp->GetHit(j)->GetLayer();
	double zz = DCConfMan::GetInstance().GetLocalZ( lnum );
	//	std::cout<<"zz="<<zz<<std::endl;
	tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
      }
    }
  }

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Deleting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      DCLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

#if 0
  for( int i=0; i<npp; ++i )
    for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
#endif
  
  return TrackCont.size();
}

#if 0
//For SDC3&4
int SdcOutLocalTrackSearch( const DCHitContainer * HC,
			    std::vector <DCLocalTrack *> &TrackCont)
		      
{
  static const std::string funcname = "[LocalTrackSearch]";

  std::vector <int> UsedPlane;
  std::vector <int> UnUsedPlane;

  double meanMultiSdc3=0.0;
  for (int i=1; i<=6; i++) {
    meanMultiSdc3 += (double)HC[i].size();
  }
  meanMultiSdc3 /= 6.;
  for (int i=1; i<=6; i++) {
    if ((HC[i].size() < meanMultiSdc3+2.) && UsedPlane.size()<3)
      UsedPlane.push_back(i);
    else
      UnUsedPlane.push_back(i);
  }

  double meanMultiSdc4=0.0;
  for (int i=7; i<=12; i++) {
    meanMultiSdc4 += (double)HC[i].size();
  }
  meanMultiSdc4 /= 6.;
  for (int i=12; i>=7; i--) {
    if ((HC[i].size() < meanMultiSdc4+2.) 
	&& (HC[i].size() > 0) 
	&& UsedPlane.size()<6)
      UsedPlane.push_back(i);
    else
      UnUsedPlane.push_back(i);
  }
  
  std::vector < std::vector <DCHit *> > CandCont;
  CandCont.resize(UsedPlane.size());

  for( int i=0; i<UsedPlane.size(); ++i ){
    MakeUnPairPlaneHitCluster( HC[UsedPlane[i]], CandCont[i] );
  }

  std::vector <int> nCombi(UsedPlane.size());
  for( int i=0; i<UsedPlane.size(); ++i ){ 
    nCombi[i]=(CandCont[i]).size();

    // If #Cluster>MaxNumerOfCluster,  error return

    if(nCombi[i]>MaxNumberOfClusters){
      for( int i=0; i<UsedPlane.size(); ++i )
	for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
      return 0;
    } 
  }

#if 0
  std::cout << funcname << ": Used Plane" << std::endl;
  for( int i=0; i<UsedPlane.size(); ++i ) std::cout << std::setw(4) << UsedPlane[i];
  std::cout << std::endl;

  std::cout << funcname << ": #Hits of each group" << std::endl;
  for( int i=0; i<UsedPlane.size(); ++i ) std::cout << std::setw(4) << nCombi[i];
  std::cout << std::endl;
  for( int i=0; i<UsedPlane.size(); ++i ){
    int n=CandCont[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((DCLTrackHit *)CandCont[i][j]->GetHit(0))->GetWire() << " ";
    }
    std::cout << std::endl;
  }
#endif

  std::vector < std::vector <int> > 
    CombiIndex = makeindex_SdcOut( UsedPlane.size(), 
				   MinNumOfHitsSdcOut, 
				   UsedPlane.size(), 
				   &(nCombi[0]) );
  int nnCombi=CombiIndex.size();
 
#if 0
  std::cout << " ===> " << nnCombi << " combinations will be checked.." 
	    << std::endl;
#endif

  for( int i=0; i<nnCombi; ++i ){
    DCLocalTrack *track = MakeTrack( CandCont, &((CombiIndex[i])[0]) );
    if( !track ) continue;
    if( track->GetNHit()>=MinNumOfHitsSdcOut && track->DoFit() &&
	track->GetChiSquare()<MaxChisquare ){
      TrackCont.push_back(track);
      double chisqr = track->GetChiSquare();
    }
    else{
      //      std::cout << "No tracks available" << std::endl;
      delete track;
    }
  }

  // Clear Flags
  int nbefore=TrackCont.size();
  for( int i=0; i<nbefore; ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
  }

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": Before Sorting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      DCLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  partial_sort( TrackCont.begin(), TrackCont.end(), 
		TrackCont.end(), DCLTrackCompSdcOut() );

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Sorting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      DCLocalTrack *track=TrackCont[i];
      if (track->GetChiSquare()<20) {
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare();
	int nh=track->GetNHit();
	for( int ih=0; ih<nh; ++ih ){
	  DCLTrackHit *hit=track->GetHit(ih);
	  int layerId=hit->GetLayer();  
	  std::cout << " " << layerId;
	}
	std::cout << std::endl;
      }
    }
    std::cout << std::endl;

  }
#endif

  // Delete Duplicated Tracks

  for( int i=0; i<int(TrackCont.size()); ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();

    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
      DCLocalTrack *tp2=TrackCont[i2];
      int nh2=tp2->GetNHit(), flag=0;
      for( int j=0; j<nh2; ++j )
	if( tp2->GetHit(j)->showFlags() ) ++flag;
      if(flag){
	delete tp2;
	TrackCont.erase(TrackCont.begin()+i2);
      }
    }      
  }

  {
    int nn=TrackCont.size();
    for(int i=0; i<nn; ++i ){
      DCLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	int lnum = tp->GetHit(j)->GetLayer();
	double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
	tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
      }
    }
  }

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Deleting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      DCLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  int nt=TrackCont.size();
  for (int it=0; it<nt; it++) {
    const DCGeomMan & geomMan = DCGeomMan::GetInstance();
    DCLocalTrack *tp=TrackCont[it];
    double x0=tp->GetX0(), y0=tp->GetY0();
    double u0=tp->GetU0(), v0=tp->GetV0();
 
    for (int ip = 0; ip < UnUsedPlane.size(); ip++) {
      int lnum = UnUsedPlane[ip]+PlOffsSdcOut;
      double angle=geomMan.GetTiltAngle( lnum );
      double z = geomMan.GetLocalZ( lnum );      
      double x = x0+u0*z;
      double y = y0+v0*z;
      double localCalPos = x*cos(angle*Deg2Rad)+y*sin(angle*Deg2Rad);
      /*
      std::cout << "Layer# = " << lnum 
		<< ", LocalCalPos = " << localCalPos
		<< std::endl;
      */
      double mindiff=1000.;
      int ihUse = -1;
      double wpUse, dlUse;
      int mUse;
      double maxDiff=1.5;
      for (int ih=0; ih<HC[UnUsedPlane[ip]].size(); ih++ ) {
	DCHit *hit = HC[UnUsedPlane[ip]][ih];
	if ( hit ) {
	  int multi = hit->GetDriftLengthSize();
	  for (int m=0; m<multi; m++) {
	    if( !(hit->rangecheck(m)) ) continue;
	    if( (hit->showFlags(m)) ) continue;
	    double wp=hit->GetWirePosition();
	    double dl=hit->GetDriftLength(m);
	    double diff1 = fabs((wp+dl) - localCalPos);
	    double diff2 = fabs((wp-dl) - localCalPos);
	    /*
	    std::cout << ih << "th Hit " << std::endl;
	    std::cout << "localPos1 = " << wp+dl 
		      << ", Diff = " << (wp+dl) - localCalPos << std::endl;
	    std::cout << "localPos2 = " << wp+dl 
		      << ", Diff = " << (wp-dl) - localCalPos << std::endl;
	    */
	    if (diff1 < maxDiff && diff1 < mindiff ) {
	      mindiff = diff1;
	      ihUse = ih; mUse = m;
	      wpUse = wp; dlUse = dl;
	    }
	    if (diff2 < maxDiff && diff2 < mindiff ) {
	      mindiff = diff2;
	      ihUse = ih; mUse = m;
	      wpUse = wp; dlUse = -dl;
	    }

	  }
	}
      }
      if (ihUse>=0 && ihUse < HC[UnUsedPlane[ip]].size()) {
	DCHit *hit = HC[UnUsedPlane[ip]][ihUse];
	DCLTrackHit *lhit = new DCLTrackHit(hit, wpUse+dlUse, mUse);
	lhit->setFlags();
	if(lhit) tp->AddHit( lhit );
      }
    }
    tp->DoFit();
  }

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Adding UnUsed Plane. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      DCLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
  getchar();
#endif

  for( int i=0; i<UsedPlane.size(); ++i )
    for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
  
  return TrackCont.size();
}
#endif

#if 0 
bool MakePairPlaneHitCluster( const DCHitContainer & HC1,
			      const DCHitContainer & HC2,
			      double CellSize,
			      std::vector <DCPairHitCluster *> & Cont )
{
  int nh1=HC1.size(), nh2=HC2.size();
  std::vector <int> UsedFlag(nh2,0);   
  
  for( int i1=0; i1<nh1; ++i1 ){
    DCHit *hit1=HC1[i1];

    double wp1=hit1->GetWirePosition();
    bool flag=false;
    for( int i2=0; i2<nh2; ++i2 ){
      DCHit *hit2=HC2[i2];
      double wp2=hit2->GetWirePosition();
      if( fabs(wp1-wp2)<CellSize ){

	int multi1 = hit1->GetDriftLengthSize();
	int multi2 = hit2->GetDriftLengthSize();
	for (int m1=0; m1<multi1; m1++) {
	  if( !(hit1->rangecheck(m1)) ) continue;
	  for (int m2=0; m2<multi2; m2++) {
	    if( !(hit2->rangecheck(m2)) ) continue;
	    double x1,x2;
	    if( wp1<wp2 ){
	      x1=wp1+hit1->GetDriftLength(m1);
	      x2=wp2-hit2->GetDriftLength(m2);
	    }
	    else {
	      x1=wp1-hit1->GetDriftLength(m1);
	      x2=wp2+hit2->GetDriftLength(m2);
	    }
	    
	    Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,x1,m1),
						  new DCLTrackHit(hit2,x2,m2) ) );
	    flag=true; ++UsedFlag[i2];
	  }
	}
      }
    }      
#if 0
    if(!flag){
      int multi1 = hit1->GetDriftLengthSize();
      for (int m1=0; m1<multi1; m1++) {
	if( !(hit1->rangecheck(m1)) ) continue;
	double dl=hit1->GetDriftLength(m1);
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,wp1+dl,m1) ) );
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,wp1-dl,m1) ) );
      }
    }
#endif
  }
#if 0
  for( int i2=0; i2<nh2; ++i2 ){
    if( UsedFlag[i2]==0 ) {
      DCHit *hit2=HC2[i2];
      int multi2 = hit2->GetDriftLengthSize();
      for (int m2=0; m2<multi2; m2++) {
	if( !(hit2->rangecheck(m2)) ) continue;
	
	double wp=hit2->GetWirePosition();
	double dl=hit2->GetDriftLength(m2);
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit2,wp+dl,m2) ) );
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit2,wp-dl,m2) ) );
      }
    }
  }
#endif

  return true;
}
#endif 

#if 0
bool MakeUnPairPlaneHitCluster( const DCHitContainer & HC,
				std::vector <DCPairHitCluster *> & Cont )
{
  int nh=HC.size();

  for( int i=0; i<nh; ++i ){
    DCHit *hit=HC[i];
    if( hit ){
      int multi = hit->GetDriftLengthSize();
      for (int m=0; m<multi; m++) {
	if( !(hit->rangecheck(m)) ) continue;

	double wp=hit->GetWirePosition();
	double dl=hit->GetDriftLength(m);
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit,wp+dl,m) ) );
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit,wp-dl,m) ) );
      }
    }
  }

  return true;
}
#endif

bool MakeUnPairPlaneHitCluster( const DCHitContainer & HC,
				std::vector <DCHit *> & Cont )
{
  int nh=HC.size();

  for( int i=0; i<nh; ++i ){
    DCHit *hit=HC[i];
    if( hit ){
	Cont.push_back( hit );
    }
  }

  return true;
}


std::vector< std::vector<int> > makeindex( int ndim, const int *index1 )
{
  if(ndim==1){
    std::vector< std::vector<int> > index2;
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem(1,i);
      index2.push_back(elem);
    }
    return index2;
  }

  std::vector< std::vector<int> > 
    index2=makeindex( ndim-1, index1+1 );

  std::vector< std::vector<int> > index;
  int n2=index2.size();
  for( int j=0; j<n2; ++j ){
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem;
      int n3=index2[j].size();
      elem.reserve(n3+1);
      elem.push_back(i);
      for( int k=0; k<n3; ++k )
	elem.push_back(index2[j][k]);
      index.push_back(elem);
      int size1=index.size();
    }
  }

  return index;
}  

std::vector< std::vector<int> > 
makeindex_SdcOut( int ndim_org, int minimumHit, int ndim, const int *index1 )
{
  if(ndim==1){
    std::vector< std::vector<int> > index2;
    for( int i=0; i<index1[0]; ++i ){
      std::vector <int> elem(1,i);
      index2.push_back(elem);
    }
    return index2;
  }

  std::vector< std::vector<int> > 
    index2=makeindex_SdcOut( ndim_org, minimumHit, ndim-1, index1+1 );
 
  std::vector< std::vector<int> > index;
  int n2=index2.size();
  for( int j=0; j<n2; ++j ){
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem;
      int validHitNum=0;
      int n3=index2[j].size();
      elem.reserve(n3+1);
      elem.push_back(i);
      if (i != -1)
	validHitNum++;
      for( int k=0; k<n3; ++k ) {
        elem.push_back(index2[j][k]);
        if (index2[j][k] != -1)
          validHitNum++;
      }
      if (ndim==ndim_org) {
        if (validHitNum >= minimumHit)
          index.push_back(elem);
      } else {
        index.push_back(elem);
      }
      int size1=index.size();
    }
  }

  return index;
}  

DCLocalTrack 
* MakeTrack(  const std::vector < std::vector <DCHit *> > &CandCont,
	      const int *combination )
{
  static const std::string funcname = "[MakeTrack]";

  DCLocalTrack *tp=new DCLocalTrack();
  if(!tp){
    std::cerr << funcname << ": new fail" << std::endl;
    return 0;
  }    
  int n=CandCont.size();

  for( int i=0; i<n; ++i ){
    int m=combination[i];
    DCHit *hit=0;
    if(m>=0) hit=CandCont[i][m];
#if 0
    std::cout << funcname << ":" << std::setw(3)
	      << i << std::setw(3) << m  << " "
	      << CandCont[i][m] << std::endl; 
#endif

    if(hit) tp->AddHit( hit );
  }

  return tp;
}
