/*
  DCAnalyzer.cc
*/

#include "DCAnalyzer.hh"
#include "DCHit.hh"
#include "DCLocalTrack.hh"
#include "SksTrack.hh"
#include "K18Track.hh"
#include "K18Parameters.hh"
#include "TemplateLib.hh"
#include "DCTrackSearch.hh"
#include "DCConfMan.hh"
#include "myfunc.hh"

//#include "Randomize.h"
//#include "SimuData.hh"

#define DefStatic
#include "DCParameters.hh"
#undef DefStatic

//#include "SksObjectId.hh"

#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <sstream>

const double MaxChiSqrSksTrack = 10000.;

DCAnalyzer::DCAnalyzer()
{}

DCAnalyzer::~DCAnalyzer()
{
  clearSksTracks();
  clearK18Tracks();
  clearTracksSdcOut();
  clearTracksSdcIn();
  clearTracksBcOut();
  clearTracksBcIn();

  clearDCHits();
}

#if 0
bool DCAnalyzer::DecodeRawHits( RawData *rawData )
{
  const std::string funcname = "[DCAnalyzer::DecodeRawHits]";

  clearDCHits();

  // BcIn
//   for( int layer=1; layer<=NumOfLayersBcIn; ++layer ){
//     const DCRHitContainer &cont=rawData->GetBcInRawHC(layer);
//     int nh=cont.size();
//     for( int i=0; i<nh; ++i ){
//       DCRawHit *rhit=cont[i];
//       //      if(!rhit) continue;
//       DCHit *hit=new DCHit(rhit->PlaneId()+PlOffsBc, rhit->WireId());
//       int nhtdc= rhit->GetTdcSize();
//       for( int j=0; j<nhtdc; ++j ){
// 	hit->SetTdcVal( rhit->GetTdc(j) );
//       }
//       if(!hit) continue;

//       if(hit->CalcObservables())
// 	BcInHC[layer].push_back(hit);
//       else
// 	delete hit;
//     }
//   } 

  // BcOut
//   for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
//     const DCRHitContainer &cont=rawData->GetBcOutRawHC(layer);
//     int nh=cont.size();
//     for( int i=0; i<nh; ++i ){
//       DCRawHit *rhit=cont[i];

//       DCHit *hit=new DCHit(rhit->PlaneId()+PlOffsBc, rhit->WireId());
//       int nhtdc= rhit->GetTdcSize();
//       for( int j=0; j<nhtdc; ++j ){
// 	hit->SetTdcVal( rhit->GetTdc(j) );
//       }
//       if(!hit) continue;

//       if(hit->CalcObservables())
// 	BcOutHC[layer].push_back(hit);
//       else
// 	delete hit;
//     }
//   }
      
  // SdcIn
  for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
    const DCRHitContainer &cont=rawData->GetSdcInRawHC(layer);
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit=cont[i];

      DCHit *hit=new DCHit( rhit->PlaneId(), rhit->WireId() );
      int nhtdc= rhit->GetTdcSize();
      for( int j=0; j<nhtdc; ++j ){
	hit->SetTdcVal( rhit->GetTdc(j) );
      }
      if(!hit) continue; 

      if(hit->CalcObservables())
	SdcInHC[layer].push_back(hit);
      else
	delete hit;
    }
  }
      
  // SdcOut
  for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
    const DCRHitContainer &cont=rawData->GetSdcOutRawHC(layer);
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit=cont[i];

      DCHit *hit=new DCHit( rhit->PlaneId(), rhit->WireId() );
      int nhtdc= rhit->GetTdcSize();
      for( int j=0; j<nhtdc; ++j ){
	hit->SetTdcVal( rhit->GetTdc(j) );
      }
      if(!hit) continue;

      if(hit->CalcObservables())
	SdcOutHC[layer].push_back(hit);
      else
	delete hit;
    }
  }
      
  return true;
}
#endif

bool DCAnalyzer::SetDCHits( int lnum, double lpos )
{
  const std::string funcname = "[DCAnalyzer::SetDCHits]";

  //std::cout << funcname << ", lnum = " << lnum 
  //<< ", lpos = " << lpos << std::endl;
  DCHit *hit=new DCHit( lnum, 0, 0 );
  if ( hit ) {
    double reso = DCConfMan::GetInstance().GetResolution(lnum);
    lpos += GaussRand(0.,reso);
    //std::cout<<"Gauss "<<GaussRand(0.,reso)<<std::endl;
    hit->SetLocalHitPos( lpos );
    hit->SetTiltAngle( DCConfMan::GetInstance().GetTiltAngle(lnum) );
    
    if ( lnum >= PlMinSdcIn && lnum <= PlMaxSdcIn ) {
      SdcInHC[lnum-PlMinSdcIn+1].push_back(hit);    
      //std::cout << "SdcInHC[" << lnum-PlMinSdcIn+1 << "]" << std::endl;
    } else if ( lnum >= PlMinSdcOut && lnum <= PlMaxSdcOut ) {
      SdcOutHC[lnum-PlMinSdcOut+1].push_back(hit);    
      //std::cout << "SdcOutHC[" << lnum-PlMinSdcOut+1 << "]" << std::endl;
    } else if ( lnum >= PlMinBcIn+PlOffsBc && lnum <= PlMaxBcIn+PlOffsBc ) {
      BcInHC[lnum-PlMinBcIn-PlOffsBc+1].push_back(hit);    
      //std::cout << "BcInHC[" << lnum-PlMinBcIn-PlOffsBc+1 << "]" << std::endl;
    } else if ( lnum == 131 || lnum == 132 ) {
      BcOutHC[lnum-118].push_back(hit);    
      //std::cout << "BcInHC[" << lnum-PlMinBcIn-PlOffsBc+1 << "]" << std::endl;
      //std::cout<<"HC="<<BcOutHC[13][0]->GetLayer()<<std::endl;
    } else if ( lnum >= PlMinBcOut+PlOffsBc && lnum <= PlMaxBcOut+PlOffsBc ) {
      BcOutHC[lnum-PlMinBcOut-PlOffsBc+1].push_back(hit);    
      //std::cout << "BcOutHC[" << lnum-PlMinBcOut-PlOffsBc+1 << "]" << std::endl;
    } else {
      std::cerr << funcname << " No such chamber plane " << lnum << std::endl;
      delete hit;
    }

  }
  return true;
}



/* New Added by K.Miwa */
// bool DCAnalyzer::DecodeSimuHits( SimuData *simuData )
// {
//   const std::string funcname = "[DCAnalyzer::DecodeSimuHits]";

//   clearDCHits();
// #if 0
//   // BcIn
//   for( int layer=1; layer<=NumOfLayersBcIn; ++layer ){
//     const DCRHitContainer &cont=rawData->GetBcInRawHC(layer);
//     int nh=cont.size();
//     for( int i=0; i<nh; ++i ){
//       DCRawHit *rhit=cont[i];
//       if(!rhit) continue;
//       DCHit *hit=new DCHit(rhit->PlaneId()+PlOffsBc,
// 			   rhit->WireId(),rhit->Tdc());
//       if(!hit) continue;

//       if(hit->CalcObservables())
// 	BcInHC[layer].push_back(hit);
//       else
// 	delete hit;
//     } /* for( int i... ) */
//  } /* for( int layer... ) */
// #endif
//   // BcOut
//   for( int layer=PlMinBcOut+PlOffsBc; layer<=PlMaxBcOut+PlOffsBc; ++layer ){
//     int nh=simuData->GetDcMulti(layer);
//     for( int i=0; i<nh; ++i ){
//       DCHit *hit=new DCHit(layer, simuData->GetDcHitWire(layer, i), -1);
//       if(!hit) continue;

//       if(hit->CalcObservablesSimulation(simuData->GetDcDrLength(layer, i)))
// 	BcOutHC[layer-(PlMinBcOut+PlOffsBc-1)].push_back(hit);
//       else
// 	delete hit;
//     } /* for( int i... ) */
//   } /* for( int layer... ) */
      
//   // SdcIn
//   for( int layer=PlMinSdcIn; layer<=PlMaxSdcIn; ++layer ){
//     int nh=simuData->GetDcMulti(layer);
//     for( int i=0; i<nh; ++i ){
//       DCHit *hit=new DCHit(layer,simuData->GetDcHitWire(layer, i),-1);

//       if(!hit) continue;

//       if(hit->CalcObservablesSimulation(simuData->GetDcDrLength(layer, i)))
// 	SdcInHC[layer-(PlMinSdcIn-1)].push_back(hit);
//       else
// 	delete hit;
//     } /* for( int i... ) */
//   } /* for( int layer... ) */
      
//   // SdcOut
//   for( int layer=PlMinSdcOut; layer<=PlMaxSdcOut; ++layer ){
//     int nh=simuData->GetDcMulti(layer);
//     for( int i=0; i<nh; ++i ){
//       DCHit *hit=new DCHit(layer,simuData->GetDcHitWire(layer, i),-1); 

//       if(!hit) continue;

//       if(hit->CalcObservablesSimulation(simuData->GetDcDrLength(layer, i)))
// 	SdcOutHC[layer-(PlMinSdcOut-1)].push_back(hit);
//       else
// 	delete hit;
//     } /* for( int i... ) */
//   } /* for( int layer... ) */
      
 
//   return true;
// }

bool DCAnalyzer::TrackSearchBcIn( void )
{
  int ntrack =
    LocalTrackSearch( BcInHC, PPInfoBcIn, NPPInfoBcIn,
		      TrackBcInCol, MinNumOfHitsBcIn );

//   int ntrack =
//     MWPCLocalTrackSearch( BcInHC, TrackBcInCol );

  return true;
}

bool DCAnalyzer::TrackSearchBcOut( void )
{
  //std::cout<<"TrackSearchBcOut before"<<std::endl;
  int ntrack =
    LocalTrackSearch( BcOutHC, PPInfoBcOut, NPPInfoBcOut,
		      TrackBcOutCol, MinNumOfHitsBcOut );
  //std::cout<<"TrackSearchBcOut after"<<std::endl;
  /*
  int n=TrackBcOutCol.size();
  for(int i=0;i<n;++i){
    std::cout<<"i="<<i<<std::endl;
    int m=TrackBcOutCol[i]->GetNHit();
    for(int j=0;j<m;++j){
      DCHit* hit = TrackBcOutCol[i]->GetHit(j);
      std::cout<<" uwaaaaa "<<hit->GetLayer()<<std::endl;
    }
  }
  std::cout<<" "<<std::endl;
  */

  return true;
}

bool DCAnalyzer::TrackSearchSdcIn( void )
{
  int ntrack =
    LocalTrackSearch( SdcInHC, PPInfoSdcIn, NPPInfoSdcIn,
		      TrackSdcInCol, MinNumOfHitsSdcIn );

  return true;
}

#if 0
bool DCAnalyzer::TrackSearchSdcOut( void )
{
  int ntrack =
    SdcOutLocalTrackSearch( SdcOutHC, TrackSdcOutCol );

  return true;
}
#endif

bool DCAnalyzer::TrackSearchSdcOut( void )
{
  int ntrack =
    LocalTrackSearch( SdcOutHC, PPInfoSdcOut, NPPInfoSdcOut,
		      TrackSdcOutCol, MinNumOfHitsSdcOut );

  return true;
}



bool DCAnalyzer::TrackSearchSks( void )
{
  static const std::string funcname = "[DCAnalyzer::TrackSearchSks]";

  clearSksTracks();

  int nIn=TrackSdcInCol.size(), nOut=TrackSdcOutCol.size();
#if 0
  std::cout<<"*********************************************"<<std::endl;
  std::cout << funcname << ": #TracksIn=" << std::setw(3) << nIn
	    << " #TracksOut=" << std::setw(3) << nOut << std::endl;
#endif

  if( nIn==0 || nOut==0 ) return true;

  for( int iIn=0; iIn<nIn; ++iIn ){
    DCLocalTrack *trIn=TrackSdcInCol[iIn];
    //if( !trIn->GoodForTracking() ) continue;
    for( int iOut=0; iOut<nOut; ++iOut ){
      DCLocalTrack *trOut=TrackSdcOutCol[iOut];
      //if( !trOut->GoodForTracking() ) continue;
      SksTrack *tp = new SksTrack( trIn, trOut );
      if(!tp) continue;
      if( tp->doFit() && tp->chisqr()<MaxChiSqrSksTrack )
	SksTrackCol.push_back(tp);
      else
	delete tp;
    }
  }

  partial_sort( SksTrackCol.begin(), SksTrackCol.end(),
		SksTrackCol.end(), SksTrackComp() );

#if 0
 {
   int nn=SksTrackCol.size();
   std::cout << funcname << ": Before Deleting. #Track="
	     << nn << std::endl;
   for( int i=0; i<nn; ++i ){
     SksTrack *tp=SksTrackCol[i];
     std::cout << std::setw(3) << i
	       << " Nitra=" << std::setw(3) << tp->Niteration()
	       << " ChiSqr=" << tp->chisqr()
	       << " P=" << tp->PrimaryMomentum().mag()
	       << " PL(TOF)=" << tp->PathLengthToTOF()
	       << std::endl;
   }
 }
#endif

  return true;
}

bool DCAnalyzer::TrackSearchK18( void )
{
  static const std::string funcname = "[DCAnalyzer::TrackSearchK18]";

  clearK18Tracks();
  int nIn=TrackBcInCol.size(), nOut=TrackBcOutCol.size();
#if 0
  std::cout<<"**************************************"<<std::endl;
  std::cout << funcname << ": #TracksIn=" << std::setw(3) << nIn
            << " #TracksOut=" << std::setw(3) << nOut << std::endl;
#endif

  if( nIn==0 || nOut==0 ) return true;

  //double pK18=ConfMan::GetConfManager()->K18Momentum();
  //  double pK18=0.3*1.15137*4.;//1.38 GeV/c
  //  double pK18 = 1.3;//[GeV/c]
  double pK18 = 1.3;//[GeV/c]

  for( int iIn=0; iIn<nIn; ++iIn ){
    DCLocalTrack *trIn=TrackBcInCol[iIn];
#if 0
    std::cout << "TrackIn  :" << std::setw(2) << iIn
              << " X0=" << trIn->GetX0() << " Y0=" << trIn->GetY0()
              << " U0=" << trIn->GetU0() << " V0=" << trIn->GetV0()
              << std::endl;
#endif
    if( /* !trIn->GoodForTracking() ||*/
        trIn->GetX0()<MinK18InX || trIn->GetX0()>MaxK18InX ||
        trIn->GetY0()<MinK18InY || trIn->GetY0()>MaxK18InY ||
        trIn->GetU0()<MinK18InU || trIn->GetU0()>MaxK18InU ||
        trIn->GetV0()<MinK18InV || trIn->GetV0()>MaxK18InV ) continue;
    for( int iOut=0; iOut<nOut; ++iOut ){
      DCLocalTrack *trOut=TrackBcOutCol[iOut];
#if 0
      std::cout << "TrackOut :" << std::setw(2) << iOut
                << " X0=" << trOut->GetX0() << " Y0=" << trOut->GetY0()
                << " U0=" << trOut->GetU0() << " V0=" << trOut->GetV0()
                << std::endl;
#endif
      if( /*!trOut->GoodForTracking() ||*/
          trOut->GetX0()<MinK18OutX || trOut->GetX0()>MaxK18OutX ||
          trOut->GetY0()<MinK18OutY || trOut->GetY0()>MaxK18OutY ||
          trOut->GetU0()<MinK18OutU || trOut->GetU0()>MaxK18OutU ||
          trOut->GetV0()<MinK18OutV || trOut->GetV0()>MaxK18OutV ) continue;

#if 0
      std::cout << funcname << ": In -> " << trIn->GetChiSquare() 
                << " (" << std::setw(2) << trIn->GetNHit() << ") "
                << "Out -> " << trOut->GetChiSquare()
                << " (" << std::setw(2) << trOut->GetNHit() << ") "
                << std::endl;
#endif

      K18Track *track=new K18Track( trIn, trOut, pK18, BcOutHC[13][0], BcOutHC[14][0]);
      if (track)
	track->CalcMomentumD2U();

      if( track && track->doFit() ) {
      //if( track && track->doFitD2U() ) {
        K18TrackCol.push_back(track);
      } else
        delete track;
    }
  }

#if 0
  std::cout<<"********************"<<std::endl;
  {
    int nn=K18TrackCol.size();
    std::cout << funcname << ": Before sorting. #Track=" 
	      << nn << std::endl;
    for( int i=0; i<nn; ++i ){
      K18Track *tp=K18TrackCol[i];

      std::cout << std::setw(3) << i 
		<< " ChiSqr=" << tp->chisquare()
		<< " Delta=" << tp->Delta()
		<< " P=" << tp->P() << "\n";
      //      std::cout<<"********************"<<std::endl;
      //      std::cout << "In :"
      //             << " X " << tp->Xin() << "(" << tp->TrackIn()->GetX0() << ")"
      //             << " Y " << tp->Yin() << "(" << tp->TrackIn()->GetY0() << ")"
      //             << " U " << tp->Uin() << "(" << tp->TrackIn()->GetU0() << ")"
      //             << " V " << tp->Vin() << "(" << tp->TrackIn()->GetV0() << ")"
      //             << "\n";
      //      std::cout << "Out:"
      //             << " X " << tp->Xout() << "(" << tp->TrackOut()->GetX0() << ")"
      //             << " Y " << tp->Yout() << "(" << tp->TrackOut()->GetY0() << ")"
      //             << " U " << tp->Uout() << "(" << tp->TrackOut()->GetU0() << ")"
      //             << " V " << tp->Vout() << "(" << tp->TrackOut()->GetV0() << ")"
      //             << std::endl;
    }
  }
#endif

  partial_sort( K18TrackCol.begin(), K18TrackCol.end(),
                K18TrackCol.end(), K18TrackComp() );

  return true;
}





void DCAnalyzer::clearDCHits( void )
{
  for( int l=0; l<=NumOfLayersBcIn; ++l ){
    for_each( BcInHC[l].begin(),  BcInHC[l].end(),  DeleteObject() );
    BcInHC[l].clear();
  }

  for( int l=0; l<=NumOfLayersBcOut; ++l ){
    for_each( BcOutHC[l].begin(), BcOutHC[l].end(), DeleteObject() );
    BcOutHC[l].clear();
  }

  for( int l=0; l<=NumOfLayersSdcIn; ++l ){
    for_each( SdcInHC[l].begin(),  SdcInHC[l].end(),  DeleteObject() );
    SdcInHC[l].clear();
  }

  for( int l=0; l<=NumOfLayersSdcOut; ++l ){
    for_each( SdcOutHC[l].begin(), SdcOutHC[l].end(), DeleteObject() );
    SdcOutHC[l].clear();
  }
}


void DCAnalyzer::clearTracksBcIn( void )
{
  for_each( TrackBcInCol.begin(), TrackBcInCol.end(), DeleteObject() );
  TrackBcInCol.clear();
}

void DCAnalyzer::clearTracksBcOut( void )
{
  for_each( TrackBcOutCol.begin(), TrackBcOutCol.end(), DeleteObject() );
  TrackBcOutCol.clear();
}

void DCAnalyzer::clearTracksSdcIn( void )
{
  for_each( TrackSdcInCol.begin(), TrackSdcInCol.end(), DeleteObject() );
  TrackSdcInCol.clear();
}

void DCAnalyzer::clearTracksSdcOut( void )
{
  for_each( TrackSdcOutCol.begin(), TrackSdcOutCol.end(), DeleteObject() );
  TrackSdcOutCol.clear();
}


void DCAnalyzer::clearSksTracks( void )
{
  for_each( SksTrackCol.begin(), SksTrackCol.end(), DeleteObject() );
  SksTrackCol.clear();
}

void DCAnalyzer::clearK18Tracks( void )
{
  for_each( K18TrackCol.begin(), K18TrackCol.end(), DeleteObject() );
  K18TrackCol.clear();
}

