/*
  UserFieldMan.cc
*/

#include "FieldMan.hh"
#include "ThreeVector.hh"
#include "RootHelper.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <string>

bool ProcessOneEvent( void );
void DefineHistograms( const char * );
void InitializeEvent( void );

#ifndef MaxHits 
#define MaxHits 22
#endif

struct Event{
  double PosX;
  double PosY;
  double PosZ;
  double Bx;
  double By;
  double Bz;
};
static Event event;


int main( int argc, char **argv )
{
  if( argc<3 ) {
    std::cerr << "Usage: " << argv[0]
              << " fieldMap RootFile" << std::endl;
    exit(-1);
  }

  FieldMan::GetInstance().Initialize(argv[1]);
  DefineHistograms(argv[3]);

  //Tree
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  InitializeEvent();


  int Nx,Ny,Nz;




  int ifl=1;

  std::cout.setf( std::ios::fixed );
  std::cout.precision(4);

  while(ifl>0){
  double x,y,z;
  ThreeVector B=
    FieldMan::GetInstance().GetField(ThreeVector(x,y,z));

  std::cout << "(x, y, z)[mm]? :";
  std::cin >> x >> y >> z;
  std::cout << "B = ( "
	    << std::setw(12) << B.x() << ","
	    << std::setw(12) << B.y() << ","
	    << std::setw(12) << B.z() << " )"
	    << std::endl;
  
  std::cout << "## ";
  std::cin >> ifl;
  }
  gFile->Write();
  gFile->Close();

  return 0;
}

void InitializeEvent( void )
{
//   //SksTrack
//   event.Hits  = -1;
//   event.chisqr= -999.0;
//   for( int it=0; it<MaxHits; it++){
//    event.DC1Posl[it] = -999.0; 
//   }
}

void DefineHistograms( const char *filename )
{
//   new TFile( filename, "recreate" );

//   HB1( 101, "stat",  20, 0., 20. );
//   HB1( 102, "status",  10, 0., 10. );
//   HB1( 103, "chisqr", 500, 0., 100. );

//   for( int i=1; i<11; ++i ){
//     std::ostringstream title1, title2, title3;
//     std::ostringstream title4, title5, title6;
    
//     title1 << "HitPos [mm] Layer#" << i;
//     HB1( 200+i, title1.str().c_str(), NBin1HPDcIn[i], 
// 	 MinHPDcIn[i], MaxHPDcIn[i] );
//  }
//   //Tree
//   HBTree("tree","tree of Sks");
//   TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

//   //Primary
//   tree->Branch("Hits",   &event.Hits,  "Hits/I");
//   tree->Branch("chisqr", &event.chisqr,"chisqr/D");

//   tree->Branch("DC1Hits", &event.DC1Hits, "DC1Hits/I");
//   tree->Branch("DC2Hits", &event.DC2Hits, "DC2Hits/I");
//   tree->Branch("DC1Posl", event.DC1Posl,"DC1Posl[DC1Hits]/D");
//   tree->Branch("DC1Resd", event.DC1Resd,"DC1Resd[DC1Hits]/D");

}
