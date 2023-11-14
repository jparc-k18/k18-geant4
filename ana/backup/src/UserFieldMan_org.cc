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

#ifndef MaxHits 
#define MaxHits 22
#endif

struct Event{
  double Bx;
  double By;
  double Bz;
};
static Event event;


int main( int argc, char **argv )
{
  if( argc<2 ) {
    std::cerr << "Usage: " << argv[0]
              << " fieldMap" << std::endl;
    exit(-1);
  }

  FieldMan::GetInstance().Initialize(argv[1]);

  int ifl=1;

  std::cout.setf( std::ios::fixed );
  std::cout.precision(4);
 
  while(ifl>0){
    double x,y,z;
    std::cout << "(x, y, z)[mm]? :";
    std::cin >> x >> y >> z;
    ThreeVector B=
      FieldMan::GetInstance().GetField(ThreeVector(x,y,z));
    std::cout << "B = ( "
	      << std::setw(12) << B.x() << ","
	      << std::setw(12) << B.y() << ","
	      << std::setw(12) << B.z() << " )"
	      << std::endl;
    cout<<"Nx="<< B.size() << endl;

    std::cout << "## ";
    std::cin >> ifl;
  }

 return 0;
 }
