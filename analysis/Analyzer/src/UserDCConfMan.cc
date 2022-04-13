/*
  DCConfMan.cc
*/

#include "DCConfMan.hh"
#include "ThreeVector.hh"

#include <iostream>

int main( int argc, char **argv )
{
  if( argc<2 ){
    std::cerr << "Usage: " << argv[0] << " confFile" << std::endl;
    exit(-1);
  }

  DCConfMan::GetInstance().Initialize(argv[1]);

  int lnum=-1;

  std::cout << "LayerNumber ? : ";
  std::cin >> lnum;

  while(lnum>=0){
    double z=DCConfMan::GetInstance().GetLocalZ(lnum);
    double resol=DCConfMan::GetInstance().GetResolution(lnum);
    double ta=DCConfMan::GetInstance().GetTiltAngle(lnum);
    double ra1=DCConfMan::GetInstance().GetRotAngle1(lnum);
    double ra2=DCConfMan::GetInstance().GetRotAngle2(lnum);
    ThreeVector pos=DCConfMan::GetInstance().GetGlobalPosition(lnum);
    ThreeVector normal=DCConfMan::GetInstance().NormalVector(lnum);
    ThreeVector unit=DCConfMan::GetInstance().UnitVector(lnum);
    std::cout << lnum << " " << pos << " " << z << " " << resol
	      << " " << ta << " " << ra1 << " " << ra2 << std::endl;
    std::cout << lnum << " Normal: " << normal << " > " << normal.mag()
	      << " Unit: " << unit << " > " << unit.mag()
	      << std::endl;

    std::cout << "LayerNumber ? : ";
    std::cin >> lnum;
  }
 
}
