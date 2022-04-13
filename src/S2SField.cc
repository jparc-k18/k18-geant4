/*
  "S2SField.cc"
  
  Modified by Toshi Gogami , 21Nov2014
*/

#include "S2SField.hh"
#include "SimpleFieldElement.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

S2SField::S2SField( const std::string &FieldMapName, double scaleFactor )
  : fMap( FieldMapName.c_str(),scaleFactor )
{
  fMap.Initialize();
}

S2SField::~S2SField()
{}

void S2SField::GetFieldValue( const double Point[4], 
			      double *Bfield ) const
{
  double X[3];
  //   X[0]=Point[0]/cm; X[1]=Point[1]/cm; X[2]=Point[2]/cm;
  X[0]=Point[0]; X[1]=Point[1]; X[2]=Point[2];
  G4ThreeVector x(X[0],X[1],X[2]);
  G4ThreeVector MagOrg(600/tan(55*deg),-600,0);
  
  if( InMagnet(X) ){
    
    //Map2 region -> conversion
    if( X[1]>-tan(70*deg)*X[0] ){
      x = x - MagOrg;
      x.rotateZ(-70*deg);
      x = MagOrg + x;
      X[0]=x.x(); X[1]=x.y(); X[2]=x.z();
    }
    
    if( fMap.GetFieldValue( X, Bfield ) ){
      Bfield[0] *= tesla;
      Bfield[1] *= tesla;
      Bfield[2] *= tesla;
    }
    
    if( X[1]>-tan(70*deg)*X[0] ){
      G4ThreeVector tmpB(Bfield[0],Bfield[1],Bfield[2]);
      tmpB.rotateZ(70*deg);
      Bfield[0] = tmpB.x();
      Bfield[1] = tmpB.y();
      Bfield[2] = tmpB.z();
    }
    
  }
  else{
    Bfield[0]=Bfield[1]=Bfield[2]=0.0;
  }
  
  
  G4ThreeVector gPos( Point[0], Point[1], Point[2] );
  G4ThreeVector B( 0., 0., 0. );
  FMIterator end=elemList_.end();
  for( FMIterator itr=elemList_.begin(); itr!=end; ++itr ){
    if( (*itr)->ExistMagneticField() )
      B += (*itr)->GetMagneticField( gPos );
  }
  
  Bfield[0] += B.x(); Bfield[1] += B.y(); Bfield[2] += B.z();
  
#if 0
  G4cout << "X=(" << X[0] << "," << X[1] << "," << X[2] << ") "
	 << "B=(" << Bfield[0]/tesla << "," << Bfield[1]/tesla
	 << "," << Bfield[2]/tesla << ")" << G4endl;
#endif
}

bool S2SField::InMagnet(double *pos) const
{

  G4ThreeVector gPos(pos[0], pos[1], pos[2]);
  G4ThreeVector MagOrg(600/tan(55*deg),-600,0);

  if( pos[1]>-tan(70*deg)*pos[0] ){
    gPos = gPos - MagOrg;
    gPos.rotateZ(-70*deg);
    gPos = MagOrg + gPos;
  }
  bool bx,by,bz;
  if(-5500<gPos.x() && gPos.x()<5400) bx=true; else bx=false;
  if(-600<gPos.y() && gPos.y()<900) by=true; else by=false;
  if(-320<gPos.z() && gPos.z()<320) bz=true; else bz=false;

  if(bx&&by&&bz) return true;
  else return false;
}


void S2SField::cleanupSimpleElementList( void )
{
  elemList_.clear();
}

void S2SField::AddSimpleElement( SimpleFieldElement *elem )
{
  elemList_.push_back( elem );
}
