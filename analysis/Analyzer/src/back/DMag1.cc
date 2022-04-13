/*
  DMag1.cc

  2004/7/12    T.Takahashi

*/

#include "DMag1.hh"

#include <cmath>

DMag1::DMag1( const char *name, const ThreeVector &pos,
	      double ta, double ra1, double ra2,
	      double halfGap, double fieldTesla )
  : FieldElements( name, pos, ta, ra1, ra2 ),
    halfGap_(halfGap), fVal_(fieldTesla)
{}

ThreeVector DMag1::GetField( const ThreeVector &gPos ) const
{
  ThreeVector lPos = Global2LocalPos( gPos );
  ThreeVector Field( 0., 0., 0. );

  if( fabs(lPos.y())<=halfGap_ && fabs(lPos.x())<=350.
      && fabs(lPos.z())<=278.-0.26*lPos.x() ){
    ThreeVector flocal( 0., fVal_, 0. );
    Field = Local2GlobalDir( flocal );
  }
  return Field;
}

bool DMag1::ExistField( const ThreeVector &gPos ) const
{
  ThreeVector lPos = Global2LocalPos( gPos );
  if( fabs(lPos.y())<=halfGap_ && fabs(lPos.x())<=350.
      && fabs(lPos.z())<=278.-0.26*lPos.x() )
    return true;
  else
    return false;
}

FldElemReg DMag1::checkRegion( const ThreeVector &gPos,
			       double Tolerance ) const
{
  double d=fabs(Tolerance), dp=d*0.26/sqrt(1.+0.26*0.26);
  ThreeVector lPos = Global2LocalPos( gPos );
  if( fabs(lPos.y())>=halfGap_+d ||
      fabs(lPos.x())>=350.+d ||
      fabs(lPos.z())>=278.-0.26*lPos.x()+dp )
    return FEROutside;
  else if( fabs(lPos.y())<=halfGap_-d &&
	   fabs(lPos.x())<=350.-d &&
	   fabs(lPos.z())<=278.-0.26*lPos.x()-dp )
    return FERInside;
  else
    return FERSurface;
}
