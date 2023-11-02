/*
  QMag.cc

  2004/7/5  T.Takahashi

*/

#include "QMag.hh"

#include <cmath>

QMag::QMag( const char *name, const ThreeVector &pos,
	    double ta, double ra1, double ra2,
	    double fieldExtX, double fieldExtY, double fieldExtZ,
	    double fValTeslaPerMM )
  : FieldElements( name, pos, ta, ra1, ra2 ),
    extX_(fieldExtX), extY_(fieldExtY), extZ_(fieldExtZ),
    fVal_(fValTeslaPerMM)
{}

ThreeVector QMag::GetField( const ThreeVector &gPos ) const
{
  ThreeVector lPos = Global2LocalPos( gPos );
  ThreeVector Field(0.,0.,0);
  if( fabs(lPos.x())<=extX_ && fabs(lPos.y())<=extY_ &&
      fabs(lPos.z())<=extZ_ ){
    ThreeVector lField( -fVal_*lPos.y(), -fVal_*lPos.x(), 0.0 );
    Field = Local2GlobalDir( lField );
  }
  return Field;
}


bool QMag::ExistField( const ThreeVector &gPos ) const
{
  ThreeVector lPos = Global2LocalPos( gPos );
  if( fabs(lPos.x())<=extX_ && fabs(lPos.y())<=extY_ &&
      fabs(lPos.z())<=extZ_ )
    return true;
  else
    return false;
}
    
FldElemReg QMag::checkRegion( const ThreeVector &gPos,
			      double Tolerance ) const
{
  double d=fabs(Tolerance);
  ThreeVector lPos = Global2LocalPos( gPos );
  if( fabs(lPos.x())>=extX_+d ||
      fabs(lPos.y())>=extY_+d ||
      fabs(lPos.z())>=extZ_+d )
    return FEROutside;
  else if( fabs(lPos.x())<=extX_-d && 
	   fabs(lPos.y())<=extY_-d &&
	   fabs(lPos.z())<=extZ_-d )
    return FERInside;
  else
    return FERSurface;
}
