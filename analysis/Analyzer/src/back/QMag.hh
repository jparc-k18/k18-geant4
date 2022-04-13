/*
  QMag.hh
*/

#ifndef QMag_h 
#define QMag_h 1

#include "FieldElements.hh"

class QMag : public FieldElements
{
public:
  QMag( const char *name, const ThreeVector &pos,
	double ta, double ra1, double ra2,
	double fieldExtX, double fieldExtY, double fieldExtZ,
	double fValTeslaPerMM );
  ~QMag(){}

private:
  double extX_, extY_, extZ_;
  double fVal_;
public:
  ThreeVector GetField( const ThreeVector &gPos ) const;
  bool ExistField( const ThreeVector &gPos ) const;
  FldElemReg checkRegion( const ThreeVector &gPos, 
			  double Tolerance ) const;
};

#endif
