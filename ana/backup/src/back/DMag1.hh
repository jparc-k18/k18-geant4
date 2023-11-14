/*
  DMag1.hh

  2004/7/12    T.Takahashi

*/

#ifndef DMag1_h

#define DMag1_h 1

#include "FieldElements.hh"

class DMag1 : public FieldElements
{
public:
  DMag1( const char *name, const ThreeVector &pos,
	 double ta, double ra1, double ra2,
	 double halfGap, double fieldTesla );
  ~DMag1(){}

private:
  double halfGap_;
  double fVal_;

public:
  ThreeVector GetField( const ThreeVector &gPos ) const;
  bool ExistField( const ThreeVector &gPos ) const;
  FldElemReg checkRegion( const ThreeVector &gPos, 
                          double Tolerance ) const;
  
};

#endif
