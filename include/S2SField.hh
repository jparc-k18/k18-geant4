/*
  "S2SField.hh"
  
  Modified by Toshi Gogami , 21Nov2014
*/

#ifndef S2SField_h 
#define S2SField_h 1

#include "S2SFieldMap.hh"
#include "G4MagneticField.hh"

#include <vector>

class SimpleFieldElement;

class S2SField : public G4MagneticField
{
public:
  S2SField();
  explicit S2SField( const std::string &FieldMapName, 
		     double scaleFactor=1.0 );
  ~S2SField();

private:
  S2SField( const S2SField & );
  S2SField & operator = ( const S2SField & );

public:
  void GetFieldValue( const double Point[4], double *Bfield ) const;
  void cleanupSimpleElementList( void );
  void AddSimpleElement( SimpleFieldElement *elem );

private:
  S2SFieldMap fMap;

  bool InMagnet(double *pos) const;

  typedef std::vector <SimpleFieldElement *> FMContainer;
  typedef std::vector <SimpleFieldElement *>
  ::const_iterator FMIterator;

  FMContainer elemList_;
};

#endif

