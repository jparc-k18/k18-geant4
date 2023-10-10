/*
  "S2STransportation.hh"
  
  Modified by Toshi Gogami , 11Nov2014
*/

#ifndef S2STransportation_h
#define S2STransportation_h 1

#include "G4Transportation.hh"

class S2STransportation : public G4Transportation
{
public:
  S2STransportation() 
    :G4Transportation()
  {}
  ~S2STransportation()
  {}
  
  G4double AlongStepGetPhysicalInteractionLength( const G4Track &track,
                                                  G4double previousStepSize, 
                                                  G4double currentMinimumStep,
                                                  G4double &currerntSafety, 
                                                  G4GPILSelection *selection );
};

#endif
