/*
  "S2STransportation.cc"
  
  Modified by Toshi Gogami , 11Nov2014
*/

#include "S2STransportation.hh"
#include "G4SystemOfUnits.hh"

const double MinimumStep = 5.0*mm;

G4double S2STransportation::
AlongStepGetPhysicalInteractionLength( const G4Track &track,
                                       G4double previousStepSize, 
                                       G4double currentMinimumStep,
                                       G4double &currentSafety, 
                                       G4GPILSelection *selection ){
  
  if( DoesGlobalFieldExist() &&currentMinimumStep>MinimumStep ){
    currentMinimumStep=MinimumStep;
  }
  
  return G4Transportation::
    AlongStepGetPhysicalInteractionLength( track, 
					   previousStepSize, 
                                           currentMinimumStep,
					   currentSafety, 
                                           selection );
  
}
