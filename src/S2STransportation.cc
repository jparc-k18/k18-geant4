// -*- C++ -*-

#include "S2STransportation.hh"

#include <G4SystemOfUnits.hh>

namespace
{
const double MinimumStep = 5.0*mm;
}

G4double S2STransportation::
AlongStepGetPhysicalInteractionLength(const G4Track &track,
                                      G4double previousStepSize,
                                      G4double currentMinimumStep,
                                      G4double &currentSafety,
                                      G4GPILSelection *selection){
  if(fAnyFieldExists && currentMinimumStep>MinimumStep){
    currentMinimumStep = MinimumStep;
  }

  return G4Transportation::
    AlongStepGetPhysicalInteractionLength(track,
                                          previousStepSize,
                                          currentMinimumStep,
                                          currentSafety,
                                          selection);
}
