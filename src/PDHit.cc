// -*- C++ -*-

#include "PDHit.hh"

#include <iomanip>

#include <G4VVisManager.hh>
#include <G4Circle.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>
#include <G4ios.hh>

G4Allocator<PDHit> PDHitAllocator;

//_____________________________________________________________________________
PDHit::PDHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
PDHit::~PDHit()
{
#if 0
  VHitInfo::Print();
#endif 
}

//_____________________________________________________________________________
void
PDHit::Draw()
{
}

//_____________________________________________________________________________
void
PDHit::Print()
{
}
