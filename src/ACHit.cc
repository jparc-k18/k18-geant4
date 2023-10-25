// -*- C++ -*-

#include "ACHit.hh"

#include <iomanip>

#include <G4VVisManager.hh>
#include <G4Circle.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>
#include <G4ios.hh>

G4Allocator<ACHit> ACHitAllocator;

//_____________________________________________________________________________
ACHit::ACHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
ACHit::~ACHit()
{
}

//_____________________________________________________________________________
void
ACHit::Draw()
{
}

//_____________________________________________________________________________
void
ACHit::Print()
{
}
