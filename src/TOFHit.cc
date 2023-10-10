// -*- C++ -*-

#include "TOFHit.hh"

#include <iomanip>

#include <G4VVisManager.hh>
#include <G4Circle.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>
#include <G4ios.hh>

G4Allocator<TOFHit> TOFHitAllocator;

//_____________________________________________________________________________
TOFHit::TOFHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
TOFHit::~TOFHit()
{
}

//_____________________________________________________________________________
void
TOFHit::Draw()
{
}

//_____________________________________________________________________________
void
TOFHit::Print()
{
}
