// -*- C++ -*-

#include "DCHit.hh"

#include <iomanip>

#include <G4VVisManager.hh>
#include <G4Circle.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>
#include <G4ios.hh>

G4Allocator<DCHit> DCHitAllocator;

//_____________________________________________________________________________
DCHit::DCHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step, true)
{
}

//_____________________________________________________________________________
DCHit::~DCHit()
{
}

//_____________________________________________________________________________
void
DCHit::Draw()
{
}

//_____________________________________________________________________________
void
DCHit::Print()
{
}
