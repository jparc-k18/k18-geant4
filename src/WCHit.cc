// -*- C++ -*-

#include "WCHit.hh"

#include <iomanip>

#include <G4VVisManager.hh>
#include <G4Circle.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>
#include <G4ios.hh>

G4Allocator<WCHit> WCHitAllocator;

//_____________________________________________________________________________
WCHit::WCHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{}


//_____________________________________________________________________________
WCHit::~WCHit()
{
}

//_____________________________________________________________________________
void
WCHit::Draw()
{
}

//_____________________________________________________________________________
void
WCHit::Print()
{
}
