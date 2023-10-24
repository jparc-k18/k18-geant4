// -*- C++ -*-

#include "VPHit.hh"

#include <iomanip>

#include <G4VVisManager.hh>
#include <G4Circle.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>
#include <G4ios.hh>

G4Allocator<VPHit> VPHitAllocator;

//_____________________________________________________________________________
VPHit::VPHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
VPHit::~VPHit()
{
}

//_____________________________________________________________________________
void
VPHit::Draw()
{
}

//_____________________________________________________________________________
void
VPHit::Print()
{
}
