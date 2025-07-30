// -*- C++ -*-

#include "RCHit.hh"

#include <iomanip>

#include <G4VVisManager.hh>
#include <G4Circle.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>
#include <G4ios.hh>

G4Allocator<RCHit> RCHitAllocator;

//_____________________________________________________________________________
RCHit::RCHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
RCHit::~RCHit()
{
#if 0
  VHitInfo::Print();
#endif 
}

//_____________________________________________________________________________
void
RCHit::Draw()
{
}

//_____________________________________________________________________________
void
RCHit::Print()
{
}
