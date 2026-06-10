// -*- C++ -*-

#include "GeHit.hh"

G4Allocator<GeHit> GeHitAllocator;

GeHit::GeHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step, true)
{
}

GeHit::~GeHit()
{
}

void
GeHit::Draw()
{
}

void
GeHit::Print()
{
}
