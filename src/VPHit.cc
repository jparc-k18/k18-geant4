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
VPHit::VPHit()
  : layerID_(0),
    wireID_(0),
    time_(-100.0),
    trackNo_(0)
{
}
