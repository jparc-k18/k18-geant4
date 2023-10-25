// -*- C++ -*-

#ifndef ACHIT_HH
#define ACHIT_HH 1

#include <G4VHit.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4Allocator.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class ACHit : public G4VHit, public VHitInfo
{
public:
  ACHit(const G4String& name, G4Step* step);
  virtual ~ACHit();

private:
  ACHit(const ACHit&);
  ACHit& operator =(const ACHit&);

public:
  void* operator new(size_t size);
  void operator delete(void* aHit);
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
using ACHitsCollection = G4THitsCollection<ACHit>;
extern G4Allocator<ACHit> ACHitAllocator;

//_____________________________________________________________________________
inline void*
ACHit::operator new(size_t)
{
  return static_cast<void*>(ACHitAllocator.MallocSingle());
}

//_____________________________________________________________________________
inline void
ACHit::operator delete(void* aHit)
{
  ACHitAllocator.FreeSingle(static_cast<ACHit*>(aHit));
}

#endif
