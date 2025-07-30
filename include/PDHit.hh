// -*- C++ -*-

#ifndef PDHit_h
#define PDHit_h 1

#include <G4VHit.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4Allocator.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class PDHit : public G4VHit, public VHitInfo
{
public:
  PDHit(const G4String& name, G4Step* step);
  virtual ~PDHit();

private:
  PDHit(const PDHit &);
  PDHit& operator =(const PDHit&);

public:
  void* operator new(size_t size);
  void operator delete(void * aHit);
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
using PDHitsCollection = G4THitsCollection<PDHit>;
extern G4Allocator<PDHit> PDHitAllocator;

//_____________________________________________________________________________
inline void*
PDHit::operator new(size_t)
{
  return PDHitAllocator.MallocSingle();
}

//_____________________________________________________________________________
inline void
PDHit::operator delete(void *aHit)
{
  PDHitAllocator.FreeSingle(static_cast<PDHit*>(aHit));
}

#endif
