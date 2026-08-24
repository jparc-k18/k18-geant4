// -*- C++ -*-

#ifndef K18_BEAMLINE_DETECTOR_BUILDER_HH
#define K18_BEAMLINE_DETECTOR_BUILDER_HH 1

#include <G4Types.hh>

class G4LogicalVolume;

// Builds the QQDQQ magnets, beam ducts, tracking chambers, and hodoscopes in
// the S-2S global coordinate system.
class K18BeamlineDetectorBuilder
{
public:
  K18BeamlineDetectorBuilder();
  ~K18BeamlineDetectorBuilder();

  void Construct(G4LogicalVolume* world_lv, G4bool check_overlaps);
};

#endif
