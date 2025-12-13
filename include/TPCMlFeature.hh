// -*- C++ -*-

#ifndef TPC_ML_FEATURE_HH
#define TPC_ML_FEATURE_HH

#include <cstddef>
#include <vector>

#include <globals.hh>
#include <G4ThreeVector.hh>

class ConfMan;
struct Event;

struct TpcMlHit
{
  G4int trackId;
  G4int pdg;
  G4int parentId;
  double time;
  double edep;
  G4ThreeVector pos;
  G4ThreeVector mom;
};

struct MlTrackFeature
{
  double ux;
  double uy;
  double uz;
  G4int pidCode; // 0: scat pi (parent=0), 1: decay proton, 2: decay pi, -1: unknown
  double dedx;
  bool operator>(const MlTrackFeature& other) const { return dedx > other.dedx; }
};

MlTrackFeature CalculateMlFeature(std::vector<TpcMlHit>& hits, double truncateRate);

void ResetMlFeatures(Event& event,
                     const ConfMan& confMan,
                     double qnan,
                     std::size_t trackCount);

double CalculateMissingMass(double p_beam,
                            double p_pi,
                            double theta_pi,
                            double phi_pi);

// Bethe-Bloch (P10) -> dE/dx [MeV/cm] without ADC conversion
double CalculateTPCDedx(double massMeV, double beta);

#endif
