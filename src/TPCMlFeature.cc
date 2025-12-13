// -*- C++ -*-

#include "TPCMlFeature.hh"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <TLorentzVector.h>
#include <G4ParticleTable.hh>

#include "ConfMan.hh"
#include "S2SAnaManager.hh"

namespace
{
constexpr double kInvalidMissingMass = -9999.0;

double GetParticleMass(const G4String& name)
{
  auto* table = G4ParticleTable::GetParticleTable();
  auto* particle = table ? table->FindParticle(name) : nullptr;
  return particle ? particle->GetPDGMass() : 0.0;
}
}


MlTrackFeature CalculateMlFeature(std::vector<TpcMlHit>& hits, double truncateRate)
{
  std::sort(hits.begin(), hits.end(),
            [](const TpcMlHit& a, const TpcMlHit& b) { return a.time < b.time; });

  /* collect and calculate de/dx */
  const auto pDir = hits.front().mom.unit();
  std::vector<double> dedxSamples;
  dedxSamples.reserve(hits.size());
  std::vector<double> dedxTrunc;
  for (const auto& h : hits) {
    if (h.edep <= 0.0) continue; // dedx is stored directly in edep [MeV/cm]
    dedxSamples.push_back(h.edep);
  }

  /* calculate de/dx truncated mean */
  double dedxVal = 0.0;
  if (!dedxSamples.empty()) {
    std::sort(dedxSamples.begin(), dedxSamples.end());
    const auto nUse = std::max<std::size_t>
      (1u, static_cast<std::size_t>(dedxSamples.size() * (1.0 - truncateRate)));
    double sum = 0.0;
    for (std::size_t i = 0; i < nUse; ++i)
      sum += dedxSamples[i];
    dedxTrunc.assign(dedxSamples.begin(), dedxSamples.begin() + nUse);
    dedxVal = sum / static_cast<double>(nUse);
  }

  /* collect direction */
  MlTrackFeature tf;
  tf.ux    = pDir.x();
  tf.uy    = pDir.y();
  tf.uz    = pDir.z();
  const auto pdg = hits.front().pdg;
  const auto parent = hits.front().parentId;
  if (std::abs(pdg) == 2212) {
    tf.pidCode = 1; // decay proton
  } else if (std::abs(pdg) == 211) {
    tf.pidCode = (parent == 0) ? 0 : 2; // scat pi if parent=0, decay pi otherwise
  } else {
    tf.pidCode = -1;
  }
  tf.dedx  = dedxVal;
  return tf;
}

double CalculateMissingMass(double p_beam,
                            double p_pi,
                            double theta_pi,
                            double phi_pi)
{
  if (p_beam <= 0.0 || p_pi <= 0.0)
    return kInvalidMissingMass;

  const double M_Kaon = GetParticleMass("kaon-");
  const double M_Deuteron = GetParticleMass("deuteron");
  const double M_PiM = GetParticleMass("pi-");
  if (M_Kaon <= 0.0 || M_Deuteron <= 0.0 || M_PiM <= 0.0)
    return kInvalidMissingMass;

  TLorentzVector beam_lv(0., 0., p_beam, std::sqrt(p_beam * p_beam + M_Kaon * M_Kaon));
  TLorentzVector target_lv(0., 0., 0., M_Deuteron);

  const double sinTheta = std::sin(theta_pi);
  const double px = p_pi * sinTheta * std::sin(phi_pi);
  const double py = p_pi * sinTheta * std::cos(phi_pi);
  const double pz = p_pi * std::cos(theta_pi);
  TLorentzVector pi_lv;
  pi_lv.SetPxPyPzE(px, py, pz, std::sqrt(p_pi * p_pi + M_PiM * M_PiM));

  const TLorentzVector X = beam_lv + target_lv - pi_lv;
  return X.M(); // [MeV]
}

void ResetMlFeatures(Event& event,
                     const ConfMan& confMan,
                     double qnan,
                     std::size_t trackCount)
{
  static const float qnanf = static_cast<float>(qnan);
  const auto generator = confMan.Get<G4int>("Generator");
  switch (generator) {
  case 9001: event.label = 1; break;
  case 9002: event.label = 2; break;
  case 9003: event.label = 3; break;
  default:   event.label = confMan.Get<G4int>("label"); break;
  }
  if (event.mlUx.size() != trackCount) {
    event.mlUx.assign(trackCount, qnanf);
    event.mlUy.assign(trackCount, qnanf);
    event.mlUz.assign(trackCount, qnanf);
    event.mlDedx.assign(trackCount, qnanf);
    // event.mlPid.assign(trackCount, -1);
  } else {
    std::fill(event.mlUx.begin(), event.mlUx.end(), qnanf);
    std::fill(event.mlUy.begin(), event.mlUy.end(), qnanf);
    std::fill(event.mlUz.begin(), event.mlUz.end(), qnanf);
    std::fill(event.mlDedx.begin(),  event.mlDedx.end(),  qnanf);
    // std::fill(event.mlPid.begin(), event.mlPid.end(), -1);
  }
  event.mlMM = qnanf;
}
