// -*- C++ -*-

#include "TPCMlFeature.hh"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "ConfMan.hh"
#include "S2SAnaManager.hh"


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
}
