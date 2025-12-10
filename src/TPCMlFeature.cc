// -*- C++ -*-

#include "TPCMlFeature.hh"

#include <algorithm>
#include <cmath>

#include "ConfMan.hh"
#include "S2SAnaManager.hh"

namespace
{

double DensityEffectCorrection(double betagamma, const double* par)
{
  const double constant = 2.0 * std::log(10.0);
  const double X = std::log10(betagamma);
  if (X <= par[2])
    return par[5] * std::pow(10.0, 2.0 * (X - par[2]));
  else if (par[2] < X && X < par[3])
    return constant * X - par[4] + par[0] * std::pow((par[3] - X), par[1]);
  else
    return constant * X - par[4];
}
}

double CalculateTPCDedx(double massMeV, double beta)
{
  double rho = std::pow(10., -3) * (0.9 * 1.662 + 0.1 * 0.6672); // g/cm3
  double ZoverA = 17.2 / 37.6;                                    // mol/g
  double I = 0.9 * 188.0 + 0.1 * 41.7;                            // eV
  double density_effect_par[6] = {
    0.9 * 0.19714 + 0.1 * 0.09253,
    0.9 * 2.9618 + 0.1 * 3.6257,
    0.9 * 1.7635 + 0.1 * 1.6263,
    0.9 * 4.4855 + 0.1 * 3.9716,
    0.9 * 11.9480 + 0.1 * 9.5243,
    0.0};

  constexpr double Z = 1.0;
  constexpr double me = 0.5109989461; // MeV
  constexpr double K = 0.307075;      // MeV cm2 mol-1

  const double constant = rho * K * ZoverA; // MeV/cm
  const double beta2 = beta * beta;
  const double gamma2 = 1.0 / (1.0 - beta2);
  const double MeVToeV = std::pow(10., 6);
  const double massRatio = me / massMeV;
  const double Wmax = 2 * me * beta2 * gamma2 /
    ((massRatio + 1.) * (massRatio + 1.) + 2 * massRatio * (std::sqrt(gamma2) - 1));
  const double delta = DensityEffectCorrection(std::sqrt(beta2 * gamma2),
                                               density_effect_par);
  const double dedx = constant * Z * Z / beta2 *
    (0.5 * std::log(2 * me * beta2 * gamma2 * Wmax * MeVToeV * MeVToeV / (I * I))
     - beta2 - 0.5 * delta);
  return dedx; // MeV/cm
}

MlTrackFeature CalculateMlFeature(std::vector<TpcMlHit>& hits, double truncateRate)
{
  std::sort(hits.begin(), hits.end(),
            [](const TpcMlHit& a, const TpcMlHit& b) { return a.time < b.time; });

  /* collect and calculate de/dx */
  const auto pDir = hits.front().mom.unit();
  std::vector<double> dedxSamples;
  dedxSamples.reserve(hits.size());
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
    event.mlPdg.assign(trackCount, -1);
    event.mlDedx.assign(trackCount, qnanf);
  } else {
    std::fill(event.mlUx.begin(), event.mlUx.end(), qnanf);
    std::fill(event.mlUy.begin(), event.mlUy.end(), qnanf);
    std::fill(event.mlUz.begin(), event.mlUz.end(), qnanf);
    std::fill(event.mlPdg.begin(), event.mlPdg.end(), -1);
    std::fill(event.mlDedx.begin(),  event.mlDedx.end(),  qnanf);
  }
}
