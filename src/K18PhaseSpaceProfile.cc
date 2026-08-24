// -*- C++ -*-

#include "K18PhaseSpaceProfile.hh"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>

#include <G4ios.hh>
#include <Randomize.hh>

#include <TBranch.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TTree.h>

namespace
{
  std::runtime_error ProfileError(const G4String& file_name,
                                  const G4String& detail)
  {
    return std::runtime_error("K18PhaseSpaceProfile " + file_name
                              + ": " + detail);
  }

  G4bool HasRootExtension(const G4String& file_name)
  {
    G4String lower = file_name;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char value){ return std::tolower(value); });
    const G4String extension(".root");
    return lower.size() >= extension.size()
      && lower.compare(lower.size()-extension.size(), extension.size(),
                       extension) == 0;
  }

  TLeaf* RequireScalarLeaf(TTree* tree, const G4String& file_name,
                           const char* name, const char* type,
                           const char* tree_name)
  {
    TBranch* branch = tree ? tree->GetBranch(name) : nullptr;
    TLeaf* leaf = branch ? branch->GetLeaf(name) : nullptr;
    if(!leaf && tree)
      leaf = tree->GetLeaf(name);
    if(!branch || !leaf)
      throw ProfileError(file_name, "tree " + G4String(tree_name)
                         + " is missing branch "
                         + G4String(name));
    if(G4String(leaf->GetTypeName()) != type){
      throw ProfileError(file_name, "branch " + G4String(name)
                         + " has type " + leaf->GetTypeName()
                         + ", expected " + type);
    }
    if(leaf->GetLeafCount() || leaf->GetLenStatic() != 1)
      throw ProfileError(file_name, "branch " + G4String(name)
                         + " must be scalar");
    return leaf;
  }
}

//_____________________________________________________________________________
K18PhaseSpaceProfile::K18PhaseSpaceProfile(
  const G4String& file_name, const G4String& root_tree,
  G4double root_z_mm, G4double root_chisqr_max,
  G4int root_charge, G4double root_p_gev, G4long root_max_tracks)
  : m_file_name(file_name),
    m_root_tree(root_tree),
    m_root_z_mm(root_z_mm),
    m_root_chisqr_max(root_chisqr_max),
    m_root_charge(root_charge),
    m_root_p_gev(root_p_gev),
    m_root_max_tracks(root_max_tracks),
    m_mode(Mode::GaussianMoments),
    m_anchors(),
    m_empirical_anchors()
{
  if(HasRootExtension(file_name)){
    m_mode = Mode::EmpiricalExact;
    LoadRoot();
  } else {
    LoadAscii();
  }
}

//_____________________________________________________________________________
void
K18PhaseSpaceProfile::LoadAscii()
{
  std::ifstream input(m_file_name);
  if(!input.is_open())
    throw ProfileError(m_file_name, "cannot open profile table");

  G4String line;
  G4int line_number = 0;
  while(std::getline(input, line)){
    ++line_number;
    const auto comment = line.find('#');
    if(comment != G4String::npos)
      line.erase(comment);

    std::istringstream stream(line);
    stream >> std::ws;
    if(stream.eof())
      continue;

    std::array<G4double, 16> value{};
    for(std::size_t i=0; i<value.size(); ++i){
      if(!(stream >> value[i])){
        std::ostringstream detail;
        detail << "line " << line_number
               << " is malformed (expected 16 numeric columns)";
        throw ProfileError(m_file_name, detail.str());
      }
      if(!std::isfinite(value[i])){
        std::ostringstream detail;
        detail << "line " << line_number << " contains a non-finite value";
        throw ProfileError(m_file_name, detail.str());
      }
    }
    G4String trailing;
    if(stream >> trailing){
      std::ostringstream detail;
      detail << "line " << line_number
             << " has extra data after the 16 expected columns";
      throw ProfileError(m_file_name, detail.str());
    }

    if(value[0] == 0.){
      std::ostringstream detail;
      detail << "line " << line_number << " has zero charge";
      throw ProfileError(m_file_name, detail.str());
    }
    if(value[1] <= 0.){
      std::ostringstream detail;
      detail << "line " << line_number
             << " has non-positive momentum " << value[1] << " GeV/c";
      throw ProfileError(m_file_name, detail.str());
    }

    Anchor anchor{};
    anchor.p_gev = value[1];
    std::copy(value.begin()+2, value.begin()+6, anchor.mean.begin());

    anchor.covariance[0][0] = value[6];
    anchor.covariance[0][1] = anchor.covariance[1][0] = value[7];
    anchor.covariance[0][2] = anchor.covariance[2][0] = value[8];
    anchor.covariance[0][3] = anchor.covariance[3][0] = value[9];
    anchor.covariance[1][1] = value[10];
    anchor.covariance[1][2] = anchor.covariance[2][1] = value[11];
    anchor.covariance[1][3] = anchor.covariance[3][1] = value[12];
    anchor.covariance[2][2] = value[13];
    anchor.covariance[2][3] = anchor.covariance[3][2] = value[14];
    anchor.covariance[3][3] = value[15];

    std::ostringstream context;
    context << "line " << line_number << " covariance";
    Cholesky(anchor.covariance, context.str());
    m_anchors[ChargeSign(value[0])].push_back(anchor);
  }

  if(m_anchors.empty())
    throw ProfileError(m_file_name, "profile table has no data rows");

  for(auto& group : m_anchors){
    auto& anchors = group.second;
    std::sort(anchors.begin(), anchors.end(),
              [](const Anchor& lhs, const Anchor& rhs)
              { return lhs.p_gev < rhs.p_gev; });
    for(std::size_t i=1; i<anchors.size(); ++i){
      if(anchors[i-1].p_gev == anchors[i].p_gev){
        std::ostringstream detail;
        detail << "duplicate q=" << (group.first > 0 ? "+" : "-")
               << " momentum anchor at " << anchors[i].p_gev << " GeV/c";
        throw ProfileError(m_file_name, detail.str());
      }
    }
  }

  G4cout << "K18PhaseSpaceProfile: loaded Gaussian-moments ASCII "
         << m_file_name;
  for(const auto& group : m_anchors){
    const auto& anchors = group.second;
    G4cout << "  q=" << (group.first > 0 ? "+" : "-")
           << ":" << anchors.size() << " anchor(s) ["
           << anchors.front().p_gev << ", " << anchors.back().p_gev
           << "] GeV/c";
  }
  G4cout << G4endl;
}

//_____________________________________________________________________________
void
K18PhaseSpaceProfile::LoadRoot()
{
  TFile input(m_file_name.c_str(), "READ");
  if(input.IsZombie())
    throw ProfileError(m_file_name, "cannot open ROOT profile");

  G4String tree_name = m_root_tree;
  if(tree_name.empty()){
    tree_name = input.Get("k18profile") ? "k18profile" : "bcout";
  }
  TTree* tree = dynamic_cast<TTree*>(input.Get(tree_name.c_str()));
  if(!tree)
    throw ProfileError(m_file_name, "ROOT file has no TTree " + tree_name);

  if(tree_name == "bcout"){
    LoadBcOutRoot(tree);
    return;
  }

  RequireScalarLeaf(tree, m_file_name, "signed_charge", "Int_t",
                    tree_name.c_str());
  for(const char* branch : {"p_anchor_gev", "x_mm", "y_mm", "u", "v"})
    RequireScalarLeaf(tree, m_file_name, branch, "Double_t",
                      tree_name.c_str());
  const G4bool has_momentum = tree->GetBranch("p_gev") != nullptr;
  if(has_momentum)
    RequireScalarLeaf(tree, m_file_name, "p_gev", "Double_t",
                      tree_name.c_str());

  G4int signed_charge = 0;
  G4double p_anchor_gev = 0.;
  G4double x_mm = 0.;
  G4double y_mm = 0.;
  G4double u = 0.;
  G4double v = 0.;
  G4double p_gev = 0.;
  tree->SetBranchAddress("signed_charge", &signed_charge);
  tree->SetBranchAddress("p_anchor_gev", &p_anchor_gev);
  tree->SetBranchAddress("x_mm", &x_mm);
  tree->SetBranchAddress("y_mm", &y_mm);
  tree->SetBranchAddress("u", &u);
  tree->SetBranchAddress("v", &v);
  if(has_momentum)
    tree->SetBranchAddress("p_gev", &p_gev);

  const Long64_t entries = tree->GetEntries();
  for(Long64_t entry=0; entry<entries; ++entry){
    if(tree->GetEntry(entry) <= 0){
      std::ostringstream detail;
      detail << "failed to read k18profile entry " << entry;
      throw ProfileError(m_file_name, detail.str());
    }
    if(signed_charge != -1 && signed_charge != 1){
      std::ostringstream detail;
      detail << "entry " << entry
             << " has signed_charge=" << signed_charge
             << ", expected -1 or +1";
      throw ProfileError(m_file_name, detail.str());
    }
    const Vector track{x_mm, y_mm, u, v};
    if(!std::isfinite(p_anchor_gev) || p_anchor_gev <= 0.){
      std::ostringstream detail;
      detail << "entry " << entry << " has invalid p_anchor_gev="
             << p_anchor_gev;
      throw ProfileError(m_file_name, detail.str());
    }
    if(!std::all_of(track.begin(), track.end(),
                   [](G4double value){ return std::isfinite(value); })){
      std::ostringstream detail;
      detail << "entry " << entry << " has a non-finite phase-space value";
      throw ProfileError(m_file_name, detail.str());
    }
    if(has_momentum && (!std::isfinite(p_gev) || p_gev <= 0.)){
      std::ostringstream detail;
      detail << "entry " << entry << " has invalid p_gev=" << p_gev;
      throw ProfileError(m_file_name, detail.str());
    }

    auto& anchors = m_empirical_anchors[signed_charge];
    auto anchor = std::find_if(
      anchors.begin(), anchors.end(),
      [&](const EmpiricalAnchor& item)
      { return SameAnchor(item.p_gev, p_anchor_gev); });
    if(anchor == anchors.end()){
      anchors.push_back(EmpiricalAnchor{p_anchor_gev, {}});
      anchor = anchors.end()-1;
    }
    anchor->tracks.push_back(
      Sample{x_mm, y_mm, u, v, has_momentum ? p_gev : 0., has_momentum});
  }

  if(m_empirical_anchors.empty())
    throw ProfileError(m_file_name, "TTree k18profile has no data rows");

  G4cout << "K18PhaseSpaceProfile: loaded empirical-exact ROOT "
         << m_file_name
         << "  associated momentum=" << (has_momentum ? "p_gev" : "none");
  for(auto& group : m_empirical_anchors){
    auto& anchors = group.second;
    std::sort(anchors.begin(), anchors.end(),
              [](const EmpiricalAnchor& lhs, const EmpiricalAnchor& rhs)
              { return lhs.p_gev < rhs.p_gev; });
    G4cout << "  q=" << (group.first > 0 ? "+" : "-") << ":";
    for(const auto& anchor : anchors)
      G4cout << " p=" << anchor.p_gev << "(" << anchor.tracks.size()
             << " tracks)";
  }
  G4cout << G4endl;
}

//_____________________________________________________________________________
void
K18PhaseSpaceProfile::LoadBcOutRoot(TTree* tree)
{
  if(m_root_charge != -1 && m_root_charge != 1)
    throw ProfileError(m_file_name,
                       "direct bcout ROOT input requires root charge -1 or +1");
  if(!std::isfinite(m_root_p_gev) || m_root_p_gev <= 0.)
    throw ProfileError(m_file_name,
                       "direct bcout ROOT input requires positive anchor momentum");
  if(!std::isfinite(m_root_z_mm))
    throw ProfileError(m_file_name, "bcout target z must be finite");
  if(!std::isfinite(m_root_chisqr_max) || m_root_chisqr_max <= 0.)
    throw ProfileError(m_file_name, "bcout chi-square maximum must be positive");
  if(m_root_max_tracks <= 0)
    throw ProfileError(m_file_name, "bcout maximum track count must be positive");

  RequireScalarLeaf(tree, m_file_name, "ntrack", "Int_t", "bcout");
  for(const char* branch : {"chisqr", "x0", "y0", "u0", "v0"}){
    TLeaf* leaf = tree->GetLeaf(branch);
    if(!leaf)
      throw ProfileError(m_file_name,
                         "tree bcout is missing branch " + G4String(branch));
    if(G4String(leaf->GetTypeName()) != "Double_t")
      throw ProfileError(m_file_name,
                         "bcout branch " + G4String(branch)
                         + " has type " + leaf->GetTypeName()
                         + ", expected Double_t");
  }

  tree->SetBranchStatus("*", 0);
  for(const char* branch : {"ntrack", "chisqr", "x0", "y0", "u0", "v0"})
    tree->SetBranchStatus(branch, 1);
  tree->SetCacheSize(64*1024*1024);
  tree->AddBranchToCache("*", true);

  TTreeReader reader(tree);
  TTreeReaderValue<Int_t> ntrack(reader, "ntrack");
  TTreeReaderArray<Double_t> chisqr(reader, "chisqr");
  TTreeReaderArray<Double_t> x0(reader, "x0");
  TTreeReaderArray<Double_t> y0(reader, "y0");
  TTreeReaderArray<Double_t> u0(reader, "u0");
  TTreeReaderArray<Double_t> v0(reader, "v0");

  auto& tracks = m_empirical_anchors[m_root_charge];
  tracks.push_back(EmpiricalAnchor{m_root_p_gev, {}});
  auto& selected = tracks.back().tracks;
  selected.reserve(static_cast<std::size_t>(m_root_max_tracks));
  std::mt19937_64 rng(0x4b313842634f7574ULL);
  G4long n_selected = 0;
  while(reader.Next()){
    if(*ntrack != 1)
      continue;
    if(chisqr.GetSize() < 1 || x0.GetSize() < 1 || y0.GetSize() < 1
       || u0.GetSize() < 1 || v0.GetSize() < 1)
      throw ProfileError(m_file_name, "short array in bcout track");
    if(!std::isfinite(chisqr[0]) || !(chisqr[0] < m_root_chisqr_max))
      continue;
    const Sample sample{
      x0[0] + u0[0]*m_root_z_mm,
      y0[0] + v0[0]*m_root_z_mm,
      u0[0], v0[0], 0., false};
    if(!std::isfinite(sample.x_mm) || !std::isfinite(sample.y_mm)
       || !std::isfinite(sample.u) || !std::isfinite(sample.v))
      continue;

    ++n_selected;
    if(static_cast<G4long>(selected.size()) < m_root_max_tracks){
      selected.push_back(sample);
    } else {
      std::uniform_int_distribution<G4long> choose(0, n_selected-1);
      const G4long slot = choose(rng);
      if(slot < m_root_max_tracks)
        selected[static_cast<std::size_t>(slot)] = sample;
    }
  }

  if(selected.empty())
    throw ProfileError(m_file_name, "bcout has no selected finite tracks");

  G4cout << "K18PhaseSpaceProfile: loaded direct bcout ROOT "
         << m_file_name << " tree=bcout charge=" << m_root_charge
         << " p_anchor=" << m_root_p_gev << " GeV/c z=" << m_root_z_mm
         << " mm chisqr<" << m_root_chisqr_max
         << " selected=" << n_selected << " stored=" << selected.size()
         << G4endl;
}

//_____________________________________________________________________________
K18PhaseSpaceProfile::Sample
K18PhaseSpaceProfile::Shoot(G4double charge, G4double p_gev,
                            G4long sample_index) const
{
  if(m_mode == Mode::EmpiricalExact)
    return ShootEmpiricalExact(charge, p_gev, sample_index);
  return ShootGaussian(charge, p_gev);
}

//_____________________________________________________________________________
K18PhaseSpaceProfile::Sample
K18PhaseSpaceProfile::ShootGaussian(G4double charge, G4double p_gev) const
{
  const G4int sign = ChargeSign(charge);
  const auto group = m_anchors.find(sign);
  if(group == m_anchors.end()){
    std::ostringstream detail;
    detail << "has no anchors for charge sign " << (sign > 0 ? "+" : "-");
    throw ProfileError(m_file_name, detail.str());
  }
  if(!std::isfinite(p_gev) || p_gev <= 0.){
    std::ostringstream detail;
    detail << "cannot sample at invalid momentum " << p_gev << " GeV/c";
    throw ProfileError(m_file_name, detail.str());
  }

  const auto& anchors = group->second;
  Vector mean{};
  Matrix covariance{};
  if(p_gev <= anchors.front().p_gev){
    mean = anchors.front().mean;
    covariance = anchors.front().covariance;
  } else if(p_gev >= anchors.back().p_gev){
    mean = anchors.back().mean;
    covariance = anchors.back().covariance;
  } else {
    const auto upper = std::upper_bound(
      anchors.begin(), anchors.end(), p_gev,
      [](G4double p, const Anchor& anchor){ return p < anchor.p_gev; });
    const auto& high = *upper;
    const auto& low = *(upper-1);
    const G4double fraction =
      (p_gev - low.p_gev)/(high.p_gev - low.p_gev);
    for(std::size_t i=0; i<mean.size(); ++i){
      mean[i] = low.mean[i] + fraction*(high.mean[i] - low.mean[i]);
      for(std::size_t j=0; j<mean.size(); ++j){
        covariance[i][j] = low.covariance[i][j]
          + fraction*(high.covariance[i][j] - low.covariance[i][j]);
      }
    }
  }

  std::ostringstream context;
  context << "interpolated covariance at q=" << (sign > 0 ? "+" : "-")
          << ", p=" << p_gev << " GeV/c";
  const Matrix lower = Cholesky(covariance, context.str());
  Vector normal{};
  Vector sampled = mean;
  for(auto& value : normal)
    value = G4RandGauss::shoot(0., 1.);
  for(std::size_t i=0; i<sampled.size(); ++i){
    for(std::size_t j=0; j<=i; ++j)
      sampled[i] += lower[i][j]*normal[j];
  }
  return Sample{sampled[0], sampled[1], sampled[2], sampled[3], 0., false};
}

//_____________________________________________________________________________
K18PhaseSpaceProfile::Sample
K18PhaseSpaceProfile::ShootEmpiricalExact(G4double charge,
                                           G4double p_gev,
                                           G4long sample_index) const
{
  const G4int sign = ChargeSign(charge);
  const auto group = m_empirical_anchors.find(sign);
  if(group == m_empirical_anchors.end()){
    std::ostringstream detail;
    detail << "empirical profile has no tracks for charge sign "
           << (sign > 0 ? "+" : "-");
    throw ProfileError(m_file_name, detail.str());
  }
  if(!std::isfinite(p_gev) || p_gev <= 0.){
    std::ostringstream detail;
    detail << "cannot sample empirical profile at invalid momentum "
           << p_gev << " GeV/c";
    throw ProfileError(m_file_name, detail.str());
  }

  const auto& anchors = group->second;
  const auto anchor = std::find_if(
    anchors.begin(), anchors.end(),
    [&](const EmpiricalAnchor& item)
    { return SameAnchor(item.p_gev, p_gev); });
  if(anchor == anchors.end()){
    std::ostringstream detail;
    detail << "empirical exact mode has no q=" << (sign > 0 ? "+" : "-")
           << " anchor at p=" << std::setprecision(12) << p_gev
           << " GeV/c; available anchors:";
    for(const auto& item : anchors)
      detail << " " << item.p_gev;
    throw ProfileError(m_file_name, detail.str());
  }
  if(anchor->tracks.empty())
    throw ProfileError(m_file_name, "empirical anchor has no stored tracks");

  const std::size_t n_tracks = anchor->tracks.size();
  const std::size_t index = sample_index >= 0
    ? static_cast<std::size_t>(sample_index)%n_tracks
    : std::min(n_tracks-1,
               static_cast<std::size_t>(G4UniformRand()*n_tracks));
  return anchor->tracks[index];
}

//_____________________________________________________________________________
G4bool
K18PhaseSpaceProfile::HasCharge(G4double charge) const
{
  if(m_mode == Mode::EmpiricalExact)
    return m_empirical_anchors.find(ChargeSign(charge))
      != m_empirical_anchors.end();
  return m_anchors.find(ChargeSign(charge)) != m_anchors.end();
}

//_____________________________________________________________________________
G4int
K18PhaseSpaceProfile::ChargeSign(G4double charge)
{
  if(!std::isfinite(charge) || charge == 0.)
    throw std::runtime_error(
      "K18PhaseSpaceProfile requires a finite non-zero particle charge");
  return charge > 0. ? 1 : -1;
}

//_____________________________________________________________________________
G4bool
K18PhaseSpaceProfile::SameAnchor(G4double lhs, G4double rhs)
{
  const G4double scale = std::max({1., std::abs(lhs), std::abs(rhs)});
  return std::abs(lhs-rhs) <= 1.e-9*scale;
}

//_____________________________________________________________________________
K18PhaseSpaceProfile::Matrix
K18PhaseSpaceProfile::Cholesky(const Matrix& covariance,
                               const G4String& context)
{
  Matrix lower{};
  for(std::size_t i=0; i<lower.size(); ++i){
    for(std::size_t j=0; j<=i; ++j){
      G4double value = covariance[i][j];
      for(std::size_t k=0; k<j; ++k)
        value -= lower[i][k]*lower[j][k];
      if(i == j){
        if(!std::isfinite(value) || value <= 0.)
          throw std::runtime_error("K18PhaseSpaceProfile " + context
                                   + " is not positive definite");
        lower[i][j] = std::sqrt(value);
      } else {
        lower[i][j] = value/lower[j][j];
        if(!std::isfinite(lower[i][j]))
          throw std::runtime_error("K18PhaseSpaceProfile " + context
                                   + " produced a non-finite Cholesky factor");
      }
    }
  }
  return lower;
}
