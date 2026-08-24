// -*- C++ -*-

#ifndef K18_PHASE_SPACE_PROFILE_HH
#define K18_PHASE_SPACE_PROFILE_HH

#include <array>
#include <map>
#include <vector>

#include <globals.hh>

class TTree;

// Samples transverse beam coordinates from Gaussian moments or an empirical
// ROOT profile, with independent anchors for charge and momentum.
class K18PhaseSpaceProfile
{
public:
  enum class Mode
  {
    GaussianMoments,
    EmpiricalExact
  };

  struct Sample
  {
    G4double x_mm;
    G4double y_mm;
    G4double u;
    G4double v;
    G4double p_gev;
    G4bool has_p;
  };

  explicit K18PhaseSpaceProfile(
    const G4String& file_name,
    const G4String& root_tree = "",
    G4double root_z_mm = 1599.6,
    G4double root_chisqr_max = 5.,
    G4int root_charge = 0,
    G4double root_p_gev = 0.,
    G4long root_max_tracks = 200000);

  Sample Shoot(G4double charge, G4double p_gev,
               G4long sample_index=-1) const;
  G4bool HasCharge(G4double charge) const;
  Mode GetMode() const { return m_mode; }
  G4bool IsEmpirical() const { return m_mode == Mode::EmpiricalExact; }

private:
  using Vector = std::array<G4double, 4>;
  using Matrix = std::array<Vector, 4>;

  struct Anchor
  {
    G4double p_gev;
    Vector mean;
    Matrix covariance;
  };

  struct EmpiricalAnchor
  {
    G4double p_gev;
    std::vector<Sample> tracks;
  };

  using AnchorMap = std::map<G4int, std::vector<Anchor>>;
  using EmpiricalAnchorMap =
    std::map<G4int, std::vector<EmpiricalAnchor>>;

  static G4int ChargeSign(G4double charge);
  static G4bool SameAnchor(G4double lhs, G4double rhs);
  static Matrix Cholesky(const Matrix& covariance,
                         const G4String& context);

  void LoadAscii();
  void LoadRoot();
  void LoadBcOutRoot(TTree* tree);
  Sample ShootGaussian(G4double charge, G4double p_gev) const;
  Sample ShootEmpiricalExact(G4double charge, G4double p_gev,
                             G4long sample_index) const;

  G4String m_file_name;
  G4String m_root_tree;
  G4double m_root_z_mm;
  G4double m_root_chisqr_max;
  G4int m_root_charge;
  G4double m_root_p_gev;
  G4long m_root_max_tracks;
  Mode m_mode;
  AnchorMap m_anchors;
  EmpiricalAnchorMap m_empirical_anchors;
};

#endif
