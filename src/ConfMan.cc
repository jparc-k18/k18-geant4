// -*- C++ -*-

#include "ConfMan.hh"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <TString.h>

#include "BeamMan.hh"
#include "DCGeomMan.hh"
#include "DetSizeMan.hh"
#include "FuncName.hh"
#include "HistMan.hh"

namespace
{
G4bool
ParseDouble(const G4String& raw, G4double& value)
{
  if(raw.empty())
    return false;
  errno = 0;
  char* end = nullptr;
  value = std::strtod(raw.c_str(), &end);
  return errno == 0 && end != raw.c_str() && *end == '\0'
    && std::isfinite(value);
}

G4bool
ParseLong(const G4String& raw, G4long& value)
{
  if(raw.empty())
    return false;
  errno = 0;
  char* end = nullptr;
  value = std::strtol(raw.c_str(), &end, 10);
  return errno == 0 && end != raw.c_str() && *end == '\0';
}

G4bool
ParseBool(const G4String& raw, G4bool& value)
{
  G4String normalized = raw;
  std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                 [](unsigned char ch){ return std::tolower(ch); });
  if(normalized == "1" || normalized == "true"){
    value = true;
    return true;
  }
  if(normalized == "0" || normalized == "false"){
    value = false;
    return true;
  }
  return false;
}

G4String
ParentDirectory(const G4String& file_name)
{
  const auto slash = file_name.find_last_of('/');
  if(slash == G4String::npos)
    return ".";
  if(slash == 0)
    return "/";
  return file_name.substr(0, slash);
}

G4String
JoinChoices(const std::vector<G4String>& choices)
{
  G4String result;
  for(std::size_t i=0; i<choices.size(); ++i){
    if(i != 0)
      result += ",";
    result += choices[i].empty() ? "<empty>" : choices[i];
  }
  return result;
}
}

//_____________________________________________________________________________
ConfMan::ConfMan()
  : m_conf_key("CONF"),
    m_conf_dir(),
    m_conf_buf(),
    m_config_warnings(),
    m_config_errors()
{
}

//_____________________________________________________________________________
ConfMan::~ConfMan()
{
}

//_____________________________________________________________________________
G4bool
ConfMan::Initialize(const G4String& file_name)
{
  m_file.clear();
  m_string.clear();
  m_double.clear();
  m_int.clear();
  m_bool.clear();
  m_config_warnings.clear();
  m_config_errors.clear();

  m_conf_dir = ParentDirectory(file_name);
  m_conf_buf = "\n";
  m_file[m_conf_key] = file_name;
  m_string[m_conf_key] = file_name;

  std::ifstream ifs(file_name);
  if(!ifs.is_open()){
    G4cerr << FUNC_NAME << " file open fail : " << file_name << G4endl;
    return false;
  }

  std::map<G4String, G4int> first_line;
  G4String line;
  G4int line_number = 0;
  while(std::getline(ifs, line)){
    ++line_number;
    const auto comment = line.find('#');
    const G4String content = line.substr(0, comment);
    std::istringstream iss(content);
    G4String key;
    if(!(iss >> key))
      continue;
    if(key.back() == ':')
      key.pop_back();
    if(key.empty())
      continue;

    G4String value;
    if(!(iss >> value))
      value.clear();

    const auto inserted = first_line.emplace(key, line_number);
    if(!inserted.second){
      std::ostringstream message;
      message << "duplicate key " << key << " at lines "
              << inserted.first->second << " and " << line_number;
      m_config_errors.push_back(message.str());
    }

    G4String trailing;
    if(iss >> trailing){
      std::ostringstream message;
      message << "line " << line_number << " key " << key
              << " has extra tokens after value; they were ignored";
      m_config_warnings.push_back(message.str());
    }

    m_conf_buf += line + "\n";
    StoreValue(key, value);
    G4cout << " key = " << std::setw(20) << std::left << key
           << " value = " << std::setw(30) << std::left
           << (value.empty() ? "<missing>" : value) << G4endl;
  }

  const G4bool config_status = ValidateConfiguration();
  PrintConfigurationDiagnostics();

  // These managers are independent during initialization.  Running each one
  // ensures a bad first file does not hide later errors.
  const G4bool parameter_status = InitializeParameterFiles();
  return config_status && parameter_status;
}

//_____________________________________________________________________________
void
ConfMan::StoreValue(const G4String& key, const G4String& value)
{
  m_file[key] = FilePath(value);
  m_string[key] = value;

  G4double double_value = 0.;
  if(!ParseDouble(value, double_value))
    double_value = 0.;
  m_double[key] = double_value;

  G4long long_value = 0;
  if(!ParseLong(value, long_value))
    long_value = 0;
  if(long_value < std::numeric_limits<G4int>::min()
     || long_value > std::numeric_limits<G4int>::max())
    long_value = 0;
  m_int[key] = static_cast<G4int>(long_value);

  G4bool bool_value = false;
  if(!ParseBool(value, bool_value))
    bool_value = long_value != 0;
  m_bool[key] = bool_value;
}

//_____________________________________________________________________________
void
ConfMan::EraseValue(const G4String& key)
{
  m_file.erase(key);
  m_string.erase(key);
  m_double.erase(key);
  m_int.erase(key);
  m_bool.erase(key);
}

//_____________________________________________________________________________
G4bool
ConfMan::ValidateConfiguration()
{
  std::set<G4String> known_keys;
  const auto configured = [this](const G4String& key){
    return m_string.find(key) != m_string.end();
  };
  const auto use_default = [this](const G4String& key,
                                  const G4String& reason,
                                  const char* fallback){
    const G4String invalid = m_string[key].empty()
      ? "<missing>" : m_string[key];
    const G4String replacement = fallback
      ? G4String(fallback) : G4String("<normal fallback chain>");
    m_config_warnings.push_back(
      key + "=" + invalid + " is invalid (" + reason
      + "); ignored and using default " + replacement);
    m_conf_buf += "# ConfMan fallback: " + key + "=" + invalid
      + " -> " + replacement + " (" + reason + ")\n";
    if(fallback)
      StoreValue(key, fallback);
    else
      EraseValue(key);
  };

  const auto validate_bool = [&](const G4String& key, const char* fallback){
    known_keys.insert(key);
    if(!configured(key))
      return;
    G4bool value = false;
    if(!ParseBool(m_string[key], value)){
      use_default(key, "expected boolean 0/1 or true/false", fallback);
      return;
    }
    m_bool[key] = value;
    m_int[key] = value ? 1 : 0;
    m_double[key] = value ? 1. : 0.;
  };

  const auto validate_integer = [&](const G4String& key,
                                    const char* fallback,
                                    const G4String& requirement,
                                    const std::function<G4bool(G4long)>& valid){
    known_keys.insert(key);
    if(!configured(key))
      return;
    G4long value = 0;
    if(!ParseLong(m_string[key], value)
       || value < std::numeric_limits<G4int>::min()
       || value > std::numeric_limits<G4int>::max()
       || !valid(value)){
      use_default(key, requirement, fallback);
      return;
    }
    m_int[key] = static_cast<G4int>(value);
    m_double[key] = static_cast<G4double>(value);
    m_bool[key] = value != 0;
  };

  const auto validate_long = [&](const G4String& key,
                                 const char* fallback,
                                 const G4String& requirement,
                                 const std::function<G4bool(G4long)>& valid){
    known_keys.insert(key);
    if(!configured(key))
      return;
    G4long value = 0;
    if(!ParseLong(m_string[key], value) || !valid(value)){
      use_default(key, requirement, fallback);
      return;
    }
    m_double[key] = static_cast<G4double>(value);
    m_bool[key] = value != 0;
    if(value >= std::numeric_limits<G4int>::min()
       && value <= std::numeric_limits<G4int>::max())
      m_int[key] = static_cast<G4int>(value);
  };

  const auto validate_double = [&](const G4String& key,
                                   const char* fallback,
                                   const G4String& requirement,
                                   const std::function<G4bool(G4double)>& valid){
    known_keys.insert(key);
    if(!configured(key))
      return;
    G4double value = 0.;
    if(!ParseDouble(m_string[key], value) || !valid(value)){
      use_default(key, requirement, fallback);
      return;
    }
    m_double[key] = value;
  };

  const auto validate_choice = [&](const G4String& key,
                                   const char* fallback,
                                   const std::vector<G4String>& choices){
    known_keys.insert(key);
    if(!configured(key))
      return;
    if(std::find(choices.begin(), choices.end(), m_string[key])
       == choices.end())
      use_default(key, "expected one of " + JoinChoices(choices), fallback);
  };

  const auto finite = [](G4double){ return true; };
  const auto positive = [](G4double value){ return value > 0.; };
  const auto nonnegative = [](G4double value){ return value >= 0.; };
  const auto nonzero = [](G4double value){ return value != 0.; };
  const auto positive_integer = [](G4long value){ return value > 0; };

  // Retired switches are structural errors, not invalid-value fallbacks.  An
  // old config must not silently resurrect a superseded geometry path.  Keep
  // checking the rest of the file so all diagnostics are still reported in
  // one startup pass.
  known_keys.insert("K18BeamPipe");
  if(configured("K18BeamPipe")){
    m_config_errors.push_back(
      "retired key K18BeamPipe is not allowed; the adopted continuous "
      "star/rectangular K1.8 beam pipe is always used");
    EraseValue("K18BeamPipe");
  }

  // Core values with an unambiguous application default.
  const std::set<G4long> experiments{10, 63, 70, 90};
  validate_integer("Experiment", "70", "unknown experiment number",
                   [&experiments](G4long value){
                     return experiments.find(value) != experiments.end();
                   });
  const std::set<G4long> generators{
    0, 1, 2, 3, 4, 5, 6, 7, 6301, 6302, 6303, 6374, 6375,
    6380, 6381, 7001, 7002, 7003, 7004, 7501,
    9000, 9001, 9002, 9003, 9004};
  validate_integer("Generator", "1", "unknown generator number",
                   [&generators](G4long value){
                     return generators.find(value) != generators.end();
                   });
  for(const char* key : {"CheckOverlap", "TargetVP", "StopAtIron",
                          "UseK18Beamline", "DisableBAC"})
    validate_bool(key, "0");
  validate_bool("TREE", "1");
  validate_bool("EM", "1");
  validate_bool("DECAY", "0");
  validate_bool("HADRON", "0");
  for(const char* key : {"Q1SCALE", "Q2SCALE", "D1SCALE"})
    validate_double(key, "1.0", "expected a finite number", finite);
  // Legacy generators use the PK18 sign to select the beam charge.
  validate_double("PK18", "1.4", "expected a finite non-zero momentum",
                  nonzero);
  validate_choice("TargetMaterial", "Air",
                  {"Be", "natLi", "Li6", "Li7", "CH2", "HeGas",
                   "AFT", "Air", "LD2", "C"});
  for(const char* key : {"TargetSizeX", "TargetSizeY", "TargetSizeZ"})
    validate_double(key, nullptr, "expected a positive finite size", positive);
  validate_double("TargetGapHalfZ", "0", "expected a non-negative size",
                  nonnegative);

  // K1.8 field and passive-material controls.
  for(const char* key : {"K18GlobalScale", "K18Q10Scale", "K18Q11Scale",
                          "K18D4Scale", "K18Q12Scale", "K18Q13Scale",
                          "K18PhaseSpaceBackPropGlobalScale",
                          "K18PhaseSpaceBackPropQ10Scale",
                          "K18PhaseSpaceBackPropQ11Scale",
                          "K18PhaseSpaceBackPropD4Scale",
                          "K18PhaseSpaceBackPropQ12Scale",
                          "K18PhaseSpaceBackPropQ13Scale",
                          "K18PhaseSpaceBackPropK18Q10Scale",
                          "K18PhaseSpaceBackPropK18Q11Scale",
                          "K18PhaseSpaceBackPropK18D4Scale",
                          "K18PhaseSpaceBackPropK18Q12Scale",
                          "K18PhaseSpaceBackPropK18Q13Scale",
                          "ReactionBeamBackPropGlobalScale",
                          "ReactionBeamBackPropQ10Scale",
                          "ReactionBeamBackPropQ11Scale",
                          "ReactionBeamBackPropD4Scale",
                          "ReactionBeamBackPropQ12Scale",
                          "ReactionBeamBackPropQ13Scale",
                          "ReactionBeamBackPropK18Q10Scale",
                          "ReactionBeamBackPropK18Q11Scale",
                          "ReactionBeamBackPropK18D4Scale",
                          "ReactionBeamBackPropK18Q12Scale",
                          "ReactionBeamBackPropK18Q13Scale"})
    validate_double(key, nullptr, "expected a finite field scale", finite);
  validate_double("K18FLDNMR", "1.0", "expected a finite number", finite);
  validate_double("K18FLDCALC", "1.0", "expected a finite non-zero number",
                  nonzero);
  validate_double("K18TargetL", "1503.0", "expected a positive length",
                  positive);
  validate_double("K18VacuumPressurePa", "0", "expected a non-negative pressure",
                  nonnegative);
  for(const char* key : {"K18BC4LongitudinalShift",
                          "K18BH2LongitudinalShift"})
    validate_double(key, "0", "expected a finite length", finite);
  for(const char* key : {"K18AllVacuum", "K18VacuumWindows"})
    validate_bool(key, "0");
  validate_bool("K18BFTRealisticMaterial", "1");
  validate_bool("K18BFTBlackSheet", "1");
  validate_bool("K18BcOutFilms", "1");
  validate_bool("K18BH2Enabled", "1");
  validate_double("K18VacuumWindowThickness", "0.5",
                  "expected a positive thickness", positive);
  validate_double("K18VacuumWindowRadius", "100.0",
                  "expected a positive radius", positive);
  validate_double("K18BFTFibreDiameter", "1.0",
                  "expected a diameter in (0,2] mm",
                  [](G4double value){ return value > 0. && value <= 2.; });
  validate_double("K18BFTBlackSheetThickness", "0.15",
                  "expected a non-negative thickness", nonnegative);
  validate_double("K18BcOutWindowMylarThickness", "0.0125",
                  "expected a positive thickness", positive);
  validate_double("K18BcOutCathodeMylarThickness", "0.012",
                  "expected a positive thickness", positive);
  validate_double("K18BcOutAnodeCathodeGap", "2.0",
                  "expected a positive gap", positive);
  validate_double("K18BcOutAlCoatingThickness", "0.000005",
                  "expected a non-negative thickness", nonnegative);
  validate_double("K18BcOutCarbonCoatingThickness", "0.000005",
                  "expected a non-negative thickness", nonnegative);
  known_keys.insert("K18VacuumWindowMaterial");
  known_keys.insert("K18BH2Material");

  // Generator 6380/6381 beam controls. nullptr erases an invalid override so
  // the existing ReactionBeam -> K18PhaseSpace -> hard-coded fallback chain
  // remains authoritative.
  const std::vector<G4String> charged_particles{
    "e-", "e+", "mu-", "mu+", "pi-", "pi+", "kaon-", "kaon+",
    "proton", "anti_proton", "deuteron", "triton", "He3", "alpha",
    "sigma-", "sigma+", "xi-", "anti_xi+", "omega-", "anti_omega+"};
  std::vector<G4String> standard_particles = charged_particles;
  standard_particles.insert(standard_particles.end(),
                            {"gamma", "pi0", "kaon0L", "kaon0S",
                             "neutron", "anti_neutron", "lambda",
                             "anti_lambda", "sigma0", "anti_sigma0",
                             "xi0", "anti_xi0"});
  validate_choice("K18PhaseSpaceParticle", "kaon-", charged_particles);
  validate_choice("ReactionBeamParticle", nullptr, charged_particles);
  validate_choice("ReactionScatteredParticle", "kaon+", standard_particles);

  for(const G4String prefix : {G4String("K18PhaseSpace"),
                               G4String("ReactionBeam")}){
    validate_double(prefix + "P", nullptr, "expected a positive momentum",
                    positive);
    validate_double(prefix + "PSigma", nullptr,
                    "expected a non-negative momentum width", nonnegative);
    validate_double(prefix + "DPOverP", nullptr,
                    "expected a non-negative relative width", nonnegative);
    validate_double(prefix + "PHalfWidth", nullptr,
                    "expected a non-negative half-width", nonnegative);
    validate_double(prefix + "PMin", nullptr,
                    "expected a non-negative lower bound", nonnegative);
    validate_double(prefix + "PMax", nullptr,
                    "expected a non-negative upper bound", nonnegative);
    for(const char* axis : {"X", "Y", "U", "V"}){
      validate_double(prefix + axis + "Mean", nullptr,
                      "expected a finite mean", finite);
      validate_double(prefix + axis + "Sigma", nullptr,
                      "expected a non-negative width", nonnegative);
      validate_double(prefix + axis + "HalfWidth", nullptr,
                      "expected a non-negative half-width", nonnegative);
    }
    validate_bool(prefix + "RejectOutsideRange", nullptr);
    validate_bool(prefix + "UseProfileMomentum", nullptr);
    validate_choice(prefix + "SourcePlane", nullptr,
                    {"", "VI", "Q10Upstream", "Q10Front", "Q10Entrance"});
    validate_double(prefix + "SourceOffset", nullptr,
                    "expected a finite source offset", finite);
    validate_double(prefix + "BackPropStep", nullptr,
                    "expected a positive integration step", positive);
    validate_double(prefix + "BackPropMaxPath", nullptr,
                    "expected a positive maximum path", positive);
    validate_integer(prefix + "BackPropMaxTrial", nullptr,
                     "expected a positive integer", positive_integer);
    validate_double(prefix + "ProfileP", nullptr,
                    "expected a positive profile momentum", positive);
    validate_double(prefix + "ProfileTargetZ", nullptr,
                    "expected a finite target z", finite);
    validate_double(prefix + "ProfileChiSqMax", nullptr,
                    "expected a positive chi-square limit", positive);
    validate_long(prefix + "ProfileMaxTracks", nullptr,
                  "expected a positive integer", positive_integer);
  }
  validate_bool("K18PhaseSpaceUseTargetPlane", "1");
  validate_bool("ReactionUseTargetPlane", "1");
  for(const char* key : {"K18PhaseSpaceZMean", "ReactionVertexZMean"})
    validate_double(key, nullptr, "expected a finite mean", finite);
  for(const char* key : {"K18PhaseSpaceZSigma",
                          "K18PhaseSpaceZHalfWidth",
                          "ReactionVertexZSigma",
                          "ReactionVertexZHalfWidth"})
    validate_double(key, nullptr, "expected a non-negative width", nonnegative);

  // Missing-mass reaction controls.
  validate_long("ReactionKinematicsSeedBase", nullptr,
                "expected a positive integer", positive_integer);
  validate_double("ReactionVertexPositionTolerance", "0.001",
                  "expected a positive tolerance", positive);
  validate_double("ReactionVertexApproachFraction", "0.9",
                  "expected a number in (0,1)",
                  [](G4double value){ return value > 0. && value < 1.; });
  validate_double("ReactionTargetContainmentTolerance", "0",
                  "expected a non-negative tolerance", nonnegative);
  validate_integer("ReactionMaxVertexApproachSteps", "16",
                   "expected a positive integer", positive_integer);
  validate_bool("ReactionRequireBeamTransportMatch", "0");
  validate_double("ReactionBeamTransportMatchMomentumTolerance", "0.000001",
                  "expected a positive tolerance", positive);
  validate_double("ReactionBeamTransportMatchPositionTolerance", "0.000001",
                  "expected a positive tolerance", positive);
  validate_choice("ReactionTargetContainmentMode", nullptr,
                  {"box", "physical_volume"});
  validate_choice("ReactionFermiModel", "none",
                  {"none", "None", "NONE", "0", "false", "False",
                   "FALSE", "gaussian", "Gaussian", "deuteron", "Deuteron"});
  validate_choice("ReactionRecoilMassModel", "fixed",
                  {"fixed", "flatte_cusp"});
  for(const char* key : {"ReactionUseFermiMotion", "ReactionUseBoundState",
                          "ReactionFlatCosTheta", "ReactionSampleCM",
                          "ReactionHighMomentumSolution"})
    validate_bool(key, nullptr);
  for(const char* key : {"ReactionRecoilMass", "ReactionTargetMass",
                          "ReactionTargetNucleusMass", "ReactionResidualMass",
                          "ReactionBoundStateMass", "ReactionFermiSigmaP",
                          "ReactionTargetBindingEnergy", "ReactionSeparationEnergy"})
    validate_double(key, nullptr, "expected a non-negative mass or momentum",
                    nonnegative);
  for(const char* key : {"ReactionBoundDMMeV",
                          "ReactionBoundExcessEnergyMeV",
                          "ReactionCuspScatteringLengthReal",
                          "ReactionCuspScatteringLengthImag"})
    validate_double(key, nullptr, "expected a finite number", finite);
  validate_double("ReactionCuspMomentumTransferMeV", "186.0",
                  "expected a positive momentum transfer", positive);
  validate_integer("ReactionTargetPDG", "0", "expected an integer PDG code",
                   [](G4long){ return true; });
  for(const char* key : {"ReactionThetaMin", "ReactionThetaMax"})
    validate_double(key, nullptr, "expected an angle in [0,180] degrees",
                    [](G4double value){ return value >= 0. && value <= 180.; });
  for(const char* key : {"ReactionPhiMin", "ReactionPhiMax"})
    validate_double(key, nullptr, "expected an angle in [0,360] degrees",
                    [](G4double value){ return value >= 0. && value <= 360.; });
  validate_double("ReactionSolvePMax", "5.0",
                  "expected a positive momentum", positive);
  validate_integer("ReactionSolveBins", "1000",
                   "expected an integer of at least 100",
                   [](G4long value){ return value >= 100; });
  validate_integer("ReactionMaxTrial", "1000",
                   "expected a positive integer", positive_integer);
  validate_double("ReactionRecoilMassMin", nullptr,
                  "expected a positive mass", positive);
  validate_double("ReactionRecoilMassMax", nullptr,
                  "expected a positive mass", positive);

  for(const char* key : {"K18FLDMAP", "K18ROOT",
                          "K18PhaseSpaceProfile", "K18PhaseSpaceProfileTree",
                          "ReactionBeamProfile", "ReactionBeamProfileTree",
                          "ReactionTargetParticle", "ReactionRecoilParticle",
                          "ReactionTargetPhysicalVolume"})
    known_keys.insert(key);

  const auto raw_value = [this](const G4String& key){
    const auto found = m_string.find(key);
    return found == m_string.end() ? G4String() : found->second;
  };
  const auto double_value = [this](const G4String& key, G4double fallback){
    const auto found = m_string.find(key);
    return found == m_string.end() || found->second.empty()
      ? fallback : m_double[key];
  };
  const auto bool_value = [this](const G4String& key, G4bool fallback){
    const auto found = m_string.find(key);
    return found == m_string.end() || found->second.empty()
      ? fallback : m_bool[key];
  };

  const auto is_standard_particle = [&standard_particles](const G4String& name){
    return std::find(standard_particles.begin(), standard_particles.end(), name)
      != standard_particles.end();
  };
  if(configured("ReactionTargetParticle")
     && !is_standard_particle(raw_value("ReactionTargetParticle"))
     && double_value("ReactionTargetMass", 0.) <= 0.)
    use_default("ReactionTargetParticle",
                "unknown target particle without a positive mass override",
                "proton");
  if(configured("ReactionRecoilParticle")
     && !is_standard_particle(raw_value("ReactionRecoilParticle"))
     && double_value("ReactionRecoilMass", 0.) <= 0.)
    use_default("ReactionRecoilParticle",
                "unknown recoil particle without a positive mass override",
                "xi-");
  const G4double configured_target_mass =
    double_value("ReactionTargetMass", 0.);
  const G4String binding_key = configured("ReactionTargetBindingEnergy")
    ? G4String("ReactionTargetBindingEnergy") : G4String("ReactionSeparationEnergy");
  if(configured_target_mass > 0. && configured(binding_key)
     && double_value(binding_key, 0.) >= configured_target_mass)
    use_default(binding_key,
                "binding energy must be below the configured target mass",
                "0");

  for(const G4String prefix : {G4String("K18PhaseSpace"),
                               G4String("ReactionBeam")}){
    const G4double minimum = double_value(prefix + "PMin", 0.);
    const G4double maximum = double_value(prefix + "PMax", 0.);
    if(minimum > 0. && maximum > 0. && maximum <= minimum)
      use_default(prefix + "PMax", "upper bound must exceed lower bound",
                  nullptr);
  }

  const G4double theta_min = double_value("ReactionThetaMin", 0.);
  const G4double theta_max = double_value("ReactionThetaMax", 16.);
  if(theta_max < theta_min){
    if(configured("ReactionThetaMin"))
      use_default("ReactionThetaMin", "theta_min must not exceed theta_max",
                  "0");
    if(configured("ReactionThetaMax"))
      use_default("ReactionThetaMax", "theta_max must not be below theta_min",
                  "16");
  }

  if(raw_value("ReactionRecoilMassModel") == "flatte_cusp"){
    const G4double mass_min = double_value("ReactionRecoilMassMin", 0.);
    const G4double mass_max = double_value("ReactionRecoilMassMax", 0.);
    if(mass_min > 0. && mass_max > 0. && mass_max <= mass_min){
      if(configured("ReactionRecoilMassMin"))
        use_default("ReactionRecoilMassMin",
                    "mass_min must be below mass_max", nullptr);
      if(configured("ReactionRecoilMassMax"))
        use_default("ReactionRecoilMassMax",
                    "mass_max must exceed mass_min", nullptr);
    }
    if(bool_value("ReactionUseBoundState",
                  double_value("ReactionBoundStateMass", 0.) > 0.))
      use_default("ReactionRecoilMassModel",
                  "flatte_cusp cannot be combined with bound-state mode",
                  "fixed");
  }

  if(bool_value("ReactionUseBoundState",
                double_value("ReactionBoundStateMass", 0.) > 0.)){
    const G4double nucleus_mass =
      double_value("ReactionTargetNucleusMass", 0.);
    const G4double bound_mass = double_value("ReactionBoundStateMass", 0.);
    const G4double residual_mass = double_value("ReactionResidualMass", 0.);
    if(nucleus_mass <= 0. || (bound_mass <= 0. && residual_mass <= 0.))
      use_default("ReactionUseBoundState",
                  "bound-state mode needs target-nucleus mass and either "
                  "bound-state or residual mass", "0");
  }

  const G4String fermi_model = raw_value("ReactionFermiModel");
  const G4bool fermi_model_enabled =
    !fermi_model.empty() && fermi_model != "none" && fermi_model != "None"
    && fermi_model != "NONE" && fermi_model != "0"
    && fermi_model != "false" && fermi_model != "False"
    && fermi_model != "FALSE";
  if(bool_value("ReactionUseFermiMotion", fermi_model_enabled)
     && !fermi_model_enabled)
    use_default("ReactionUseFermiMotion",
                "Fermi motion was enabled without a sampling model", "0");

  if(raw_value("ReactionTargetContainmentMode") == "physical_volume"
     && raw_value("ReactionTargetPhysicalVolume").empty())
    use_default("ReactionTargetContainmentMode",
                "physical_volume mode needs ReactionTargetPhysicalVolume",
                nullptr);

  for(const G4String prefix : {G4String("K18PhaseSpace"),
                               G4String("ReactionBeam")}){
    if(!bool_value(prefix + "UseProfileMomentum", false))
      continue;
    G4String profile = raw_value(prefix + "Profile");
    if(profile.empty() && prefix == "ReactionBeam")
      profile = raw_value("K18PhaseSpaceProfile");
    if(profile.empty())
      use_default(prefix + "UseProfileMomentum",
                  "profile momentum was requested without a profile", "0");
  }

  // Missing optional K1.8 resources have established alternatives: analytic
  // field for K18FLDMAP and Gaussian sampling for beam profiles.
  for(const char* key : {"K18FLDMAP", "K18PhaseSpaceProfile",
                          "ReactionBeamProfile"}){
    if(!configured(key) || m_string[key].empty()
       || m_string[key] == "none" || m_string[key] == "skip"
       || m_string[key] == "0")
      continue;
    const G4String resolved = FilePath(m_string[key]);
    std::ifstream resource(resolved);
    if(!resource.good()){
      G4String reason = "configured optional resource cannot be opened";
      if(G4String(key) != "K18FLDMAP")
        reason += "; Gaussian fallback is debug-only and invalid for "
          "data-anchored production";
      use_default(key, reason, nullptr);
    } else if(resolved != m_string[key]) {
      StoreValue(key, resolved);
    }
  }

  // Resource validation above can remove a profile that was present during
  // the first cross-field pass.  Re-evaluate this dependency before startup.
  for(const G4String prefix : {G4String("K18PhaseSpace"),
                               G4String("ReactionBeam")}){
    if(!bool_value(prefix + "UseProfileMomentum", false))
      continue;
    G4String profile = raw_value(prefix + "Profile");
    if(profile.empty() && prefix == "ReactionBeam")
      profile = raw_value("K18PhaseSpaceProfile");
    if(profile.empty())
      use_default(prefix + "UseProfileMomentum",
                  "profile momentum was requested without an open profile",
                  "0");
  }

  std::vector<G4String> unsupported_keys;
  for(const auto& item : m_string){
    const G4String& key = item.first;
    const G4bool managed_namespace = key.compare(0, 3, "K18") == 0
      || key.compare(0, 8, "Reaction") == 0;
    if(managed_namespace && known_keys.find(key) == known_keys.end())
      unsupported_keys.push_back(key);
  }
  for(const auto& key : unsupported_keys)
    use_default(key, "unknown or unsupported configuration key", nullptr);

  return m_config_errors.empty();
}

//_____________________________________________________________________________
void
ConfMan::PrintConfigurationDiagnostics() const
{
  if(!m_config_warnings.empty()){
    G4cerr << "#W [ConfMan::Initialize] configuration warnings ("
           << m_config_warnings.size() << "):" << G4endl;
    for(std::size_t i=0; i<m_config_warnings.size(); ++i)
      G4cerr << "#W   [" << i + 1 << "/" << m_config_warnings.size()
             << "] " << m_config_warnings[i] << G4endl;
  }
  if(!m_config_errors.empty()){
    G4cerr << "#E [ConfMan::Initialize] structural configuration errors ("
           << m_config_errors.size() << "):" << G4endl;
    for(std::size_t i=0; i<m_config_errors.size(); ++i)
      G4cerr << "#E   [" << i + 1 << "/" << m_config_errors.size()
             << "] " << m_config_errors[i] << G4endl;
  }
}

//_____________________________________________________________________________
G4bool
ConfMan::InitializeParameterFiles()
{
  G4bool status = true;
  if(!InitializeParameter<DCGeomMan>("DCGEO"))
    status = false;
  if(!InitializeParameter<BeamMan>("BEAM"))
    status = false;
  if(!InitializeParameter<DetSizeMan>("DSIZE"))
    status = false;
  if(!InitializeParameter<HistMan>("HIST"))
    status = false;
  // InitializeParameter<JamMan>("JAM");
  // InitializeParameter<IncMan>("INC");
  if(!status)
    G4cerr << "#E [ConfMan::InitializeParameterFiles] checked all parameter "
           << "managers; one or more failed" << G4endl;
  return status;
}

//_____________________________________________________________________________
G4String
ConfMan::FilePath(const G4String& src) const
{
  if(src.empty())
    return src;
  std::ifstream tmp(src);
  if (tmp.good())
    return src;
  else
    return m_conf_dir + "/" + src;
}
