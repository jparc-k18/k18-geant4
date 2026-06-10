// -*- C++ -*-

#ifndef GENERATORPARTICLEBRANCHES_HH
#define GENERATORPARTICLEBRANCHES_HH

#include <vector>

#include <G4String.hh>
#include <globals.hh>

namespace GeneratorParticleBranches
{
inline constexpr char kPrimPi[]   = "PrimPi";
inline constexpr char kLambda[]   = "Lambda";
inline constexpr char kPrimProt[] = "PrimProt";
inline constexpr char kPrimNeut[] = "PrimNeut";
inline constexpr char kDecPi[]    = "DecPi";
inline constexpr char kDecProt[]  = "DecProt";
inline constexpr char kSigma0[]   = "Sigma0";
inline constexpr char kSigmaP[]   = "SigmaP";
inline constexpr char kHypNucleus[]  = "HypNucleus";
inline constexpr char kHypFragment[] = "HypFragment";
inline constexpr char kRecoilIon[]   = "RecoilIon";
inline constexpr char kDecGamma[]    = "DecGamma";

enum class ParticleId : G4int
{
  None     = -1,
  PrimPi   = 0,
  Lambda   = 1,
  PrimProt = 2,
  PrimNeut = 3,
  DecPi    = 4,
  DecProt  = 5,
  Sigma0   = 6,
  SigmaP   = 7,
  HypNucleus  = 8,
  HypFragment = 9,
  RecoilIon   = 10,
  DecGamma    = 11
};

inline std::vector<G4String> BranchList(G4int generator)
{
  switch(generator){
  case 6376:
    return {kPrimPi};
  case 6375:
    return {kPrimPi, kHypNucleus, kHypFragment, kRecoilIon, kDecPi, kDecGamma};
  case 9001:
    return {kPrimPi, kLambda, kPrimProt, kDecPi, kDecProt};
  case 9002:
    return {kPrimPi, kLambda, kPrimProt, kDecPi, kDecProt};
  case 9003:
    return {kPrimPi, kSigma0, kLambda, kPrimProt, kDecPi, kDecProt};
  case 9004:
    return {kPrimPi, kSigmaP, kPrimNeut, kDecPi, kDecProt};
  case 6370:
  case 6371:
  case 6374:
    return {kDecGamma};
  default:
    return {};
  }
}

} // namespace GeneratorParticleBranches

#endif
