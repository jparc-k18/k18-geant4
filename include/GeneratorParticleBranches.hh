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
inline constexpr char kReactionBeam[] = "ReactionBeam";
inline constexpr char kReactionBeamVertexTransport[] =
  "ReactionBeamVertexTransport";
inline constexpr char kReactionBeamVertex[] = "ReactionBeamVertex";
inline constexpr char kReactionTarget[] = "ReactionTarget";
inline constexpr char kReactionResidual[] = "ReactionResidual";
inline constexpr char kReactionScat[] = "ReactionScat";
inline constexpr char kReactionRecoil[] = "ReactionRecoil";

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
  DecGamma    = 11,
  ReactionBeam = 12,
  ReactionBeamVertex = 13,
  ReactionTarget = 14,
  ReactionResidual = 15,
  ReactionScat = 16,
  ReactionRecoil = 17,
  ReactionBeamVertexTransport = 18
};

inline std::vector<G4String> BranchList(G4int generator)
{
  switch(generator){
  case 6376:
  case 6380:
    return {kPrimPi};
  case 6381:
    return {kPrimPi, kReactionBeam, kReactionBeamVertexTransport,
            kReactionBeamVertex,
            kReactionTarget, kReactionResidual,
            kReactionScat, kReactionRecoil};
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
