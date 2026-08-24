# K1.8 upstream transport and missing-mass generator

## Summary

This working tree adds a reusable K1.8 QQDQQ simulation path to G4S2S.  It
contains the upstream beam-line geometry and field, phase-space beam sampling,
target-plane transport, configurable missing-mass kinematics, truth output,
and reproducible run controls.

The implementation has two public generator IDs:

- `Generator: 6380` launches a K1.8 phase-space beam.  It accepts either a
  source-plane state or a target-plane state that is back-propagated to the
  source.
- `Generator: 6381` launches the incident beam upstream and creates the
  configured reaction when the transported track reaches the sampled target
  position.  Target containment is checked using a physical volume or a box.

The QQDQQ field can use the analytic magnet model or a regular three-dimensional
field map.  When an S-2S map is also enabled, `K18CombinedField` exposes both
fields through one Geant4 field manager.

## Configuration overview

### Geometry and field

- `UseK18Beamline` enables the QQDQQ geometry, K1.8 detector branches, and
  upstream field.
- `K18FLDMAP`, `K18FLDNMR`, and `K18FLDCALC` select and scale a regular field
  map.  Without `K18FLDMAP`, the analytic QQDQQ field is used.
- `K18GlobalScale` and `K18Q10Scale` through `K18Q13Scale`/`K18D4Scale` scale
  the analytic magnets or back-propagation field.
- `K18BeamPipe`, `K18VacuumWindows`, `K18BFTRealisticMaterial`, `K18BcOutFilms`,
  and the associated thickness/material keys control passive material.
- `DisableBAC` allows upstream-only transport without constructing BAC.

### Phase-space beam

`K18PhaseSpace*` keys configure generator 6380.  The principal groups are:

- particle and momentum: `Particle`, `P`, `PSigma`, `DPOverP`, `PHalfWidth`,
  `PMin`, `PMax`, and `RejectOutsideRange`;
- transverse state: `XMean/XSigma/XHalfWidth`, the corresponding `Y`, `U`, and
  `V` keys, plus `ZMean/ZSigma/ZHalfWidth` for target-plane sampling;
- source and propagation: `SourcePlane`, `SourceOffset`, `UseTargetPlane`,
  `BackPropStep`, `BackPropMaxPath`, and `BackPropMaxTrial`;
- profile input: `K18PhaseSpaceProfile`, `K18PhaseSpaceProfileP`, and
  `K18PhaseSpaceUseProfileMomentum`.

An ASCII profile stores charge, momentum, four means `(x,y,u,v)`, and the ten
independent covariance elements.  A ROOT profile uses tree `k18profile` with
`signed_charge`, `p_anchor_gev`, `x_mm`, `y_mm`, `u`, and `v`; optional `p_gev`
provides an event-associated momentum.

For direct use of a calibrated `BcOutTracking` ROOT file, set
`K18PhaseSpaceProfile` to that file.  If the ROOT file has no `k18profile`
tree, the loader automatically selects its `bcout` tree; the selection can be
made explicit with `K18PhaseSpaceProfileTree: bcout`.  Direct mode reads
`ntrack`, `chisqr`, `x0`, `y0`, `u0`, and `v0`, applies `ntrack == 1` and
`chisqr < 5` by default, and extrapolates to the target plane as
`x=x0+u0*z`, `y=y0+v0*z`.

One direct `BcOutTracking` file is one charge/momentum anchor.  The charge is
derived from `K18PhaseSpaceParticle`, and the anchor momentum is selected by
`K18PhaseSpaceProfileP` (falling back to `K18PhaseSpaceP`).  The direct-mode
keys are:

```text
K18PhaseSpaceProfile:          /path/to/run74894_BcOutTracking.root
K18PhaseSpaceProfileTree:      bcout       # optional; auto-detected
K18PhaseSpaceProfileTargetZ:   1503.0      # mm
K18PhaseSpaceProfileChiSqMax:  5.0
K18PhaseSpaceProfileMaxTracks: 200000
K18PhaseSpaceProfileP:         1.4         # GeV/c
```

Generator 6381 accepts the corresponding `ReactionBeamProfile*` overrides.
When the selected track count exceeds `*MaxTracks`, a deterministic reservoir
sample is retained, preserving empirical correlations and tails without the
Gaussian-moments `.dat` conversion.  Existing ASCII and `k18profile` ROOT
inputs remain supported.

### Missing-mass reaction

Generator 6381 reads `ReactionBeam*` keys first and falls back to matching
`K18PhaseSpace*` beam settings.  Reaction products, angular acceptance, target
and recoil masses, Fermi motion, bound-state mass, and recoil line shape are
configured by the `Reaction*` keys in `K18MissingMassReaction.cc`.

The output tree records `PrimPi`, `ReactionBeam`,
`ReactionBeamVertexTransport`, `ReactionBeamVertex`, `ReactionTarget`,
`ReactionResidual`, `ReactionScat`, and `ReactionRecoil`.  For generator 6381,
`ReactionBeamVertexTransport` and `ReactionBeamVertex` are filled from the same
live transported state.  With `ReactionRequireBeamTransportMatch: 1`, the
event writer checks their PDG, four-vector, and position before filling the
tree.  The default momentum and position tolerances are `1e-6` in Geant4 units
and can be set with `ReactionBeamTransportMatchMomentumTolerance` and
`ReactionBeamTransportMatchPositionTolerance`.  Events in which the
transported beam does not satisfy target containment are not stored as
reaction entries.

`TargetSizeX`, `TargetSizeY`, and `TargetSizeZ` override the E63 box target
dimensions used by both the detector geometry and reaction containment.
`TargetMaterial: Li6` and `TargetMaterial: Li7` select enriched lithium at
`0.534 g/cm3`. `TargetMaterial: C` uses bulk graphite at `1.8 g/cm3`; convert
an intended areal density `T` to `TargetSizeZ` with
`TargetSizeZ[mm] = 10*T[g/cm2]/1.8`.

Five ready-to-run reaction configurations are provided:

| Configuration | Reaction | Beam momentum |
| --- | --- | ---: |
| `k18_6381_e90_cusp_example.conf` | `d(K-,pi-)Sigma-N` cusp | 1.4 GeV/c |
| `k18_6381_e75_phase1_li7_example.conf` | `7Li(K-,K+)7_XiH` | 1.8 GeV/c |
| `k18_6381_kpi_c12lambda_example.conf` | `12C(K-,pi-)12_LambdaC` | 1.4 GeV/c |
| `k18_6381_kk_c12xibe_example.conf` | `12C(K-,K+)12_XiBe` | 1.8 GeV/c |
| `k18_6381_pik_c12lambda_example.conf` | `12C(pi+,K+)12_LambdaC` | 1.4 GeV/c |

The E90 line-shape name in this tree is `flatte_cusp`; the E90 configuration
uses that name directly.

### Reproducibility

- `K18G4_SEED` fixes the run-level Geant4 seed.
- `K18G4_EVENT_SEED_BASE` derives an independent seed for each global event.
- `K18G4_EVENT_OFFSET` preserves event numbering and random streams across
  split jobs.
- `ReactionKinematicsSeedBase` can isolate reaction sampling from the
  geometry-dependent random consumption of upstream transport.

Cross-geometry regression tests should use the event seed because later events
otherwise inherit the random consumption of all preceding detector transport.
All integer environment controls are parsed strictly and checked for overflow.

## File guide

### New files

| File | Purpose |
| --- | --- |
| `include/K18BeamlineDetectorBuilder.hh` | Public interface for constructing the K1.8 beam line. |
| `src/K18BeamlineDetectorBuilder.cc` | Data-driven QQDQQ iron/coil geometry, ducts, windows, BFT, BC3/BC4, BH1, and BH2 placement. |
| `include/K18BeamlineField.hh` | Coordinate-frame, QQDQQ field, global adapter, and combined-field interfaces. |
| `src/K18BeamlineField.cc` | Analytic magnet field, mapped-field selection, coordinate conversion, and S-2S field composition. |
| `include/K18FieldMap.hh` | Regular-grid Geant4 field-map interface. |
| `src/K18FieldMap.cc` | Strict map loading, grid validation, scaling, and trilinear interpolation. |
| `include/K18MissingMassPrimaryGenerator.hh` | Incident beam generator interface for IDs 6380 and 6381. |
| `src/K18MissingMassPrimaryGenerator.cc` | Beam/profile sampling, source-plane construction, target-state back-propagation, and truth registration. |
| `include/K18MissingMassReaction.hh` | Reaction arming API and Geant4 discrete-process declaration. |
| `src/K18MissingMassReaction.cc` | Live-track vertex approach, target containment, two-body kinematics, recoil models, and truth output. |
| `include/K18PhaseSpaceProfile.hh` | Gaussian-moment and empirical-profile sampler interface. |
| `src/K18PhaseSpaceProfile.cc` | Validated ASCII, empirical `k18profile` ROOT, and direct `BcOutTracking` ROOT loading, plus anchor sampling. |
| `include/K18ReactionModelSampler.hh` | Interfaces for recoil-line-shape and Fermi-momentum sampling. |
| `src/K18ReactionModelSampler.cc` | Rejection samplers for the configured recoil line shape and Gaussian/deuteron Fermi motion. |
| `include/K18RunControl.hh` | Strict environment parsing, global event indexing, and checked seed arithmetic. |

### Updated files

| File | Change |
| --- | --- |
| `include/FdComplex.hh` | Makes header-defined coefficient storage safe across translation units. |
| `include/FermiMotion.hh` | Makes header-defined sampler storage safe across translation units. |
| `include/GeneratorParticleBranches.hh` | Registers K1.8 beam and reaction truth branches and particle IDs. |
| `include/S2SAnaManager.hh` | Adds target-plane truth fields and accessors. |
| `include/S2SDetectorConstruction.hh` | Declares K1.8 beam-line construction. |
| `include/S2SPrimaryGeneratorAction.hh` | Owns and dispatches the K1.8 generator module. |
| `src/MaterialList.cc` | Adds enriched Li6/Li7 target materials and the BC3/BC4 Ar:isobutane:methylal gas mixture. |
| `src/S2SAnaManager.cc` | Creates K1.8 hit/truth branches, records global event IDs, filters target misses, validates the generator-6381 beam-state match, and guards optional histograms. |
| `src/S2SDetectorConstruction.cc` | Builds the upstream geometry, applies configured E63 target dimensions/material, extends the world, and installs accurate long-baseline field integration. |
| `src/S2SPhysicsList.cc` | Attaches the missing-mass reaction process to the configured incident particle. |
| `src/S2SPrimaryGeneratorAction.cc` | Dispatches generators 6380/6381 and applies optional event-level seeds. |
| `src/S2SRunAction.cc` | Provides fixed or high-resolution run seeds. |
| `src/S2SSteppingAction.cc` | Interpolates the transported primary state at the nominal target plane for K1.8 runs. |

## Verification

- The normal working simulation directory and this staging tree both completed
  a full CMake build with Geant4 11.2.2.
- A 100-event missing-mass regression used `K18G4_SEED=314159` and
  `K18G4_EVENT_SEED_BASE=271828`.  Reference and candidate each stored 51
  events with identical event IDs.  All seven truth branches agreed exactly:
  `mismatches=0`, `max_abs_dp4_mev=0`, and `max_abs_dx4=0`.
- A 1.8 GeV/c central-ray field-map test produced one G4 entry and one K1.8
  tracking entry on distinct canonical paths.  Truth values agreed with the
  reference, including `p0=pTgtTruth=1800 MeV/c`; the only coordinate-level
  difference was numerical roundoff below `1e-12 mm`.

This tree is prepared for review only.  No commit or push has been made.
