# K18 Geant4 / analyzer: quick start

Clone `k18-geant4` and `k18-analyzer` (branch `k18dev`) side by side, e.g. under
`~/pro`. The analyzer configs read `../k18-geant4/param/DCGEO/...`, so keep the
directory name `k18-geant4`. Run G4 in `k18-geant4` and reconstruction in
`k18-analyzer`. Detailed commands and the five reaction pairs:
[README_K18_GEANT4.md](../k18-analyzer/runmanager/runlist/README_K18_GEANT4.md).

Not in git: create `k18-geant4/fieldmap` (link to `/group/had/sks/fieldmap/S2S`)
and the analyzer `param` and `param_release` links before running.

1. BT: `param/conf/k18_6380_gaussian_example.conf` or
   `k18_6380_bcout_example.conf` in Geant4. For direct bcout, supply a calibrated
   ROOT profile; its momentum is configured separately.
2. Tracking: analyzer `runmanager/runlist/k18_6380_{k18,s2s}_example.conf`.
   Validate K18 with `bin/DstK18TrackingGeant4 --check-config CONFIG`.
3. Reaction: use the matching Geant4 and analyzer tag from the five-pair table.
4. Missing mass: `bin/DstMissingMassGeant4 --preset=TAG --g4=G4 --k18=K18
   --s2s=S2S --root=MM --csv=CSV`.

For registered reactions, the convenience pipeline runs all four stages:

```bash
cd ~/pro/k18-analyzer
OUT=/group/had/sks/Users/$USER/k18_runs/kpi_test   # a new path every run
python3 scripts/run_k18_pipeline.py --list-reactions
python3 scripts/run_k18_pipeline.py --reaction kpi_c12lambda -n 1000 \
  -o "$OUT" --g4-dir "$HOME/pro/k18-geant4"
```

The pipeline uses shipped configs; use manual stages for custom configs.
Inspect startup warnings (`#W`), entry alignment, the missing-mass CSV cutflow,
histogram underflow/overflow and fit quality before interpreting results.

K18 configuration has one current format: version 2. BC resolution is DCGEO
`Res`; BFT cluster sigma is separate. The S-2S config needs only the field map
and scales, `S2sGeant4AcceptedPdg` (signed PDG of the tracked particle; its
charge sets the RK polarity) and optionally `DstS2sTrackingSeedMomentum`.
Hit smearing with DCGEO `Res` is applied in the analyzer, not in Geant4. New G4
ROOT uses `G4DCReadoutFrame=chamber-local-v1`. The detector-material toggles are
obsolete; windows, fibres, sheets and chamber films are fixed nominal geometry.

For the full current K18 input contract, see
[config/README.k18_tracking.md](../k18-analyzer/config/README.k18_tracking.md).
