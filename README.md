k18geant4
=========

K1.8 geant4 simulation tool.

## Platform

This tool is developed on the platform of KEKCC, Red Hat Enterprise Linux release 9.3 (Plow).
- g++ (GCC) 11.4.1
- ROOT 6.32.04
- Geant4 11.2.2

## Anaconda setting

To use Python,
it is necessary to build the Anaconda local environment once using the `conda` command as follows.
Note that it is recommended to use `conda install` instead of `pip install` in the anaconda environment.

```sh
$ conda create -n myenv python=3.9 # myenv is an example name
$ conda activate myenv
$ conda install numpy psutil pyyaml rich
```

Add the following line in .bashrc to activate your environment.

```sh
conda activate myenv
```

If the prompt header of conda is annoying, add the following line in .condarc.

```yaml
changeps1: False
```

## How to install

Set environment variables.

```shell
. /sw/packages/root/6.32.04/bin/thisroot.sh
. /sw/packages/geant4/11.2.2/bin/geant4.sh
. /sw/packages/geant4/11.2.2/share/Geant4/geant4make/geant4make.sh
conda activate myenv
export G4WORKDIR=$HOME/work/geant4 # set as you like
export PATH=/group/had/sks/software/unpacker/e70/bin:$PATH
export PATH=$G4WORKDIR/bin/Linux-g++:$PATH
```

then

```shell
git clone git@github.com:jparc-k18/k18-geant4.git k18geant4
cd k18geant4
git checkout k18dev
./build.sh
```



## How to use

Arguments of ConfFile and OutputName are necessary.
G4Macro is an optional argument.

```shell
./bin/G4S2S [ConfFile] [OutputName] (G4Macro)
./bin/G4S2S param/conf/default.conf foo.root
./bin/G4S2S param/conf/default.conf foo.root bar.mac
```

### K18 upstream beam-through

Use `Generator: 6380` for a K18 beam-through sample.  The generator samples
the target-plane phase space, back-propagates it through the K18 field to find
an upstream starting state, and then transports the particle forward once in
Geant4.  `Generator: 6381` adds the configured missing-mass reaction when the
transported beam reaches the sampled target vertex.

The phase-space source may be an idealized Gaussian configuration or a measured
`BcOutTracking` ROOT file.  A calibrated `BcOutTracking` file can be used
directly; no `.dat` conversion is required:

```text
Generator:                    6380
UseK18Beamline:               1
K18PhaseSpaceParticle:        kaon+
K18PhaseSpaceP:               1.4
K18PhaseSpaceProfile:         /path/to/run74894_BcOutTracking.root
K18PhaseSpaceProfileTree:     bcout
K18PhaseSpaceProfileP:        1.4
K18PhaseSpaceProfileTargetZ:  1503.0
K18PhaseSpaceProfileChiSqMax: 5.0
K18PhaseSpaceProfileMaxTracks: 200000
```

For direct `bcout` input, the loader requires `ntrack == 1` and
`chisqr < K18PhaseSpaceProfileChiSqMax`, then evaluates
`x=x0+u0*z` and `y=y0+v0*z` at the configured target-plane `z`.  One raw
`BcOutTracking` ROOT file supplies one charge/momentum anchor; use a separate
configuration or file for each momentum point.  The particle charge selects
the charge sign and `K18PhaseSpaceProfileP` labels the anchor momentum.

The existing ASCII moments format and the empirical ROOT format with a
`k18profile` tree remain supported.  See `README_UPSTREAM.md` for the full
configuration reference and the `6381` reaction overrides.

Minimal configurations are provided in `param/conf/`:

- `k18_6380_gaussian_example.conf`: generator 6380 without an external profile;
- `k18_6380_bcout_example.conf`: generator 6380 with direct `bcout` input;
- `k18_6381_bcout_reaction_example.conf`: generator 6381 with direct `bcout`
  input and an illustrative two-body reaction;
- `k18_6381_e90_cusp_example.conf`: E90
  `d(K-,pi-)Sigma-N` cusp at 1.4 GeV/c;
- `k18_6381_e75_phase1_li7_example.conf`: E75 phase-1
  `7Li(K-,K+)7_XiH` at 1.8 GeV/c;
- `k18_6381_kpi_c12lambda_example.conf`: `12C(K-,pi-)12_LambdaC`;
- `k18_6381_kk_c12xibe_example.conf`: `12C(K-,K+)12_XiBe`;
- `k18_6381_pik_c12lambda_example.conf`: `12C(pi+,K+)12_LambdaC`.

For the two direct-input examples, replace the `/path/to/...BcOutTracking.root`
entry with one calibrated ROOT file.  No converted `.dat` file is needed.
The five reaction examples contain only repository-relative Geant4 paths and
can be launched directly from the repository root.  Their matching K18 and
S-2S analyzer configurations are under
`k18-analyzer/runmanager/runlist/k18_6381_*_{k18,s2s}.conf`.
