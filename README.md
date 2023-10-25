k18geant4
=========

K1.8 geant4 simulation tool.

**Note that this README needs to be modified.**

## Platform

This tool is developed on the platform of KEKCC, CentOS 7.9.2009.
- g++ (GCC) 8.3.0
- ROOT 6.22/08
- Geant4 11.0.2

## How to install

Set environment variables.

```shell
module load git/2260
export G4WORKDIR=$HOME/work/geant4
export MAKEFLAGS="-j40"
. /group/had/sks/software/root/6.22.08/bin/thisroot.sh
. /sw/packages/geant4/11.0.2/bin/geant4.sh
. /sw/packages/geant4/11.0.2/share/Geant4-11.0.2/geant4make/geant4make.sh
```

then

```shell
git clone ssh://sks@www-online.kek.jp:8022/~/public_html/git/k18geant4.git
cd k18geant4
git checkout s2s-dev
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



## Parameters

Some parameter files that are out of the git control should be linked.

```shell
ln -s /group/had/sks/fieldmap/S2S/E70_2023May fieldmap
```
