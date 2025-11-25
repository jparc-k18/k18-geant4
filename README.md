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
git clone git@github.com:jparc-k18/k18-geant4.git
cd k18geant4
git checkout s2s
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
