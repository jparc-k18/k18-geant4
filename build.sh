#!/bin/sh

set -e
type cmake3 >/dev/null 2>&1 && alias cmake=cmake3

main_dir=$(dirname `readlink -f $0`)
src_dir=$main_dir/src
linkdef_dir=$main_dir/linkdef
obj_dir=$main_dir/obj

if [ "$(uname)" == 'Darwin' -a ! -e "$src_dir/Dict.cc" ]; then
    rootcling -f $src_dir/Dict.cc -c TVector3.h TParticle.h \
	      $linkdef_dir/LinkDef.h
fi

mkdir -pv $obj_dir
cd $obj_dir
#cmake .. -DCMAKE_INSTALL_PREFIX=$G4WORKDIR
cmake .. -DCMAKE_INSTALL_PREFIX=$main_dir
cmake --build . -- -j4
#cmake --install .

if [ "$(uname)" == 'Darwin' -a ! -e "$obj_dir/Dict_rdict.pcm" ]; then
    cp $src_dir/Dict_rdict.pcm $obj_dir
fi
