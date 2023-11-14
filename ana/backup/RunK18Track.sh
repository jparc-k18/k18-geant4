#!/bin/sh

Ana=bin/analysK18Tracking
GeoFile=param/DCgeom.param_0
#MatrixFile=param/K18MatrixN1370FFpFocus.param
#MatrixFile=matrix.D2U.848
MatrixFile=matrix.U2D.848
#MatrixFile=matrix

echo $1 $2
if [ $# -ne 2 ]
    then
    echo Usage: RunK18Track.sh [output root file] [input dat file]
    exit 0
fi

$Ana $GeoFile $MatrixFile $1 $2