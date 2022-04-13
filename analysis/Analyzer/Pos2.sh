#!/bin/sh

#nice -n 15 ./bin/Linux-g++/S-2S mac/kaon.mac mac/run160000.mac 450  1100 Conf/k18.conf.default Rootfiles/Acc45Mom1100.root Datfiles/Acc45Mom1100.dat > /dev/null &

command="./bin/StsTracking"
conf1="param/DCgeom.StsDC"
conf2="param/StsMap2500A.dat"
rt1="Rootfile/MomRes"
#rt2="Res"
rt3=".root"
dt1="../Datfiles/MomRes"
#dt2="Res"
dt3=".dat"
throw="/dev/null"
ComNo=0

for plane in `seq 1 6`
  do
  for res in `seq 1 3`
    do
    conf1_="${conf1}${plane}${res}"
    root="${rt1}22_${plane}${res}${rt3}"
    dat="${dt1}22${dt3}"
    docom="nice -n 15 ${command} ${conf1_} ${conf2} ${root} ${dat}"
    ComNo=`ps -ef | grep -v grep | grep "${command}" | wc -l`
    while [ ${ComNo} -gt 2 ]
      do
      sleep 10
      ComNo=`ps -ef | grep -v grep | grep "${command}" | wc -l`
    done
    echo `date` ${docom}
    ${docom} &#> ${throw} 2>&1 &
  done
done

