#!/bin/sh

#nice -n 15 ./bin/Linux-g++/S-2S mac/kaon.mac mac/run160000.mac 450  1100 Conf/k18.conf.default Rootfiles/Acc45Mom1100.root Datfiles/Acc45Mom1100.dat > /dev/null &

command="./bin/analysSksTracking"
conf1="../Param/DCgeom.param_sts300"
conf2="param/StsMap2500A.dat"
rt1="Mom"
rt2="Res"
rt3=".root"
dt1="../Datfiles/Mom"
dt2="Res"
dt3=".dat"
throw="/dev/null"
ComNo=0

for dist in `seq 0 0`
  do
  for p in `seq 0 5`
    do
    t2q1=`expr 600 + ${dist} \* 50`
    mom=`expr 1200 + ${p} \* 50`
    root="${rt1}${mom}${rt2}3${rt3}"
    dat="${dt1}${mom}${dt2}3${dt3}"
    docom="nice -n 15 ${command} ${conf1} ${conf2} ${root} ${dat}"
    ComNo=`ps -ef | grep -v grep | grep "${command}" | wc -l`
    while [ ${ComNo} -gt 6 ]
      do
      sleep 10
      ComNo=`ps -ef | grep -v grep | grep "${command}" | wc -l`
    done
    echo `date` ${docom}
    ${docom} > ${throw} 2>&1 &
  done
done

