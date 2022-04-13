#!/bin/sh

nice -n 15 ./bin/StsTracking  ../Param/DCgeom.StsDC1 ../Param/StsMap2500A.dat MomRes2_1.root ../Datfiles/MomRes2.dat &
nice -n 15 ./bin/StsTracking  ../Param/DCgeom.StsDC2 ../Param/StsMap2500A.dat MomRes2_2.root ../Datfiles/MomRes2.dat &
nice -n 15 ./bin/StsTracking  ../Param/DCgeom.StsDC3 ../Param/StsMap2500A.dat MomRes2_3.root ../Datfiles/MomRes2.dat &
nice -n 15 ./bin/StsTracking  ../Param/DCgeom.StsDC4 ../Param/StsMap2500A.dat MomRes2_4.root ../Datfiles/MomRes2.dat &
nice -n 15 ./bin/StsTracking  ../Param/DCgeom.StsDC5 ../Param/StsMap2500A.dat MomRes2_5.root ../Datfiles/MomRes2.dat &
nice -n 15 ./bin/StsTracking  ../Param/DCgeom.StsDC6 ../Param/StsMap2500A.dat MomRes2_6.root ../Datfiles/MomRes2.dat &
nice -n 15 ./bin/StsTracking  ../Param/DCgeom.StsDC7 ../Param/StsMap2500A.dat MomRes2_7.root ../Datfiles/MomRes2.dat &
nice -n 15 ./bin/StsTracking  ../Param/DCgeom.StsDC8 ../Param/StsMap2500A.dat MomRes2_8.root ../Datfiles/MomRes2.dat &

