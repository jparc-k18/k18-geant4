#!/bin/sh

nice -n 15 ./bin/StsTracking  ../Param/DCgeom.StsDC1 ../Param/StsMap2500A.dat MomRes1_1.root ../Datfiles/MomRes1.dat &
nice -n 15 ./bin/StsTracking  ../Param/DCgeom.StsDC2 ../Param/StsMap2500A.dat MomRes1_2.root ../Datfiles/MomRes1.dat &
nice -n 15 ./bin/StsTracking  ../Param/DCgeom.StsDC3 ../Param/StsMap2500A.dat MomRes1_3.root ../Datfiles/MomRes1.dat &
#nice -n 15 ./bin/StsTracking  ../Param/DCgeom.StsDC4 ../Param/StsMap2500A.dat MomRes1_4.root ../Datfiles/MomRes1.dat &
# nice -n 15 ./bin/StsTracking  ../Param/DCgeom.StsDC5 ../Param/StsMap2500A.dat MomRes1_5.root ../Datfiles/MomRes1.dat &
# nice -n 15 ./bin/StsTracking  ../Param/DCgeom.StsDC6 ../Param/StsMap2500A.dat MomRes1_6.root ../Datfiles/MomRes1.dat &
# nice -n 15 ./bin/StsTracking  ../Param/DCgeom.StsDC7 ../Param/StsMap2500A.dat MomRes1_7.root ../Datfiles/MomRes1.dat &
# nice -n 15 ./bin/StsTracking  ../Param/DCgeom.StsDC8 ../Param/StsMap2500A.dat MomRes1_8.root ../Datfiles/MomRes1.dat &
nice -n 15 ./bin/StsTracking  ../Param/DCgeom.StsDC11 ../Param/StsMap2500A.dat MomRes1_11.root ../Datfiles/MomRes1.dat &
nice -n 15 ./bin/StsTracking  ../Param/DCgeom.StsDC12 ../Param/StsMap2500A.dat MomRes1_12.root ../Datfiles/MomRes1.dat &
nice -n 15 ./bin/StsTracking  ../Param/DCgeom.StsDC13 ../Param/StsMap2500A.dat MomRes1_13.root ../Datfiles/MomRes1.dat &

