#!/bin/sh

bsub ./bin/s2sTracking param/DCgeom.param_mod_001 Matrix/D2Umatrix_s2s_46cm_orbit data/D2UShort46cm00_001.root ../Short46cm00.dat
bsub ./bin/s2sTracking param/DCgeom.param_mod_001 Matrix/D2Umatrix_s2s_46cm_orbit data/D2UShort46cm01_001.root ../Short46cm01.dat
bsub ./bin/s2sTracking param/DCgeom.param_mod_001 Matrix/D2Umatrix_s2s_46cm_orbit data/D2UShort46cm02_001.root ../Short46cm02.dat
bsub ./bin/s2sTracking param/DCgeom.param_mod_001 Matrix/D2Umatrix_s2s_46cm_orbit data/D2UShort46cm03_001.root ../Short46cm03.dat

bsub ./bin/s2sTracking param/DCgeom.param_mod_200 Matrix/D2Umatrix_s2s_46cm_orbit data/D2UShort46cm00_200.root ../Short46cm00.dat
bsub ./bin/s2sTracking param/DCgeom.param_mod_200 Matrix/D2Umatrix_s2s_46cm_orbit data/D2UShort46cm01_200.root ../Short46cm01.dat
bsub ./bin/s2sTracking param/DCgeom.param_mod_200 Matrix/D2Umatrix_s2s_46cm_orbit data/D2UShort46cm02_200.root ../Short46cm02.dat
bsub ./bin/s2sTracking param/DCgeom.param_mod_200 Matrix/D2Umatrix_s2s_46cm_orbit data/D2UShort46cm03_200.root ../Short46cm03.dat

bsub ./bin/s2sTracking param/DCgeom.param_mod_020 Matrix/D2Umatrix_s2s_46cm_orbit data/D2UShort46cm00_020.root ../Short46cm00.dat
bsub ./bin/s2sTracking param/DCgeom.param_mod_020 Matrix/D2Umatrix_s2s_46cm_orbit data/D2UShort46cm01_020.root ../Short46cm01.dat
bsub ./bin/s2sTracking param/DCgeom.param_mod_020 Matrix/D2Umatrix_s2s_46cm_orbit data/D2UShort46cm02_020.root ../Short46cm02.dat
bsub ./bin/s2sTracking param/DCgeom.param_mod_020 Matrix/D2Umatrix_s2s_46cm_orbit data/D2UShort46cm03_020.root ../Short46cm03.dat
