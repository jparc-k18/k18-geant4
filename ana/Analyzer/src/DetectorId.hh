/*
  DetectorId.hh

  2004/12/16

*/

#ifndef DetectorId_h

#define DetectorId_h 1

const int DetIdSdc =  1;
const int DetIdBdc =  2;
const int DetIdBH1 =  3;
const int DetIdBH2 =  4;
const int DetIdTof =  5;
const int DetIdAc  =  6;
const int DetIdLc  =  7;
const int DetIdBGc =  8;
const int DetIdBAc =  9;
const int DetIdRC  = 10;
const int DetGe    = 11;
const int DetSupr  = 12;
const int DetMisc  = 13;

const int PlMinBcIn  =  1;
const int PlMaxBcIn  = 12;
const int PlMinBcOut = 13;
const int PlMaxBcOut = 26;

const int PlMinSdcIn  =  1;
const int PlMaxSdcIn  = 12;

const int PlMinSdcOut = 41;
const int PlMaxSdcOut = 52;

const int PlOffsBc = 100;

const int PlMinAc =  1;
const int PlMaxAc =  2;
const int PlMinRC =  1;
const int PlMaxRC = 10;
const int PlMinMisc = 1;
const int PlMaxMisc = 4;


const int NumOfLayersBcIn  = PlMaxBcIn  - PlMinBcIn  + 1 + 2;
const int NumOfLayersBcOut = PlMaxBcOut - PlMinBcOut + 1;
const int NumOfLayersSdcIn  = PlMaxSdcIn  - PlMinSdcIn  + 1;
const int NumOfLayersSdcOut = PlMaxSdcOut - PlMinSdcOut + 1;

const int NumOfLayersAc = PlMaxAc - PlMinAc + 1;
const int NumOfLayersRC = PlMaxRC - PlMinRC + 1;
const int NumOfLayersMisc = PlMaxMisc - PlMinMisc + 1;

const int NumOfSegRC = 6;

const int IdK18Target = 130;


const int NumOfSegBH1 =  7;
const int NumOfSegBH2 =  4;
const int NumOfSegTof = 15;
const int NumOfSegAc  = 20;
const int NumOfSegLc  = 14;

const int MaxWireBDC    =  48;
const int MaxWireSDC1X  = 128;
const int MaxWireSDC1Y  =  96;
const int MaxWireSDC1U  =  80;
const int MaxWireSDC2   =  80; 
const int MaxWireSDC3X  =  24; 
const int MaxWireSDC3Y  =  16;
const int MaxWireSDC4X  =  24; 
const int MaxWireSDC4Y  =  24; 


#define E559Analyzer 1

#endif
