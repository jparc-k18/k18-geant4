// -*- C++ -*-

#ifndef DETECTOR_ID_HH
#define DETECTOR_ID_HH 1

#include <vector>

#include <TString.h>

const Int_t NumOfSegTOF = 18;
const Int_t NumOfSegWC = 12;
const Int_t NumOfLayersSDC1 = 6;
const Int_t NumOfLayersSDC2 = 4;
const Int_t NumOfLayersSDC3 = 4;
const Int_t NumOfLayersSDC4 = 4;
const Int_t NumOfLayersSDC5 = 4;
const Int_t NumOfLayersSSD1 = 4;
const Int_t NumOfLayersSSD2 = 4;
const Int_t NumOfLayersVP = 11;

enum ETriggerFlag
  {
    kSSD1,
    kSSD2,
    kSFTX,
    kSFTUV,
    kVP1,
    kVP2,
    kVP3,
    kVP4,
    kVP5,
    kVP6,
    kVP7,
    kVP8,
    kVP9,
    kVP10,
    kVP11,
    kSDC3,
    kSDC4,
    kSDC5,
    kTOF,
    kAC1,
    kWC,
    kTriggerFlagSize
  };

const std::vector<TString> TriggerFlag =
  {
    "SSD1",
    "SSD2",
    "SFTX",
    "SFTUV",
    "VP1",
    "VP2",
    "VP3",
    "VP4",
    "VP5",
    "VP6",
    "VP7",
    "VP8",
    "VP9",
    "VP10",
    "VP11",
    "SDC3",
    "SDC4",
    "SDC5",
    "TOF",
    "AC1",
    "WC",
  };

#endif
