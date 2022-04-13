/*
  DCParameters.hh
*/

#ifndef DCParameters_h 
#define DCParameters_h 1 

struct DCPairPlaneInfo 
{
  bool flag;
  int id1, id2;
  double CellSize;
};

extern const DCPairPlaneInfo PPInfoBcOut[], PPInfoSdcIn[];
extern const int NPPInfoBcOut, NPPInfoSdcIn;

#ifdef DefStatic
const DCPairPlaneInfo PPInfoBcIn[] = {
  { false,  1,  0,  0.0 },{ false,  2,  0,  0.0 },{ false,  3,  0,  0.0 },
  { false,  4,  0,  0.0 },{ false,  5,  0,  0.0 },{ false,  6,  0,  0.0 },
  { false,  7,  0,  0.0 },{ false,  8,  0,  0.0 },{ false,  9,  0,  0.0 },
  { false,  10,  0,  0.0 },{ false,  11,  0,  0.0 },{ false,  12,  0,  0.0 }
};


const DCPairPlaneInfo PPInfoBcOut[] = {
  { false,  1,  0,  0.0 },{ false,  2,  0,  0.0 },{ false,  3,  0,  0.0 },
  { false,  4,  0,  0.0 },{ false,  5,  0,  0.0 },{ false,  6,  0,  0.0 },
  { false,  7,  0,  0.0 },{ false,  8,  0,  0.0 },{ false,  9,  0,  0.0 },
  { false,  10,  0,  0.0 },{ false,  11,  0,  0.0 },{ false,  12,  0,  0.0 },
  { false,  13,  0,  0.0 },{ false,  14,  0,  0.0 }
};

/*
const DCPairPlaneInfo PPInfoBcOut[] = {
  { false,  1,  0,  0.0 },{ false,  2,  0,  0.0 },{ false,  3,  0,  0.0 },
  { false,  4,  0,  0.0 },{ false,  5,  0,  0.0 },{ false,  6,  0,  0.0 },
  { false,  7,  0,  0.0 },{ false,  8,  0,  0.0 },{ false,  9,  0,  0.0 },
  { false,  10,  0,  0.0 },{ false,  11,  0,  0.0 },{ false,  12,  0,  0.0 }
};
*/
/*
const DCPairPlaneInfo PPInfoSdcIn[] = {
  { false,  1,  0,  0.0 },{ false,  2,  0,  0.0 },{ false,  3,  0,  0.0 },
  { false,  4,  0,  0.0 },{ false,  5,  0,  0.0 },{ false,  6,  0,  0.0 },
  { false,  7,  0,  0.0 },{ false,  8,  0,  0.0 },{ false,  9,  0,  0.0 },
  { false,  10,  0,  0.0 },{ false,  11,  0,  0.0 },{ false,  12,  0,  0.0 }
};
*/

const DCPairPlaneInfo PPInfoSdcIn[] = {
  { false,  1,  0,  0.0 },{ false,  2,  0,  0.0 },{ false,  3,  0,  0.0 },
  { false,  4,  0,  0.0 },{ false,  5,  0,  0.0 },{ false,  6,  0,  0.0 },
  { false,  7,  0,  0.0 },{ false,  8,  0,  0.0 },{ false,  9,  0,  0.0 },
  { false,  10,  0,  0.0 }
};

const DCPairPlaneInfo PPInfoSdcOut[] = {
  { false,  1,  0,  0.0 },{ false,  2,  0,  0.0 },{ false,  3,  0,  0.0 },
  { false,  4,  0,  0.0 },{ false,  5,  0,  0.0 },{ false,  6,  0,  0.0 },
  { false,  7,  0,  0.0 },{ false,  8,  0,  0.0 },{ false,  9,  0,  0.0 },
  { false,  10,  0,  0.0 },{ false,  11,  0,  0.0 },{ false,  12,  0,  0.0 }
};

const int NPPInfoBcIn = sizeof(PPInfoBcIn)/sizeof(DCPairPlaneInfo);
const int NPPInfoBcOut = sizeof(PPInfoBcOut)/sizeof(DCPairPlaneInfo);
const int NPPInfoSdcIn  = sizeof(PPInfoSdcIn) /sizeof(DCPairPlaneInfo);
const int NPPInfoSdcOut  = sizeof(PPInfoSdcOut) /sizeof(DCPairPlaneInfo);

#endif

const int MinNumOfHitsBcIn   = 8;
const int MinNumOfHitsBcOut  = 10;
const int MinNumOfHitsSdcIn  = 8;
const int MinNumOfHitsSdcOut = 10;

// DL Ranges
const double MinDLBc[13] = {
   0.0,
  -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
  -0.5, -0.5, -0.5, -0.5, -0.5, -0.5
};

const double MaxDLBc[13] = {
  0.0,
  1.8, 1.8, 1.8, 1.8, 1.8, 1.8,
  1.8, 1.8, 1.8, 1.8, 1.8, 1.8
};

const double MinDLSdc[43] = {
  0.0,
  // SDC1
  -0.5, -0.5, -0.5, -0.5,
  // SDC2
  -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
  // Dummy Id=11-30
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  // SDC3
  -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
  // SDC4
  -0.5, -0.5, -0.5, -0.5, -0.5, -0.5
}; 

const double MaxDLSdc[43] = {
  0.0,
  // SDC1
  1.8, 1.8, 1.8, 1.8,
  // SDC2
  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 
  // Dummy Id=11-30
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  // SDC 3
  13.0, 13.0, 13.0, 13.0, 13.0, 13.0,
  // SDC4
  13.0, 13.0, 13.0, 13.0, 13.0, 13.0
};

// extern const int PlIdSdcOutX[], PlIdSdcOutY[];
// extern const int NPlIdSdcOutX, NPlIdSdcOutY;
// const int PlIdSdcOutX[] = { 31, 32, 36, 37, 38, 39 };
// const int PlIdSdcOutY[] = { 33, 34, 42, 43, 44, 45 };
// const int NPlIdSdcOutX = sizeof(PlIdSdcOutX)/sizeof(int);
// const int NPlIdSdcOutY = sizeof(PlIdSdcOutY)/sizeof(int);


#endif
