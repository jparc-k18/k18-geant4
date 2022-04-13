/*
  K18TrackD2UFCN.hh
*/

#ifndef K18TrackD2UFCN_h
#define K18TrackD2UFCN_h 1

#include "Minuit.hh"

class K18Track;
class K18TransMatrix;

class K18TrackD2UFCN : public MinuitFCN
{
  explicit K18TrackD2UFCN( K18Track *track, K18TransMatrix *trMatrix );
  ~K18TrackD2UFCN();

public:
  double operator()( int np, double *g, double *u, int flag );

private:  
  K18Track *Tr_;
  K18TransMatrix *trM_;

  friend class K18Track;
};


#endif
