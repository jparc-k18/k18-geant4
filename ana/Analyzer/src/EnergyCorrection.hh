#ifndef EnergyCorrection_hh 
#define EnergyCorrection_hh  1

#define PION   1
#define KAON   2
#define PROTON 3

#include "ThreeVector.hh"

int CorrElossOut(double *mom_new, double *energy_new, double mom, int particle, ThreeVector dir, ThreeVector vertex);
double diffE(int particle, double E, double length1, double length2, double Elast);
int caldE(double momentum, double mass, double distance, double *momentum_cor, double *energy_cor, int tag);
double calc_dE_dx(double beta);
double calc_dE_dx2(double beta);
double mygamma(double beta);
double mybeta(double energy,double mormentum);
double calcLengthInTarget(ThreeVector dir, ThreeVector vertex);
double calcLengthInWindow(ThreeVector dir, ThreeVector vertex);

#endif 
