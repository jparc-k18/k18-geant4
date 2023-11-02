#ifndef INCLUDE_MYFUNC
#define INCLUDE_MYFUNC

#include<cmath>
#include<cstdlib>
#include<fstream>
using namespace std;

double GaussRand(double mean,double sigma){
  if(sigma){
    //srand(55);
    //srand((unsigned int)time(NULL));
    double r;
    while(1){
      r = ((double)rand()/RAND_MAX -0.5)*8*sigma + mean;
      if(exp(-(r-mean)*(r-mean)/2/sigma/sigma) >= (double)rand()/RAND_MAX) break;
    }
    return r;
  }
  else{
    return mean;
  }
}

#endif
