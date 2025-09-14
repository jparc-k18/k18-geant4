#ifndef _FERMIMOTION_H
#define _FERMIMOTION_H
// #ifndef FERMIMOTION_HH
// #define FERMIMOTION_HH

#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"

namespace FermiMotion
{
    double pmag;
    double cos_theta;
    double sin_theta;
    double phi;
    TVector3 p;
    TRandom* R;

    inline double func(double x)
    {
        // Deuteron wave function in p-space (Phys Rept. 149, 1 (1987) Table11)
        double p = x;
        double r = p / 0.197;

        const double sqrtpi2 = 0.79788456;
        const double alpha = 0.231609;

        double m[11], m2[11];
        for (int i = 0; i < 11; i++)
        {
            m[i] = alpha + i * 0.9;
            m2[i] = m[i] * m[i];
        }

        double c[11] = {0.90457337e0, -0.35058661e0, -0.17635927e0, -0.10418261e2,
                        0.45089439e2, -0.14861947e3, 0.31779642e3, -0.37496518e3,
                        0.22560032e3, -0.54858290e2, 0.};

        double d[11] = {0.24133026e-1, -0.64430531e0, 0.51093352e0, -0.54419065e1,
                        0.15872034e2, -0.14742981e2, 0.44956539e1, -0.71152863e-1,
                        0., 0., 0.};

        c[10] = 0.;
        for (int i = 0; i < 10; i++) c[10] -= c[i]; // normalize c[] properly

        int n = 10, n1 = 9, n2 = 8;

        double sum1 = 0.;
        double sum2 = 0.;
        double sum3 = 0.;
        double rtemp;
        for (int i = 0; i < 8; i++)
        {
            rtemp = d[i] / m2[i];
            sum1 = sum1 + rtemp;
            rtemp = d[i];
            sum2 = sum2 + rtemp;
            rtemp = d[i] * m2[i];
            sum3 = sum3 + rtemp;
        }

        for (int i = 0; i < 3; i++)
        { // normalize d[] properly
            d[n2] = -m2[n1] * m2[n] * sum1 + (m2[n1] + m2[n]) * sum2 - sum3;
            d[n2] = d[n2] * m2[n2] / (m2[n] - m2[n2]) / (m2[n1] - m2[n2]);
            int cycle = n2;
            n2 = n1;
            n1 = n;
            n = cycle;
        }

        double U = 0.;
        double W = 0.;
        for (int i = 0; i < 11; i++)
        {
            U += c[i] / (r * r + m2[i]);
            W += d[i] / (r * r + m2[i]);
        }
        U = sqrtpi2 * U; // s wave
        W = sqrtpi2 * W; // d wave

        return r * r * (U * U + W * W); // total probability to have momentum p
    }

    // FermiMotion();
    // void Generate();
    // TVector3 GetMomentum(){return p ;}
    // double GetpMag(){return pmag ;}
    // double GetpX(){return p.x() ;}
    // double GetpY(){return p.y() ;}
    // double GetpZ(){return p.z() ;}
    // double Getfunc(double x){return func(x);}

    inline TVector3 GetMomentum()
    {
        p.SetX(0);
        p.SetY(0);
        p.SetZ(0);
        R = new TRandom(0);

        double tmp = 0;
        while (1)
        {
            //    pmag=gRandom->Uniform(0.,750.);
            //    pmag=gRandom->Uniform(0.,350.);
            //    tmp=gRandom->Uniform(0.,1.);
            pmag = gRandom->Uniform(0, 0.50);
            tmp = gRandom->Uniform(0, 2.5);
            if (tmp < func(pmag)) break;
        }
        cos_theta = R->Uniform(-1., 1.);
        sin_theta = sqrt(1. - cos_theta * cos_theta);
        phi = R->Uniform(0., 2. * acos(-1.));

        // set Fermimom = 0.;
        // pmag=0.;
        // pmag =pmag/30.;
        p.SetX(pmag * sin_theta * TMath::Cos(phi));
        p.SetY(pmag * sin_theta * TMath::Sin(phi));
        p.SetZ(pmag * cos_theta);
        return p;
    }

}; // namespace FermiMotion
#endif