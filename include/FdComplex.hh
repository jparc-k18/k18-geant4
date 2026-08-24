#ifndef FERMIMOTION_HH
#define FERMIMOTION_HH

#include <TComplex.h>
#include <TMath.h>

namespace FdComplex
{
    const double msp = 1189.37;
    const double ms0 = 1192.642;
    const double mp = 938.2720813;
    const double mn = 939.5654133;
    const double hbarc = 197.327;
    inline double Cj[11] = {0.88472985,
                     -0.26408759,
                     -0.44114404 * 0.1,
                     -0.14397512 * 100.,
                     0.85591256 * 100.,
                     -0.31876761 * 1000.,
                     0.70336701 * 1000.,
                     -0.90049586 * 1000.,
                     0.66145441 * 1000.,
                     -0.25958894 * 1000.,
                     0.};
    const double m0 = 0.9;          //[fm-1]
    const double gamma = 0.2315380; //[fm-1]

    inline TComplex GetFd_Spn(double x, double q)
    {
        double Q = q * mp / (msp + mp);
        Cj[10] = 0.;
        for (int i = 0; i < 10; ++i)
        {
            Cj[10] += Cj[i];
        }
        Cj[10] = -1. * Cj[10];
        double mj[11];
        for (int i = 0; i < 11; ++i)
        {
            mj[i] = (gamma + ((double)i * m0)) * hbarc;
        }
        TComplex Fd_c = TComplex(0., 0.);

        if (x > msp + mn)
        {
            double k = sqrt((x * x - pow(msp + mn, 2)) * (x * x - pow(msp - mn, 2))) / (2. * x);

            TComplex Fd_cj[11];
            for (int i = 0; i < 11; ++i)
            {
                TComplex coeff = TComplex(Cj[i] / (sqrt(4. * acos(-1)) * Q), 0.);
                TComplex atan_c = TComplex((Q / mj[i]) / (1 + (k * k / (mj[i] * mj[i]))),
                                           (Q * k / (mj[i] * mj[i])) / (1 + (k * k / (mj[i] * mj[i]))));
                // std::cout<<"a, "<<i<<" atan_c="<<atan_c.Rho2()<<std::endl;
                // getchar();
                Fd_cj[i] = coeff * TComplex::ATan(atan_c);
                Fd_c += Fd_cj[i];
            }
        }
        else if (x < msp + mn)
        {
            double k = sqrt(-1. * (x * x - pow(msp + mn, 2)) * (x * x - pow(msp - mn, 2))) / (2. * x);
            double k_fm = k * hbarc;
            TComplex Fd_cj[11];
            for (int i = 0; i < 11; ++i)
            {
                TComplex coeff = TComplex(Cj[i] / (sqrt(4. * acos(-1)) * Q), 0.);
                TComplex atan_c = TComplex((Q / mj[i]) / (1 + (k / mj[i])), 0.);
                Fd_cj[i] = coeff * TComplex::ATan(atan_c);
                Fd_c += Fd_cj[i];
            }
        }
        return Fd_c;
    }

    inline TComplex GetFd_S0p(double x, double q)
    {
        double Q = q * mn / (ms0 + mn);
        Cj[10] = 0.;
        for (int i = 0; i < 10; ++i)
        {
            Cj[10] += Cj[i];
        }
        Cj[10] = -1. * Cj[10];
        double mj[11];
        for (int i = 0; i < 11; ++i)
        {
            mj[i] = (gamma + ((double)i * m0)) * hbarc;
        }
        TComplex Fd_c = TComplex(0., 0.);

        if (x > ms0 + mp)
        {
            double k = sqrt((x * x - pow(ms0 + mp, 2)) * (x * x - pow(ms0 - mp, 2))) / (2. * x);

            TComplex Fd_cj[11];
            for (int i = 0; i < 11; ++i)
            {
                TComplex coeff = TComplex(Cj[i] / (sqrt(4. * acos(-1)) * Q), 0.);
                TComplex atan_c = TComplex((Q / mj[i]) / (1 + (k * k / (mj[i] * mj[i]))),
                                           (Q * k / (mj[i] * mj[i])) / (1 + (k * k / (mj[i] * mj[i]))));
                // std::cout<<"a, "<<i<<" atan_c="<<atan_c.Rho2()<<std::endl;
                // getchar();
                Fd_cj[i] = coeff * TComplex::ATan(atan_c);
                Fd_c += Fd_cj[i];
            }
        }
        else if (x < ms0 + mp)
        {
            double k = sqrt(-1. * (x * x - pow(ms0 + mp, 2)) * (x * x - pow(ms0 - mp, 2))) / (2. * x);
            double k_fm = k * hbarc;
            TComplex Fd_cj[11];
            for (int i = 0; i < 11; ++i)
            {
                TComplex coeff = TComplex(Cj[i] / (sqrt(4. * acos(-1)) * Q), 0.);
                TComplex atan_c = TComplex((Q / mj[i]) / (1 + (k / mj[i])), 0.);
                Fd_cj[i] = coeff * TComplex::ATan(atan_c);
                Fd_c += Fd_cj[i];
            }
        }
        return Fd_c;
    }

    inline double f_single(double x, double a, double b, double q)
    {
        TComplex k1_c, k2_c;
        double Rs = 0.; //, k1=0., k2=0.;
        double A = a / hbarc;
        double B = b / hbarc;
        TComplex A_c = TComplex(A, B);
        TComplex i_c = TComplex(0, 1.);

        if (x > msp + mn)
        {
            double k1 = sqrt((x * x - pow(msp + mn, 2)) * (x * x - pow(msp - mn, 2))) / (2. * x);
            if (x > ms0 + mp)
            {
                double k2 = sqrt((x * x - pow(ms0 + mp, 2)) * (x * x - pow(ms0 - mp, 2))) / (2. * x);
                //      Rs = B / (pow(1 + sqrt(2./3.)*k1*B + sqrt(1./3.)*k2*B,2) + pow(sqrt(2./3.)*k1*A + sqrt(1./3.)*k2*A,2));
                k1_c = TComplex(k1, 0.);
                k2_c = TComplex(k2, 0.);
            }
            else if (x < ms0 + mp)
            {
                double k2 = sqrt(-1. * (x * x - pow(ms0 + mp, 2)) * (x * x - pow(ms0 - mp, 2))) / (2. * x);
                k1_c = TComplex(k1, 0.);
                k2_c = TComplex(0, k2);
            }
        }
        else if (x < msp + mn)
        {
            double k1 = sqrt(-1. * (x * x - pow(msp + mn, 2)) * (x * x - pow(msp - mn, 2))) / (2. * x);
            double k2 = sqrt(-1. * (x * x - pow(ms0 + mp, 2)) * (x * x - pow(ms0 - mp, 2))) / (2. * x);
            k1_c = TComplex(0., k1);
            k2_c = TComplex(0., k2);
        }

        if (x == msp + mn || x == ms0 + mp)
            return 0.;
        else
        {
            // TComplex f_Spn = sqrt((2./3.)*B)/(1.-((2./3.)*i_c*k1_c + (1./3.)*i_c*k2_c)*A_c);
            // TComplex f_S0p = sqrt((1./3.)*B)/(1.-((2./3.)*i_c*k1_c + (1./3.)*i_c*k2_c)*A_c);
            // TComplex Fd_Spn = FdComplex::GetFd_Spn(x,q);
            // TComplex Fd_S0p = FdComplex::GetFd_S0p(x,q);
            // TComplex Tele_Sp = TComplex(sqrt(585.),0.);
            // TComplex Tele_S0 = TComplex(sqrt(68.),0.);

            TComplex f_Spn = sqrt(B) / (1. - (i_c * k1_c) * A_c);
            TComplex Fd_Spn = FdComplex::GetFd_Spn(x, q);
            TComplex Tele_Sp = TComplex(sqrt(585.), 0.);
            TComplex f_tot = (Fd_Spn * Tele_Sp * f_Spn);
            //    TComplex f_tot = (Tele_Sp*f_Spn)+(Tele_S0*f_S0p);
            Rs = f_tot.Rho2();
            return Rs * 10000.;
        }
    }

    inline double f(double x, double a, double b, double q)
    {
        TComplex k1_c, k2_c;
        double Rs = 0.; //, k1=0., k2=0.;
        double A = a / hbarc;
        double B = b / hbarc;
        TComplex A_c = TComplex(A, B);
        TComplex i_c = TComplex(0, 1.);

        if (x > msp + mn)
        {
            double k1 = sqrt((x * x - pow(msp + mn, 2)) * (x * x - pow(msp - mn, 2))) / (2. * x);
            if (x > ms0 + mp)
            {
                double k2 = sqrt((x * x - pow(ms0 + mp, 2)) * (x * x - pow(ms0 - mp, 2))) / (2. * x);
                //      Rs = B / (pow(1 + sqrt(2./3.)*k1*B + sqrt(1./3.)*k2*B,2) + pow(sqrt(2./3.)*k1*A + sqrt(1./3.)*k2*A,2));
                k1_c = TComplex(k1, 0.);
                k2_c = TComplex(k2, 0.);
            }
            else if (x < ms0 + mp)
            {
                double k2 = sqrt(-1. * (x * x - pow(ms0 + mp, 2)) * (x * x - pow(ms0 - mp, 2))) / (2. * x);
                k1_c = TComplex(k1, 0.);
                k2_c = TComplex(0, k2);
            }
        }
        else if (x < msp + mn)
        {
            double k1 = sqrt(-1. * (x * x - pow(msp + mn, 2)) * (x * x - pow(msp - mn, 2))) / (2. * x);
            double k2 = sqrt(-1. * (x * x - pow(ms0 + mp, 2)) * (x * x - pow(ms0 - mp, 2))) / (2. * x);
            k1_c = TComplex(0., k1);
            k2_c = TComplex(0., k2);
        }

        if (x == msp + mn || x == ms0 + mp)
            return 0.;
        else
        {
            TComplex f_Spn = sqrt((2. / 3.) * B) / (1. - ((2. / 3.) * i_c * k1_c + (1. / 3.) * i_c * k2_c) * A_c);
            TComplex f_S0p = sqrt((1. / 3.) * B) / (1. - ((2. / 3.) * i_c * k1_c + (1. / 3.) * i_c * k2_c) * A_c);
            TComplex Fd_Spn = FdComplex::GetFd_Spn(x, q);
            TComplex Fd_S0p = FdComplex::GetFd_S0p(x, q);
            TComplex Tele_Sp = TComplex(sqrt(585.), 0.);
            TComplex Tele_S0 = TComplex(sqrt(68.), 0.);

            TComplex f_tot = (Fd_Spn * Tele_Sp * f_Spn) + (Fd_S0p * Tele_S0 * f_S0p);
            //    TComplex f_tot = (Tele_Sp*f_Spn)+(Tele_S0*f_S0p);

            Rs = f_tot.Rho2();

            return Rs * 10000.;
        }
    }

}; // namespace FdComplex
#endif
