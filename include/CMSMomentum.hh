#ifndef CMSMOMENTUM_HH
#define CMSMOMENTUM_HH

#include <G4Types.hh>

namespace CMSMomentum
{
    inline G4double GetCMSMomentum(G4double M, G4double m1, G4double m2)
    {
        if (M < m1 + m2)
            return 0.0;
        G4double M2 = M * M;
        G4double m12 = m1 * m1;
        G4double m22 = m2 * m2;
        G4double p2 = (M2 * M2 - 2 * M2 * (m12 + m22) + (m12 - m22) * (m12 - m22)) / (4. * M2);
        return std::sqrt(p2);
    }
}

#endif