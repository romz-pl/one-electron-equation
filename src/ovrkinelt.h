#ifndef RSCHR_OVRKINELT__H
#define RSCHR_OVRKINELT__H


//
// AUTHOR: Zbigniew Romanowski [romz@wp.pl]
//

//
// Returns element of overlap (S) and kinetic (K) matrix between two B-splines .
//
// By definition, we have:
//
//   S_{i,j} = \int B_i(r) B_j(r) dr
//
//   K_{i,j} = \int \nabla B_i(r) \cdot \nabla B_j(r) dr
//


#include <vector>
#include "paramdb.h"


class OvrKinElt
{
public:
    explicit OvrKinElt(const ParamDb& db);

    void OvrKin(const int* idx, const int* jdx, double& ovr, double& kin) const;

private:
    void Precalc(const ParamDb& db);

private:
    // One-dimensional overlap integrals
    std::vector<double> m_ovr;

    // One-dimensional derivatives used to calculated kinetic integrals
    std::vector<double> m_der;

};

//
// Returns element of overlap (ovr) and kinetic (kin) matrix between two B-splines .
//
inline
void OvrKinElt::OvrKin(const int* idx, const int* jdx, double& ovr, double& kin) const
{
const int jx = abs(idx[0] - jdx[0]);
const int jy = abs(idx[1] - jdx[1]);
const int jz = abs(idx[2] - jdx[2]);

    ovr = m_ovr[jx] * m_ovr[jy] * m_ovr[jz];

    const double dx = m_der[jx] * m_ovr[jy] * m_ovr[jz];
    const double dy = m_ovr[jx] * m_der[jy] * m_ovr[jz];
    const double dz = m_ovr[jx] * m_ovr[jy] * m_der[jz];

    kin = 0.5 * (dx + dy + dz);
}


#endif
