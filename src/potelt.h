#ifndef RSCHR_POTELT__H
#define RSCHR_POTELT__H

//
// AUTHOR: Zbigniew Romanowski [romz@wp.pl]
//

//
// 1. Calculates potential integrals:
//         U_{i,j} = \int B_i(r) U(r) B_j(r) dr
//
// 2. Integral is evaluated numericaly in cubature points {q}:
//         U_{i,j} \approx \sum_q  B_i(r_q)  U(r_q)  B_j(r_q) dr
//
// 3. In order to seep up calculations, it uses:
//       a) cached values of potential.
//       b) cached value of B-splies.
//




#include "cubat.h"
#include "recmtx.h"
#include "potcache.h"
#include "poten.h"
#include "graphbspl.h"

class PotElt : public GraphBspl
{
public:
    PotElt(Poten& pot, const ParamDb& db);
    ~PotElt();

    // void PrecalcPot(const Poten& pot);
    double Pot(const int* idx, const int* jdx);

    // int PotCacheSize(void) const { return m_potCache->Size(); }

    void FreePotCache(void)
    {
        delete m_potCache;
        m_potCache = NULL;
    }

private:
    void PrecalcBspl(void);
    static int IdxBsplQ(int px, int py, int pz, int w);


private:
    // Potemtial
    Poten& m_pot;

    // Begining of the mesh (cube)
    const double m_x0, m_y0, m_z0;

    // Mesh parameter
    const double m_h;

    // B-spline degree
    const int m_deg;

    // Tabulated (cached) values of "interaction potential" in quadrature nodes
    PotCache *m_potCache;

    // Three-dimensional quadrature over cube [-1,1]^3
    Cubat m_cubat;

    // Tabulated (cached) values of B-spline in quadrature nodes
    RecMtx m_bspl;

    double* m_potbuf;
};

//
// Returns index of pre-calculated B-spline at cubature points
//
inline
int PotElt::IdxBsplQ(int px, int py, int pz, int w)
{
    // return w * w * px + w * py + pz;
    return w * (w * px + py) + pz;
}


#endif
