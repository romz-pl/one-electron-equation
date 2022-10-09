#include "potelt.h"
#include "bspline.h"
#include "paramdb.h"
#include "util.h"
#include "potcachehash.h"
#include "potcachetree.h"

//
// Constructor
//
PotElt::PotElt( Poten& pot, const ParamDb& db )
    : GraphBspl( db )
    , m_pot      ( pot )
    , m_x0       ( db.GetDbl("Domain_X0") )
    , m_y0       ( db.GetDbl("Domain_Y0") )
    , m_z0       ( db.GetDbl("Domain_Z0") )
    , m_h        ( db.GetDbl("Domain_h" ) )
    , m_deg      ( db.GetInt("General_BsplDeg") )
    , m_potCache ( nullptr )
    , m_cubat    ( db.GetInt("General_CubatOrder") )
{
    // Maximum index in "X, Y"
    const int nx = db.GetInt("Domain_Nx");
    const int ny = db.GetInt("Domain_Ny");
    const int nz = db.GetInt("Domain_Nz");

    if( db.GetStr("General_PotCache") == "Tree" )
    {
        m_potCache = new PotCacheTree(nx, ny, m_cubat.Q());
    }
    else if( db.GetStr("General_PotCache") == "Hash" )
    {
        m_potCache = new PotCacheHash(nx, ny, nz, m_cubat.Q());
    }
    else
    {
        Throw_invalid_argument("ERROR-INPUT-PARAMETER: Wrong value of 'General_PotCache'. Currently supported values: Tree, Hash.");
    }

    m_potbuf = new double [m_cubat.Q()];

    PrecalcBspl();
}

//
// Destructor
//
PotElt::~PotElt()
{
    delete m_potCache;
    delete [] m_potbuf;
}

/*
//
// Precalculates (caches) potential values in qubature points
//
void PotElt::PrecalcPot(const Poten& pot)
{
// Minimum "X, Y, Z" coordinates
const double x0 = Db_GetDbl("Domain_X0");
const double y0 = Db_GetDbl("Domain_Y0");
const double z0 = Db_GetDbl("Domain_Z0");

// Maximum index in "X, Y"
const int nx = Db_GetInt("Domain_Nx");
const int ny = Db_GetInt("Domain_Ny");

const int P = m_deg + 1;
const int Q = m_cubat.Q();

std::vector<double> vecq(Q);
int px, py, pz, ix, iy, iz, q;

    // Create cache with sufficiently large arguments of constructor!
    delete m_potCache;
    m_potCache = new PotCache(nx + 1, ny + 1);

    for(int m = 0; m < NodeNo(); ++m)
    {
        const int *idx = GetIdx(m, 0);

        for(ix = 0; ix < P; ++ix)
        {
            px = idx[0] + ix;
            for(iy = 0; iy < P; ++iy)
            {
                py = idx[1] + iy;
                for(iz = 0; iz < P; ++iz)
                {
                    pz = idx[2] + iz;

                    if(m_potCache->Is(px, py, pz))
                        continue;

                    for(q = 0; q < Q; q++)
                    {
                        const double r = m_cubat.R(q);
                        const double s = m_cubat.S(q);
                        const double t = m_cubat.T(q);
                        const double w = m_cubat.W(q);

                        const double x = x0 + m_h * (0.5 * (r + 1) + px);
                        const double y = y0 + m_h * (0.5 * (s + 1) + py);
                        const double z = z0 + m_h * (0.5 * (t + 1) + pz);

                        vecq[q] = w * pot.Get(x, y, z);
                    }
                    m_potCache->Insert(px, py, pz, vecq);
                }
            }
        }
        // if((m+1) % 100 == 0) { fprintf(stdout, "x"); fflush(stdout); }
    }

}
*/


//
// Returns numerical approximation integral potential
//    U_{i,j} \approx \sum_q  B_i(r_q)  U(r_q)  B_j(r_q) dr
//
// idx - defines basis function B_i
// jdx - defines basis function B_j
//
double PotElt::Pot(const int* idx, const int* jdx)
{
const int P = m_deg + 1;
const int Q = m_cubat.Q();
int beg[3], end[3];
double potTot = 0, cellPot;
double * pot;


    for(int i = 0; i < 3; i++)
    {
        beg[i] = std::max(idx[i], jdx[i]);
        end[i] = std::min(idx[i], jdx[i]) + m_deg;
    }

    // Over all mesh-cell in range
    for(int px = beg[0]; px <= end[0]; px++)
    {
        for(int py = beg[1]; py <= end[1]; py++)
        {
            for(int pz = beg[2]; pz <= end[2]; pz++)
            {
                cellPot = 0;
                const int ki = IdxBsplQ(px - idx[0], py - idx[1], pz - idx[2], P); // Index of B-Spline B_i
                const int kj = IdxBsplQ(px - jdx[0], py - jdx[1], pz - jdx[2], P); // Index of B-Spline B_j
                // const int kp = IdxPotQ(px, py, pz, m_ny, m_nz); // Index of pre-calculated potential

                pot = m_potCache->Get(px, py, pz);

                if(!pot)
                {
                    for(int q = 0; q < Q; q++)
                    {
                        const double r = m_cubat.R(q);
                        const double s = m_cubat.S(q);
                        const double t = m_cubat.T(q);
                        const double w = m_cubat.W(q);

                        const double x = m_x0 + m_h * (0.5 * (r + 1) + px);
                        const double y = m_y0 + m_h * (0.5 * (s + 1) + py);
                        const double z = m_z0 + m_h * (0.5 * (t + 1) + pz);

                        // m_potEval++;
                        m_potbuf[q] = w * m_pot.Get(x, y, z);
                    }
                    pot = m_potCache->Insert(px, py, pz, m_potbuf);
                }


                // p = m_potCache->Get(px, py, pz);
                for(int q = 0; q < Q; q++)
                {
                    cellPot += m_bspl(q, ki) * m_bspl(q, kj) * pot[q];
                }

                potTot += cellPot;
            }
        }
    }
    const double h2 = 0.5 * m_h;
    return potTot * h2 * h2 * h2;
}

/*
//
// Returns numerical approximation integral potential
//    U_{i,j} \approx \sum_q  B_i(r_q)  U(r_q)  B_j(r_q) dr
//
// idx - defines basis function B_i
// jdx - defines basis function B_j
//
double PotElt::PotOld(const int* idx, const int* jdx) const
{
const int P = m_deg + 1;
const int Q = m_cubat.Q();
int beg[3], end[3];
double potVal = 0, cellPot;
int q, px, py, pz;

    for(int i = 0; i < 3; i++)
    {
        beg[i] = std::max(idx[i], jdx[i]);
        end[i] = std::min(idx[i], jdx[i]) + m_deg;
    }

    // Over all mesh-cell in range
    for(px = beg[0]; px <= end[0]; px++)
    {
        for(py = beg[1]; py <= end[1]; py++)
        {
            for(pz = beg[2]; pz <= end[2]; pz++)
            {
                const int ki = IdxBsplQ(px - idx[0], py - idx[1], pz - idx[2], P); // Index of B-Spline B_i
                const int kj = IdxBsplQ(px - jdx[0], py - jdx[1], pz - jdx[2], P); // Index of B-Spline B_j

                assert(m_potCache);
                assert(m_potCache->Is(px, py, pz));
                const std::vector<double>& vec = m_potCache->Get(px, py, pz);

                // Three dimensional Gauss quadrature
                cellPot = 0.0;
                for(q = 0; q < Q; ++q)
                {
                    cellPot += m_bspl(q, ki) * m_bspl(q, kj) * vec[q];
                }


                potVal += cellPot;
            }
        }
    }
    const double h2 = 0.5 * m_h;
    return potVal * h2 * h2 * h2;
}
*/

//
// Precalculates tensor product of B-Splines at cubature nodes
//
void PotElt::PrecalcBspl(void)
{
const int P = m_deg + 1;
const int Q = m_cubat.Q();

    // Printf("B-spline precalculation...");

    m_bspl.Malloc(Q, P * P * P);

    for(int px = 0; px < P; px++)
    {
        for(int py = 0; py < P; py++)
        {
            for(int pz = 0; pz < P; pz++)
            {
                for(int q = 0; q < Q; q++)
                {
                    const double r = m_cubat.R(q);
                    const double s = m_cubat.S(q);
                    const double t = m_cubat.T(q);
                    // const double w = m_cubat.W(q);

                    const double x = 0.5 * (r + 1) + px;
                    const double y = 0.5 * (s + 1) + py;
                    const double z = 0.5 * (t + 1) + pz;

                    const double ax = BSpline::Value(x, m_deg);
                    const double ay = BSpline::Value(y, m_deg);
                    const double az = BSpline::Value(z, m_deg);

                    const int k = IdxBsplQ(px, py, pz, P);
                    m_bspl(q, k) = ax * ay * az;
                }
            }
        }
    }

    // Printf("Ok\n\n");
}

