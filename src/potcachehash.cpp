#include "potcachehash.h"

//
// Constructor
//
PotCacheHash::PotCacheHash(int nx, int ny, int nz, int q)
    : m_nx(nx)
    , m_ny(ny)
    , m_nz(nz)
    , m_q(q)
{
    m_value = new double* [nx * ny * nz];

    for(int i = 0; i < nx * ny * nz; i++)
        m_value[i] = nullptr;
}

//
// Destructor
//
PotCacheHash::~PotCacheHash()
{
    for(int i = 0; i < m_nx * m_ny * m_nz; i++)
    {
        if(m_value[i])
            delete [] m_value[i];
    }
    delete [] m_value;
}

//
//
//
double* PotCacheHash::Insert(int ix, int iy, int iz, double* val)
{
const int id = Id(ix, iy, iz);
double *p;

    assert(m_value[id] == nullptr); // It must be the empty slot

    m_value[id] = new double [m_q];

    p = m_value[id];
    for(int i = 0; i < m_q; i++)
        p[i] = val[i];

    return p;
}


