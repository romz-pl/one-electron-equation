#ifndef RSCHR_POTCACHEHASH__H
#define RSCHR_POTCACHEHASH__H

//
// AUTHOR: Zbigniew Romanowski [romz@wp.pl]
//

//
// 1. Cache for potential values
//
// 2. Cache is represented as a map, where triple (ix,iy,iz) is a key
//    and array of doubles are values.
//
// 3. Triple (ix,iy,iz) defines point in space.
//
// 4. Array of doubles represent values of potential in cubature points.
//
// 5. Values of potential stored in the cache are used for evaluation
//    of potential integrals:
//         U_{i,j} = \int B_i(r) U(r) B_j(r) dr
//
// 6. For proper working of the cache, two integers (nxMax, nyMax) must be provided in constructor.
//    These integers determine mapping from (ix,iy,iz) into key used in chache.
//    The maaping will be correct if:
//            nxMax > max{ix} and
//            nyMax > max{iy}
//



#include <cassert>
#include "potcache.h"

class PotCacheHash : public PotCache
{
public:
    PotCacheHash(int nx, int ny, int nz, int q);
    virtual ~PotCacheHash();

    virtual double* Get(int ix, int iy, int iz);

    virtual double* Insert(int ix, int iy, int iz, double* val);

private:

    int Id(int ix, int iy, int iz) const;

private:
    // Maximal index in "X" and "Y" and "Z" direction respectively
    const int m_nx, m_ny, m_nz;

    // Cubature nodes
    const int m_q;

    // Array of pointers to values of potential
    double** m_value;
};

//
// Returns unique indetifier for triple (ix, iy, iz)
//
inline
int PotCacheHash::Id(int ix, int iy, int iz) const
{
    assert(ix >= 0 && ix < m_nx);
    assert(iy >= 0 && iy < m_ny);
    assert(iz >= 0 && iy < m_nz);

    return (ix + (iy + iz * m_ny) * m_nx);
}


//
// Returns pointer to values associated with the triple (ix, iy, iz)
//
inline
double* PotCacheHash::Get(int ix, int iy, int iz)
{
    return m_value[ Id(ix, iy, iz) ];
}

#endif
