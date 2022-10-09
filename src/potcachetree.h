#ifndef RSCHR_POTCACHETREE__H
#define RSCHR_POTCACHETREE__H

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


#include <map>
#include <cassert>


#include "potcache.h"

class PotCacheTree : public PotCache
{
public:
    PotCacheTree(int nxMax, int nyMax, int q);
    virtual ~PotCacheTree(void);

    virtual double* Get(int ix, int iy, int iz);

    virtual double* Insert(int ix, int iy, int iz, double* val);

private:
    int Id(int ix, int iy, int iz) const;

private:
    // Maximal index in "X" and "Y" direction respectively
    const int m_nxMax, m_nyMax;

    // Number of cubature nodes. The length of arrays in the tree.
    const int m_q;

    // Thee of pairs (index, pointer to array)
    std::map<int, double* > m_set;

};

//
// Returns unique indetifier for triple (ix, iy, iz)
//
inline
int PotCacheTree::Id(int ix, int iy, int iz) const
{
    assert(ix >= 0 && ix < m_nxMax);
    assert(iy >= 0 && iy < m_nyMax);

    return ix + (iy + iz * m_nyMax) * m_nxMax;
}

//
// Returns values associated with the triple (ix, iy, iz)
//
inline
double* PotCacheTree::Get(int ix, int iy, int iz)
{
    std::map<int, double*>::iterator it = m_set.find( Id(ix, iy, iz) );

    if( it != m_set.end() )
        return it->second;

    return nullptr;
}


#endif
