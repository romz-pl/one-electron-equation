#ifndef RSCHR_POTCACHE__H
#define RSCHR_POTCACHE__H

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


class PotCache
{
public:
    PotCache(void) { }
    virtual ~PotCache(void) { }

    virtual double* Get(int ix, int iy, int iz) = 0;

    virtual double* Insert(int ix, int iy, int iz, double* val) = 0;
};


#endif
