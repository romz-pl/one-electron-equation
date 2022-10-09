#include "potcachetree.h"

PotCacheTree::PotCacheTree(int nxMax, int nyMax, int q) : m_nxMax(nxMax), m_nyMax(nyMax), m_q(q)
{
}

PotCacheTree::~PotCacheTree(void)
{
std::map<int, double*>::iterator it;

    for(it = m_set.begin(); it != m_set.end(); ++it)
    {
        delete [] it->second;
    }
    m_set.clear();
}


//
// Inserts values associated with the triple (ix, iy, iz) into the cache
//
double* PotCacheTree::Insert(int ix, int iy, int iz, double* val)
{
const int id = Id(ix, iy, iz);
double *p = new double[m_q];

    for(int i = 0; i < m_q; i++)
        p[i] = val[i];

    m_set.insert( std::pair<int, double*>(id, p) );

    return p;
}
