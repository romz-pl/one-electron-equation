#ifndef RSCHR_GRAPHBSPL__H
#define RSCHR_GRAPHBSPL__H


//
// AUTHOR: Zbigniew Romanowski [romz@wp.pl]
//

//
// 1. Generates graph of basis functions, which are tensor product of B-splines.
//
// 2. We use the fact, that doamin is a brick.
//
// 3. Initialy generated graph is partitioned by ParMETIS library.
//
// 4. After partitioning, the new graph is generated with properly re-number nodes.
//
// 5. For re-numbering of nodes, sorting is applied.
//
// 6. Renumbering procedure is quite complicated, since it requires communication
//    with other proceses within MPI library.
//
// 7. After re-numbering new DOFs are consequitive, as required by HYPRE library.
//


#include "graph.h"
#include <mpi.h>
#include "paramdb.h"


class GraphBspl : public Graph
{
public:
    explicit GraphBspl(const ParamDb& db);

    void Generate(MPI_Comm comm);

    int Ilower(void) const { return m_ilower; }
    int Iupper(void) const { return m_iupper; }

private:
//	int NxAux(void) const { return (m_nx - m_deg); }
//	int NyAux(void) const { return (m_ny - m_deg); }
//	int NzAux(void) const { return (m_nz - m_deg); }


private:
    int Dof(int ix, int iy, int iz) const;
    int GenNgb(const int* idx, std::vector<int>& ngb, const int* map) const;
    void Adj2Csr(int*& xadj, int*& adjncy);

    void Split(MPI_Comm comm, int *vtxdist);
    void GenPre(MPI_Comm comm, int minM, int maxM);
    int* Partition(MPI_Comm comm, int *vtxdist);
    void GenPost(MPI_Comm comm, const int *part, const int *map);

    static int* DistFrag(int* arrlocal, int nlocal, MPI_Comm comm);
    static int* GenMap(int* part, int nlocal, int minM, int& mapSize, MPI_Comm comm);

private:
    // Lower and upper limit of DOF
    int m_ilower, m_iupper;

    // Auxiliary domain N_x, N_y, N_z
    const int m_nxAux, m_nyAux, m_nzAux;

    // B-spline degree
    const int m_deg;
};

//
// Returns DOF for triple (ix,iy,iz).
// DOF-s strats from zero.
// DOF-s are consequtive integer numbers.
//
inline
int GraphBspl::Dof(int ix, int iy, int iz) const
{
    return ix + (iy + iz * m_nyAux) * m_nxAux;
}


#endif
