#include <cstdio>
#include <cassert>
#include <parmetis.h>
#include "graphbspl.h"
#include "sortmpi.h"
#include "paramdb.h"
#include "util.h"


GraphBspl::GraphBspl(const ParamDb& db) :
    m_ilower(-1),
    m_iupper(-1),
    m_nxAux (db.GetInt("Domain_Nx") - db.GetInt("General_BsplDeg")),
    m_nyAux (db.GetInt("Domain_Ny") - db.GetInt("General_BsplDeg")),
    m_nzAux (db.GetInt("Domain_Nz") - db.GetInt("General_BsplDeg")),
    m_deg   (db.GetInt("General_BsplDeg"))
{
}

//
// Generates and partitions the graph (mesh) required for calculations
//
void GraphBspl::Generate(MPI_Comm comm)
{
int rank, npes;
int mapSize;
int *vtxdist;          // initial distribution of verexes in graph
int *partlocal, *part; // partitions
int *maplocal , *map ; // mapping


    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &rank);

    vtxdist = new int[npes + 1];

    Split(comm, vtxdist);
    const int nlocal = vtxdist[rank + 1] - vtxdist[rank];
    const int minM = vtxdist[rank];
    const int maxM = vtxdist[rank + 1] - 1;

    GenPre(comm, minM, maxM);
    partlocal = Partition(comm, vtxdist);

    part = DistFrag(partlocal, nlocal, comm);
    maplocal = GenMap(partlocal, nlocal, minM, mapSize, comm);
    map = DistFrag(maplocal, mapSize, comm);

    GenPost(comm, part, map);


    delete [] partlocal;
    delete [] vtxdist;
}

//
// Splits graph for generating initial partitions (PRE-partitions)
//
void GraphBspl::Split(MPI_Comm comm, int *vtxdist)
{
const int M = m_nxAux * m_nyAux * m_nzAux;
int npes, local_size, extra, iupper;

    MPI_Comm_size(comm, &npes);

    local_size =  M / npes;
    extra = M - local_size * npes;

    vtxdist[0] = 0;
    for(int rank = 1; rank <= npes; rank++)
    {
        iupper = local_size * rank;
        iupper += std::min(rank, extra);
        vtxdist[rank] = iupper;
        // printf("%6d", iupper);
    }
}

//
// Generate graph on initial partitioning
//
void GraphBspl::GenPre(MPI_Comm comm, int minM, int maxM)
{
int idx[3];
int npes, j, ngbNo, dof;
const int w = 2 * m_deg + 1;
std::vector<int> ngb(4 * w * w * w); // Maksimum number of neigbours

    MPI_Comm_size(comm, &npes);
    // MPI_Comm_rank(comm, &rank);

    assert(minM < maxM);
    // printf("(minM, maxM) = (%d, %d)\n", minM, maxM); fflush(stdout);

    Resize(maxM - minM + 1);

    j = 0;
    for(int ix = 0; ix < m_nxAux; ix++)
    {
        for(int iy = 0; iy < m_nyAux; iy++)
        {
            for(int iz = 0; iz < m_nzAux; iz++)
            {
                dof = Dof(ix, iy, iz);

                if(dof < minM || dof > maxM) // Select part of the graph related to this processor
                    continue;

                // printf("%d\n", dof);

                idx[0] = ix;
                idx[1] = iy;
                idx[2] = iz;
                ngbNo = GenNgb(idx, ngb, NULL);

                Copy(j, &(ngb[0]), ngbNo);
                j++;
            }
        }
    }
    // printf("j = %d\n", j);

    Sort();


#if 0
    int rank;
    MPI_Comm_rank(comm, &rank);
    char path[30];
    sprintf(path, "./graph/a-part.%04d", rank);
    WriteTxt(path);
#endif
}


//
// Generates neigbours of B-Splines with index "idx"
// "map" is a mapping array from old-numbering to new-numbering
// If "map==NULL" no re-numbering is made.
//
int GraphBspl::GenNgb(const int* idx, std::vector<int>& ngb, const int* map) const
{
const int aux[3] = {m_nxAux, m_nyAux, m_nzAux};
int k, beg[3], end[3];
int dofJ, dofI;

    for(int i = 0; i < 3; i++)
    {
        beg[i] = std::max(idx[i] - m_deg, 0);
        end[i] = std::min(idx[i] + m_deg, aux[i] - 1);
    }

    dofJ = Dof(idx[0], idx[1], idx[2]);
    if(map)
        dofJ = map[dofJ];

    ngb[0] = dofJ;
    ngb[1] = idx[0];
    ngb[2] = idx[1];
    ngb[3] = idx[2];

    k = 4;
    for(int ix = beg[0]; ix <= end[0]; ix++)
    {
        for(int iy = beg[1]; iy <= end[1]; iy++)
        {
            for(int iz = beg[2]; iz <= end[2]; iz++)
            {
                dofI = Dof(ix, iy, iz);
                if(map)
                    dofI = map[dofI];

                if(dofJ != dofI)
                {
                    assert(static_cast<int>(ngb.size()) > k + 3);

                    // printf("%d ", dofI);

                    ngb[k + 0] = dofI;
                    ngb[k + 1] = ix;
                    ngb[k + 2] = iy;
                    ngb[k + 3] = iz;
                    k += 4;
                }
            }
        }
    }
    // printf("\n");
    return k;
}

//
// Graph partitioning by ParMETIS library
// "vtxdist" - initial distrubution of graph nodes
//
int* GraphBspl::Partition(MPI_Comm comm, int *vtxdist)
{
// int wgtflag = 0; // No weights (vwgt and adjwgt are both NULL).
int wgtflag = 2; // Weights on the vertices only (adjwgt is NULL).
int numflag = 0; // C-style numbering that starts from 0.
int ncon = 1; // Number of weights that each vertex has.
float *tpwgts; // An array  to specify the fraction of vertex weight
float *ubvec; // An array of size ncon that is used to specify the imbalance tolerance
int options[10] = {1, 3, 7}; // additional parameters for the routine
int edgecut; // Upon successful completion, the number of edges that are cut by the partitioning
int *part, nparts;
int *xadj, *adjncy; // CSR (compressed Sparse Row) representation of the graph
int *vwgt; // Weights on the vertices
int npes, rank;



    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &rank);

    const int ni = vtxdist[rank + 1] - vtxdist[rank];

    // Number of partitions generated by ParMETIS is equal to number of processes!
    nparts = npes;

    vwgt = new int[NodeNo()];
    for(int i = 0; i < NodeNo(); ++i)
    {
        vwgt[i] = 1.0;
        // vwgt[i] = EltNo(i) - 1;
    }

    Adj2Csr(xadj, adjncy);

    Clear(); // Graph (in adjacency structure) no more needed

    tpwgts = new float[ncon * nparts];
    for(int i = 0; i < nparts * ncon; ++i)
        tpwgts[i] = 1.0 / (float)nparts;

    ubvec = new float[ncon];
    for(int i = 0; i < ncon; ++i)
        ubvec[i] = 1.05;

    part = new int [ni];

    assert(vtxdist);
    ParMETIS_V3_PartKway(vtxdist, xadj, adjncy,
                         vwgt, NULL, &wgtflag, &numflag, &ncon,
                         &nparts, tpwgts, ubvec, options, &edgecut, part, &comm);

#if 0
    char name[300];
    sprintf(name, "part/part.%04d", mype);
    FILE* out = fopen(name, "wt");
    fprintf(out, "ni = %d\n", ni);
    for(int i = 0; i < ni; ++i)
        fprintf(out, "%d\n", part[i]);
    fclose(out);
#endif

    delete [] vwgt;
    delete [] xadj;
    delete [] adjncy;
    delete [] ubvec;
    delete [] tpwgts;

    return part;
}

//
// Conversion of graph from adjacency representation to CSR representation
//
void GraphBspl::Adj2Csr(int*& xadj, int*& adjncy)
{
    xadj = new int[NodeNo() + 1];
    adjncy = new int[EltNo() - NodeNo()];

    xadj[0] = 0;
    for(int i = 0; i < NodeNo(); i++)
    {
        // The first element must be skipped!
        xadj[i + 1] = xadj[i] + EltNo(i) - 1;
        // printf("%d\n", xadj[i + 1]);
        for(int k = 0, j = xadj[i]; j < xadj[i + 1]; j++, k++)
            adjncy[j] = GetDof(i, k + 1);
            //adjncy[j] = m_node[i][4 * (k + 1)];
    }
    // printf("** %d\n", xadj[NodeNo()]);
}

//
// Distributes array "arrlocal" to all procesors from group "comm".
// Returns pointer to gathered global array;
//
int* GraphBspl::DistFrag(int* arrlocal, int nlocal, MPI_Comm comm)
{
int npes, rank, shift, totn;
int *counts, *displs, *arr;

    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &rank);

    counts = new int[npes];
    displs = new int[npes];

    MPI_Allreduce(&nlocal, &totn, 1, MPI_INT, MPI_SUM, comm);
    arr = new int[totn];

    MPI_Allgather(&nlocal, 1, MPI_INT, counts, 1, MPI_INT, comm);
    MPI_Scan(&nlocal, &shift, 1, MPI_INT, MPI_SUM, comm);
    shift -= nlocal;

    MPI_Allgather(&shift, 1, MPI_INT, displs, 1, MPI_INT, comm);

#if 0
    // if(rank == 0)
    {
        for(int i = 0; i < npes; i++)
            printf("rank=%d, %d: shift=%d, counts=%d, displs=%d\n", rank, i, shift, counts[i], displs[i]);
    }
#endif

    MPI_Allgatherv(arrlocal, nlocal, MPI_INT, arr, counts, displs, MPI_INT, comm);

#if 0
    if(rank == 0)
    {
        printf("--------------------------------\n");
        for(int i = 0; i < totn; i++)
            printf("%d %d\n", i, arr[i]);
        printf("--------------------------------\n");
        fflush(stdout);
    }
    //MPI_Barrier(comm);
#endif

    delete [] counts;
    delete [] displs;

    return arr;
}

//
// Generates maping from DOF-old to DOF-new.
// Returns pointer to global mapping array.
//
int* GraphBspl::GenMap(int* part, int nlocal, int minM, int& mapSize, MPI_Comm comm)
{
PairInt *dataF, *dataS;
int i, sizeF, sizeS, shift;

int rank;

    MPI_Comm_rank(comm, &rank);


    { // local scope for "data"
    PairInt *data = new PairInt[nlocal];
    for(i = 0; i < nlocal; i++)
        data[i] = PairInt(part[i], minM + i); // (part, DOF-old)

    // Sort by the first member, i.e. sort by part
    dataF = sortmpi(nlocal, (PairF*)data, &sizeF, comm);
    delete [] data;
    }

    // WriteData(dataF, sizeF, comm);

    { // local scope for "shift"
    MPI_Scan(&sizeF, &shift, 1, MPI_INT, MPI_SUM, comm);

    shift -= sizeF;
    for(i = 0; i < sizeF; i++)
        dataF[i].m_first = shift + i; //(DOF-new, DOF-old)
    }

    // Sort by the second member, i.e. sort by DOF-old
    dataS = sortmpi(sizeF, (PairS*)dataF, &sizeS, comm);
    delete [] dataF;

    // WriteData(dataS, sizeS, comm);

    // Mapping from DOF-old to DOF-new: dofnew = map[dofold]
    int *map = new int[sizeS];
    for(i = 0; i < sizeS; i++)
        map[i] = dataS[i].m_first;

#if 0
    char name[300];
    sprintf(name, "map/part.%04d", rank);
    FILE* out = fopen(name, "wt");
    fprintf(out, "sizeS = %d\n", sizeS);
    for(int i = 0; i < sizeS; ++i)
        fprintf(out, "%d %d\n", dataS[i].m_second, dataS[i].m_first);
    fclose(out);
#endif

    delete [] dataS;

    mapSize = sizeS;
    return map;
}

//
// Generates graph based on partiotion "part" and mapping "map"
//
void GraphBspl::GenPost(MPI_Comm comm, const int *part, const int *map)
{
int idx[3];
int j, k, rank, npes, ngbNo, dof;
const int w = 2 * m_deg + 1;
std::vector<int> ngb(4 * w * w * w); // Maksimum number of neigbours

    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &rank);


    int cnt = 0;
    const int M = m_nxAux * m_nyAux * m_nzAux;
    for(k = 0; k < M; k++)
    {
        if(part[k] == rank)
            cnt++;
    }
    Resize(cnt);

    assert(map);
    assert(part);
    j = 0;
    for(int ix = 0; ix < m_nxAux; ix++)
    {
        for(int iy = 0; iy < m_nyAux; iy++)
        {
            for(int iz = 0; iz < m_nzAux; iz++)
            {
                dof = Dof(ix, iy, iz);
                assert(dof < M);

                // Select DOFs belonging to part corresponding to proper PE rank
                if(part[dof] != rank)
                    continue;

                dof = map[dof];
                assert(dof < M);

                idx[0] = ix;
                idx[1] = iy;
                idx[2] = iz;
                ngbNo = GenNgb(idx, ngb, map);

                Copy(j, &(ngb[0]), ngbNo);

                j++;
            }
        }
    }
    // printf("rank=%d, GraphBspl::NodeNo=%d\n", rank, NodeNo());

    Sort();

    // Lower and upper limitis in new graph partitions
    int nlocal = NodeNo();
    MPI_Scan(&nlocal, &m_ilower, 1, MPI_INT, MPI_SUM, comm);
    m_ilower -= nlocal;
    m_iupper = m_ilower + nlocal - 1;

#if 0
    // int rank;
    MPI_Comm_rank(comm, &rank);
    char path[30];
    sprintf(path, "./graph/b-part.%04d", rank);
    WriteTxt(path);
#endif

}

