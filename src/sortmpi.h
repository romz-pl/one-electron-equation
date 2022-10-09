#ifndef RSCHR_SORTMPI__H
#define RSCHR_SORTMPI__H

//
// AUTHOR: Zbigniew Romanowski [romz@wp.pl]
//

//
// 1. Distributed sorting algorithm for array of pairs of integers.
//
// 2. Implemented algorithm: sampling algorithm.
//
// 3. Name of sorting function: sortmpi().
//
// 4. It returns pointer to sorted array.
//
// 5. Sorted array may have different length then input array.
//
// 6. The length of sorted array is returned in argument "bucketSize".
//
// 7. The function "sortmpi" works only for defined types "PairF" and "PairS".
//    "PairF" and "PairS" is a pair of integers.
//
// 8. Type "PairF" sorts by the first integer.
//
// 9. Type "PairS" sorts by the second integer.
//
//10. This procedure is dedicated for RSchr program, howevet it is
//    strightforward to extend its functionality for other problems
//


#include <mpi.h>
#include <cstdlib>
#include <climits>
#include <cassert>
#include <algorithm>


//
// Pair of integers
//
class PairInt
{
public:
    PairInt() : m_first(-1), m_second(-1) { }
    PairInt(int f, int s) : m_first(f), m_second(s) { }

public:
    int m_first;
    int m_second;
};

//
// Sorting by the first memeber
//
class PairF : public PairInt
{
public:
    PairF() : PairInt(-1, -1) { }
    PairF(int f, int s) : PairInt(f, s) { }

public:
    bool operator< (const PairF& a) const
    {
        if(m_first ==  a.m_first)
            return  m_second <  a.m_second;
        return m_first <  a.m_first;
    }

    bool operator<=(const PairF& a) const
    {
        if(m_first ==  a.m_first)
            return  m_second <=  a.m_second;
        return m_first <=  a.m_first;
    }

    bool operator> (const PairF& a) const { return a <  *this; }
    bool operator>=(const PairF& a) const { return a <= *this; }
};

//
// Sorting by the second memeber
//
class PairS : public PairInt
{
public:
    PairS() : PairInt(-1, -1) { }
    PairS(int f, int s) : PairInt(f, s) { }

public:
    bool operator< (const PairS& a) const
    {
        if(m_second ==  a.m_second)
            return  m_first <  a.m_first;
        return m_second <  a.m_second;
    }

    bool operator<=(const PairS& a) const
    {
        if(m_second ==  a.m_second)
            return  m_first <=  a.m_first;
        return m_second <=  a.m_second;
    }

    bool operator> (const PairS& a) const { return a <  *this; }
    bool operator>=(const PairS& a) const { return a <= *this; }
};


//
// ---------------------------------------------------------------------------
//

template <typename T>
T* sortmpi(int nlocal, const T *data, int *bucketSize, MPI_Comm comm);

template <typename T>
T* SampleSort(int nlocal, T *data, int *bucketSize, MPI_Comm comm);


MPI_Datatype CreateType(PairInt* data);


//
//     IT WORKS FOR "T=PairF" or "T=PairS" ONLY!!!
//
// It is strightforward to extent for other types.
// See definition of type PairF and PairS.
//
// Sorts data in MPI framework.
// Sorted data array is unchanged during sorting.
// Returns pointer to sorted array.
// The lenght of sorted array is "bucketSize".
//
// nlocal     [IN]  - number of elements to sort. Length of array "data"
// data       [IN]  - array of elements to sort
// bucketSize [OUT] - number of elements in returned sorted array
// comm       [IN]  - MPI communication handler
//
template <typename T>
T* sortmpi(int nlocal, const T *data, int *bucketSize, MPI_Comm comm)
{
int rank, nproc;
T *sortedData;

    MPI_Comm_size(comm, &nproc);
    MPI_Comm_rank(comm, &rank);

    if(nproc == 1)
    {
//		if(rank ==0)
//		{
//			printf("Sequential sort!\n");
//			fflush(stdout);
//		}
        // Copy to "result" table.
        sortedData = new T[nlocal];
        for(int i = 0; i < nlocal; ++i)
            sortedData[i] = data[i];

        // Sort "in place"
        std::sort(sortedData, sortedData + nlocal);
        *bucketSize = nlocal;
    }
    else
    {
//		if(rank ==0)
//		{
//			printf("MPI-parallel sort!\n");
//			fflush(stdout);
//		}
        sortedData = SampleSort(nlocal, (T*)data, bucketSize, comm);
    }
    return sortedData;
}



//
// Sample Sort - MPI implementation
//
template <typename T>
T* SampleSort(int nlocal, T *data, int *bucketSize, MPI_Comm comm)
{
int i, j, p, rank;
T *sortedData;
T *splitters, *allpicks = NULL;
int *scounts, *sdispls, *rcounts, *rdispls;
MPI_Datatype PairIntMPI;

    PairIntMPI = CreateType(data);

    /* Get communicator-related information */
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    assert(p > 1);

    /* allocate memory for the arrays that will store the splitters */
    splitters = new T[p];
    if (rank == 0)
        allpicks = new T[p * (p - 1)];

    /* sort local array */
    std::sort(data, data + nlocal);

    /* select local p-1 equally spaced elements */
    for(i = 1; i < p; i++)
        splitters[i - 1] = data[i * (nlocal / p)];


    /* gather the samples into the master */
    MPI_Gather(splitters, p - 1, PairIntMPI, allpicks, p - 1, PairIntMPI, 0, comm);

    /* master selects the splitters among allpicks and broadcasts them */
    if (rank == 0)
    {
        std::sort(allpicks, allpicks + p * (p - 1));

        for(i = 0; i < p - 1; i++)
            splitters[i] = allpicks[i * p];
        splitters[p - 1] = T(INT_MAX, INT_MAX);
    }

    /* now the splitters array contains the global splitters */
    MPI_Bcast(splitters, p, PairIntMPI, 0, comm);

    /* compute the number of elements that belong to each bucket */
    scounts = new int[p];
    for(i = 0; i < p; i++)
        scounts[i] = 0;

    j = 0;
    while(data[0] >= splitters[j]) j++;

    for(i = 0; i < nlocal; i++)
    {
        if (data[i] < splitters[j])
        {
            scounts[j]++;
        }
        else
        {
            while(data[i] >= splitters[j]) j++;
            scounts[j]++;
        }
    }

    /* determine the starting location of each bucket's elements in the data array */
    sdispls = new int[p];
    sdispls[0] = 0;
    for(i = 1; i < p; i++)
    {
        sdispls[i] = sdispls[i - 1] + scounts[i - 1];
    }


    /* Perform an all2all communication to inform the corresponding processes */
    /* of the number of elements they are going to receive. */
    /* This information is stored in rcounts array */
    rcounts = new int[p];
    MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm);


    /* Based on rcounts determines where in the local array the data from each processor */
    /* will be stored. This array will store the received elements as well as the final */
    /* sorted sequence.*/
    rdispls = new int[p];
    rdispls[0] = 0;
    for(i = 1; i < p; i++)
    {
        rdispls[i] = rdispls[i - 1] + rcounts[i - 1];
    }

    /* how much data elements I will get */
    *bucketSize = rdispls[p - 1] + rcounts[p - 1];
    sortedData = new T[*bucketSize];

    /* Each process sends and receives the corresponding elements, using the MPI__Alltoallv */
    /* operation. The arrays scounts and sdispls are used to specify the number of elements */
    /* to be sent and where these elements are stored, respectively. The arrays rcounts */
    /* and rdispls are used to specify the number of elements to be received, and where these */
    /* elements will be stored, respectively. */
    MPI_Alltoallv(data, scounts, sdispls, PairIntMPI, sortedData, rcounts, rdispls, PairIntMPI, comm);

    /* perform the final local sort */
    std::sort(sortedData, sortedData + *bucketSize);

    /* cleanup */
    if (rank == 0)
        delete [] allpicks;

    delete [] splitters;
    delete [] scounts;
    delete [] sdispls;
    delete [] rcounts;
    delete [] rdispls;

    return sortedData;
}




//
// Ceates MPI compliant type for "PairInt" used during sorting
//
inline
MPI_Datatype CreateType(PairInt* data)
{
MPI_Datatype PairIntMPI;
MPI_Datatype type[3] = {MPI_INT, MPI_INT, MPI_UB};
int len[3] = {1, 1, 1};
MPI_Aint disp[3], base;

    MPI_Address(data, disp + 0);
    MPI_Address(&(data[0].m_second), disp + 1);
    MPI_Address(data + 1, disp + 2);

    base = disp[0];
    for(int i = 0; i < 3; i++)
        disp[i] -= base;

    MPI_Type_struct(3, len, disp, type, &PairIntMPI);
    MPI_Type_commit(&PairIntMPI);
    return PairIntMPI;
}



#endif

