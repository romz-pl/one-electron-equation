#include <cassert>
#include "schr.h"
#include "paramdb.h"
#include "util.h"


//
// Constructor
//
Schr::Schr(Poten& pot, const ParamDb& db) : PotElt(pot, db), OvrKinElt(db), m_db(db)
{
}


//
// Allocates spce for matrices H and S.
//
int Schr::Alloc(Mat& mtxH, Mat& mtxS, int& offdiagAll, int& diagAll, MPI_Comm comm) const
{
const int minM = Ilower();
const int maxM = Iupper();
const int n = maxM - minM + 1;
int nnz;
// const int w = 2 * Db_GetInt("General_BsplDeg") + 1;
// const int w = 2 * 3 + 1;
// const int wmax = w * w * w;
PetscErrorCode ierr;

/*
    // printf("(minM, maxM) = (%d, %d)\n", minM, maxM); fflush(stdout);

    // Create the matrix.
    // Note that this is a square matrix, so we indicate the row partition
    // size twice (since number of rows = number of cols)
    HYPRE_IJMatrixCreate(comm, minM, maxM, minM, maxM, &mtxH);
    HYPRE_IJMatrixCreate(comm, minM, maxM, minM, maxM, &mtxS);

    // Choose a parallel csr format storage
    HYPRE_IJMatrixSetObjectType(mtxH, HYPRE_PARCSR);
    HYPRE_IJMatrixSetObjectType(mtxS, HYPRE_PARCSR);

    nnz = DefPattern(mtxH, mtxS, offdiagAll, diagAll);

    // Initialize before setting coefficients
    HYPRE_IJMatrixInitialize(mtxH);
    HYPRE_IJMatrixInitialize(mtxS);

    return nnz;
*/


    ierr = MatCreate(PETSC_COMM_WORLD, &mtxH);
    CHKERRQ(ierr);

    ierr = MatCreate(PETSC_COMM_WORLD, &mtxS);
    CHKERRQ(ierr);

    ierr = MatSetSizes(mtxH, n, n, PETSC_DETERMINE, PETSC_DETERMINE);
    CHKERRQ(ierr);

    ierr = MatSetSizes(mtxS, n, n, PETSC_DETERMINE, PETSC_DETERMINE);
    CHKERRQ(ierr);

    ierr = MatSetType(mtxH, MATAIJ);
    CHKERRQ(ierr);

    ierr = MatSetType(mtxS, MATAIJ);
    CHKERRQ(ierr);

    //ierr = MatMPIAIJSetPreallocation(mtxH, wmax, PETSC_NULL, wmax, PETSC_NULL); CHKERRQ(ierr);
    //ierr = MatMPIAIJSetPreallocation(mtxS, wmax, PETSC_NULL, wmax, PETSC_NULL); CHKERRQ(ierr);

    //ierr = MatSeqAIJSetPreallocation(mtxH, wmax, PETSC_NULL); CHKERRQ(ierr);
    //ierr = MatSeqAIJSetPreallocation(mtxS, wmax, PETSC_NULL); CHKERRQ(ierr);

    nnz = DefPattern(mtxH, mtxS, offdiagAll, diagAll);

    ierr = MatSetUp(mtxH);
    CHKERRQ(ierr);

    ierr = MatSetUp(mtxS);
    CHKERRQ(ierr);

    return nnz;

}


//
// Defines pattern of sparse matrices.
// It speed up assembling process sagnificantly!
//
int Schr::DefPattern(Mat& mtxH, Mat& mtxS, int& offdiagAll, int& diagAll) const
{
const int minM = Ilower();
const int maxM = Iupper();

int rowNo = maxM - minM + 1;
int nnz = 0;  // Number of nonzero elements
int *row; // contains estimated sizes for each row on this process.
int *diag, *offdiag; // Contain estimated sizes for each row of the diagonal and off-diagonal blocks, respectively

    row     = new int[rowNo];
    diag    = new int[rowNo];
    offdiag = new int[rowNo];

    diagAll    = 0;
    offdiagAll = 0;
    for(int m = 0; m < NodeNo(); m++)
    {
#ifndef NDEBUG
        const int dofI = GetDof(m);
        assert(dofI == m + minM); // DOFs must be consequtive
        assert(dofI >= minM && dofI <= maxM); // DOFs must be given in specified range
#endif

        row[m]     = 0;
        diag[m]    = 0;
        offdiag[m] = 0;

        for(int j = 0; j < EltNo(m); j++)
        {
            row[m]++;
            const int dofJ = GetDof(m, j);

            if(dofJ < minM || dofJ > maxM)
                offdiag[m]++;
            else
                diag[m]++;
            nnz++;
        }
        diagAll    += diag[m];
        offdiagAll += offdiag[m];

    }
    // printf("(diag,offdiag) = (%d, %d)\n", diagAll, offdiagAll);

    //HYPRE_IJMatrixSetRowSizes(mtxH, row);
    //HYPRE_IJMatrixSetRowSizes(mtxS, row);

    //HYPRE_IJMatrixSetDiagOffdSizes(mtxH, diag, offdiag);
    //HYPRE_IJMatrixSetDiagOffdSizes(mtxS, diag, offdiag);


    PetscErrorCode ierr;
    ierr = MatMPIAIJSetPreallocation(mtxH, 0, diag, 0, offdiag);
    CHKERRQ(ierr);

    ierr = MatMPIAIJSetPreallocation(mtxS, 0, diag, 0, offdiag);
    CHKERRQ(ierr);



    for(int m = 0; m < rowNo; m++)
        diag[m] += offdiag[m];

    ierr = MatSeqAIJSetPreallocation(mtxH, 0, diag);
    CHKERRQ(ierr);

    ierr = MatSeqAIJSetPreallocation(mtxS, 0, diag);
    CHKERRQ(ierr);


    delete [] diag;
    delete [] offdiag;

    return nnz;
}



//
// Calculates matrix elements
//
int Schr::Assem(Mat& mtxH, Mat& mtxS, MPI_Comm comm)
{
// const double alpha = Db_GetDbl("General_Alpha");
const int w = 2 * m_db.GetInt("General_BsplDeg") + 1;
const int W = w * w * w; // The maximum number of nonzero elemenst in one row
const int nodeNo = NodeNo();

double *ss = new double[W];
double *hh = new double[W];
int* cols = new int[W];

int nnz;
double ovr, kin;

PetscErrorCode ierr;

#define DOTSNO 50
int rank, k, dots[DOTSNO];

    MPI_Comm_rank(comm, &rank);

    for(k = 0; k < DOTSNO; ++k)
        dots[k] = k * (nodeNo / DOTSNO);

    k = 0;
    for(int m = 0; m < nodeNo; m++)
    {
        const int dofI = GetDof(m, 0);
        const int *idx = GetIdx(m, 0);
        nnz = 0;

        for(int j = 0; j < EltNo(m); j++)
        {
            const int dofJ = GetDof(m, j);
            const int *jdx = GetIdx(m, j);

            assert(nnz < W);

            OvrKin(idx, jdx, ovr, kin);
            const double pot = Pot(idx, jdx);

            // Properly chosen "alpha" parameter guarantees that matrix H is positive definite!
            // const double aux = kin + pot + alpha * ovr;

            cols[nnz] = dofJ;
            ss[nnz] = ovr;
            hh[nnz] = kin + pot;
            nnz++;

        }
        // This is NOT correct!!!
        //HYPRE_IJMatrixSetValues(mtxH, 1, &nnz, &dofI, cols, hh);
        //HYPRE_IJMatrixSetValues(mtxS, 1, &nnz, &dofI, cols, ss);


        // HYPRE_IJMatrixAddToValues(mtxH, 1, &nnz, &dofI, cols, hh);
        // HYPRE_IJMatrixAddToValues(mtxS, 1, &nnz, &dofI, cols, ss);

        ierr = MatSetValues(mtxH, 1, &dofI, nnz, cols, hh, INSERT_VALUES);
        CHKERRQ(ierr);

        ierr = MatSetValues(mtxS, 1, &dofI, nnz, cols, ss, INSERT_VALUES);
        CHKERRQ(ierr);

        if(m == dots[k])
        {
            ++k;
            Printf(rank, ".");
            if(k % 10 == 0)
            {
                char str[10];
                sprintf(str, "%d", k);
                Printf(rank, str);
            }
        }
    }

    // Assemble after setting the coefficients
    // HYPRE_IJMatrixAssemble(mtxH);
    // HYPRE_IJMatrixAssemble(mtxS);

    ierr = MatAssemblyBegin(mtxH, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    ierr = MatAssemblyBegin(mtxS, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    ierr = MatAssemblyEnd(mtxH, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    ierr = MatAssemblyEnd(mtxS, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    ierr = MatSetOption(mtxH, MAT_SYMMETRIC, PETSC_TRUE);
    CHKERRQ(ierr);

    ierr = MatSetOption(mtxS, MAT_SPD, PETSC_TRUE);
    CHKERRQ(ierr);

    // ierr = PetscPrintf(PETSC_COMM_WORLD, "\nLoading completed.\n\n"); CHKERRQ(ierr);

    delete [] ss;
    delete [] hh;
    delete [] cols;

    FreePotCache();

    /*
    MatPartitioning part;
    IS              is;

    MatPartitioningCreate(comm, &part);
    MatPartitioningSetAdjacency(part, mtxH);
    // MatPartitioningSetFromOptions(part);
    MatPartitioningSetType(part, MATPARTITIONINGPARMETIS);
    MatPartitioningApply(part, &is);
    // ISView(is, PETSC_VIEWER_STDOUT_WORLD);
    ISDestroy(&is);
    MatPartitioningDestroy(&part);


    MatPartitioningCreate(comm, &part);
    MatPartitioningSetAdjacency(part, mtxS);
    // MatPartitioningSetFromOptions(part);
    MatPartitioningSetType(part, MATPARTITIONINGPARMETIS);
    MatPartitioningApply(part, &is);
    // ISView(is, PETSC_VIEWER_STDOUT_WORLD);
    ISDestroy(&is);
    MatPartitioningDestroy(&part);
    */

    return 0;
}

