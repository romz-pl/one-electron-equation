#ifndef RSCHR_EIGSOL__H
#define RSCHR_EIGSOL__H


//
// AUTHOR: Zbigniew Romanowski [romz@wp.pl]
//

//
// Findes a few the lowest eigenvalues of generalized eigenvalue problem
//     H x = \lambda S x
//
// Matrices H, S are represented as sprse matrices from HYPRE library.
//
// For solving generalized eigenvalue problem LOBPCG algorithm is from SLEPc
//



/* for "__builtin_expect" */
#ifdef __xlC__
#  include <builtins.h>
#endif

#include <slepceps.h>
#include <vector>
#include "paramdb.h"

class EigSol
{
public:
    explicit EigSol(const ParamDb& db);
    ~EigSol();

    int Run(Mat& mtxH, Mat& mtxS, int ilower, int iupper, MPI_Comm comm);
    // void WriteRes(int rank) const;

    int EigValNo() const;
    double EigVal(int i) const;
    double EigErr(int i) const;

private:
    //void CreatePrecon(HYPRE_Solver& precon, MPI_Comm comm);
    //void UsePrecon(HYPRE_Solver& solver, HYPRE_Solver& precon);
    //void DestroyPrecon(HYPRE_Solver& precon);

    int ConfigEPS(EPS& eps, Mat mtxH, Mat mtxS);
    int ConfigST(ST& st);
    int ConfigKSP(KSP& ksp);
    int ConfigPC(PC& pc);

    int ConfigASM(PC& pc);
    int ConfigGAMG(PC& pc);

    int WriteInfo(EPS eps);
    int SaveEigen(EPS eps);

private:
    // Calculated eigenvalues
    std::vector< double > m_eigenvalues;


    // Errors of calculated eigenvalues
    std::vector< double > m_err;

    const ParamDb& m_db;
};

#endif

