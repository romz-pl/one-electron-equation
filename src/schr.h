#ifndef RSCHR_SCHR__H
#define RSCHR_SCHR__H

//
// AUTHOR: Zbigniew Romanowski [romz@wp.pl]
//

//
// 1. Allocated and evaluaest matrices H and S for generalized eigenvalue problem
//      H x = \lambda S x
//
// 2. Matrices H and S are sparse from HYPRE library
//
// 3. For speeding up assembling proces, the patern for matrics H and S are evaluated
//    before assembling.
//
// 4. Matrices H and S are distrubuted over all processes in the group.
//
// 5. Only local part (i.e. assigne to process) are evaluated.
//
// 6. Evaluation of matrix elements is based on the graph creted and partitioned by ParMETIS.
//



/* for "__builtin_expect" */
#ifdef __xlC__
#  include <builtins.h>
#endif

#include <slepceps.h>
#include "poten.h"
#include "potelt.h"
#include "ovrkinelt.h"
#include "paramdb.h"



class Schr : public PotElt, public OvrKinElt
{
public:
    Schr(Poten& pot, const ParamDb& db);
    ~Schr() = default;

    int Alloc(Mat& mtxH, Mat& mtxS, int& offdiagAll, int& diagAll, MPI_Comm comm) const;
    int Assem(Mat& mtxH, Mat& mtxS, MPI_Comm comm);

private:

    int DefPattern(Mat& mtxH, Mat& mtxS, int& offdiagAll, int& diagAll) const;

private:
    const ParamDb& m_db;

};

#endif

