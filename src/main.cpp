//
// AUTHOR: Zbigniew Romanowski [romz@wp.pl]
//
/* for "__builtin_expect" */
#ifdef __xlC__
#  include <builtins.h>
#endif

#include <sstream>
#include <fstream>
#include <string>
#include <iostream>

#include <slepceps.h>

#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include "eigsol.h"
#include "paramdb.h"
#include "poten.h"
#include "schr.h"
#include "util.h"

static void RunSolver( const std::string& paramInp, const std::string& potDef, MPI_Comm comm );
static void CreateGraph( Schr* schr, MPI_Comm comm );
static void MatrixAllocation( const Schr* schr, Mat& mtxH, Mat& mtxS, MPI_Comm comm );
static void Assembling( Schr* schr, Mat& mtxH, Mat& mtxS, MPI_Comm comm );
static void Solving( Mat& mtxH, Mat& mtxS, int ilower, int iupper, MPI_Comm comm, const ParamDb& db );


static void InitParamDb( ParamDb& db, const std::string& path, MPI_Comm comm );
static void WriteParamDb( const ParamDb& db );

static void InitPoten( Poten& pot, const std::string& potDefPath, MPI_Comm comm );

static int ReadInFile( const std::string& path, std::string& buf );

template<typename T>
static void BcastFile( const std::string& file, MPI_Comm comm, T& a );


//
// MAIN
//
int main(int argc, char* argv[])
{
time_t beg, end;
int rank, npes;
MPI_Comm comm;
char paramInp[PETSC_MAX_PATH_LEN], potDef[PETSC_MAX_PATH_LEN];
PetscBool flg;

    SlepcInitialize(&argc, &argv, (char*)0, (char*)0);

    PetscOptionsGetString(PETSC_NULL, "-a", paramInp, PETSC_MAX_PATH_LEN, &flg);
    if (!flg)
    {
        fprintf(stdout, "%s", "Missing -a option\n");
        exit(1);
    }

    PetscOptionsGetString(PETSC_NULL, "-b", potDef, PETSC_MAX_PATH_LEN, &flg);
    if (!flg)
    {
        fprintf(stdout, "%s", "Missing -b option\n");
        exit(1);
    }



    // MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);

    if(rank == 0)
    {
        Intro(stdout);
        ::time(&beg);
        fprintf(stdout, "Start of calculation: %s\n", ctime(&beg) );
        fprintf(stdout, "PEs = %d\n\n", npes);
        fflush(stdout);
    }

    ::time(&beg);
    RunSolver(paramInp, potDef, comm);
    ::time(&end);

    if(rank == 0)
    {
        fprintf(stdout, "\n\n End of calculation: %s\n", ctime(&end) );
        printf("\n\n CALCULATION TIME = %d [seconds] = %.2lf [hours]\n", static_cast<int>(end - beg), static_cast<float>(end - beg) / 3600.);
    }

    // MPI_Finalize();
    SlepcFinalize();
    return 0;
}



//
// Generates matrices H and S, as a discretization of one-electron Schroedinger equation
//
void RunSolver( const std::string& paramPath, const std::string& potenPath, MPI_Comm comm )
{
Poten poten;
Schr* schr;
Mat mtxH, mtxS;
ParamDb param;

    InitParamDb(param, paramPath, comm);

    InitPoten(poten, potenPath, comm);

    schr = new Schr(poten, param);
    assert(schr);

    CreateGraph(schr, comm);

    MPI_Barrier(comm); // ?????

    MatrixAllocation(schr, mtxH, mtxS, comm);
    Assembling(schr, mtxH, mtxS, comm);

    // Store needed parameters
    const int ilower = schr->Ilower();
    const int iupper = schr->Iupper();

    delete schr; // This object is no longer needed. It frees memory!

    MPI_Barrier(comm); // This barier is required by PETSc (?)

    Solving(mtxH, mtxS, ilower, iupper, comm, param);
}

//
// Initialize parameters read from file 'paramDbPath'
//
void InitParamDb( ParamDb& db, const std::string& path, MPI_Comm comm )
{
std::string file;
int rank;

    MPI_Comm_rank(comm, &rank);

    Printf(rank, "----INITIALIZATION----========================================================\n");

    if(rank == 0) // Read parameters from the file on the ROOT
    {
        if( ReadInFile(path, file) != 0)
        {
            MPI_Finalize();
            exit(1);
        }
        // printf("RANK=%d,'%s'\n", rank, file.c_str());
    }

    BcastFile(file, comm, db);

    if(rank == 0)
        WriteParamDb(db);
}

//
// Writes information about the parameters
//
void WriteParamDb(const ParamDb& db)
{
const double x0  = db.GetDbl("Domain_X0" );
const double y0  = db.GetDbl("Domain_Y0" );
const double z0  = db.GetDbl("Domain_Z0" );
const double h   = db.GetDbl("Domain_h"  );
const int nx     = db.GetInt("Domain_Nx" );
const int ny     = db.GetInt("Domain_Ny" );
const int nz     = db.GetInt("Domain_Nz" );
const int deg    = db.GetInt("General_BsplDeg");
const int ord    = db.GetInt("General_CubatOrder");

const double xmax = x0 + h * nx;
const double ymax = y0 + h * ny;
const double zmax = z0 + h * nz;

    fprintf(stdout,
        "DOMAIN:\n"
        "    BsplDeg = %d, CubatOrder = %d, CubatNodeNo = %d\n"
        "    h = %E\n"
        "    Nx = %d,    Ny = %d,   Nz = %d\n"
        "    X0 = %E,    Xmax = %E\n"
        "    Y0 = %E,    Ymax = %E\n"
        "    Z0 = %E,    Zmax = %E\n\n",
        deg, ord, ord * ord * ord,
        h,
        nx, ny, nz,
        x0, xmax,
        y0, ymax,
        z0, zmax);

    // const int eigNo = Db_GetInt("Lobpcg_BlockSize");
    // fprintf(stdout, "NUMBER OF SEARCHED EIGENVALUES = %d\n\n", eigNo);

    fflush(stdout);
}

//
// Initializes the potential
//
void InitPoten( Poten& poten, const std::string& path, MPI_Comm comm )
{
int rank;
std::string file;

    MPI_Comm_rank(comm, &rank);

    if(rank == 0) // Read parameters from the file on the ROOT
    {
        if( ReadInFile(path, file) != 0)
        {
            MPI_Finalize();
            exit(1);
        }
        // printf("###### '%s'\n", file.c_str());
    }

    BcastFile(file, comm, poten);



    if(rank == 0)
        poten.Info(stdout);
}

//
// Reading graph
//
void CreateGraph( Schr* schr, MPI_Comm comm )
{
int rank, npes;

    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &rank);

    Printf(rank, "----GRAPH-CREATION-and-PARTITIONING----=======================================\n");

    schr->Generate(comm);

    double nodeNo = schr->NodeNo();
    double eltNo = schr->EltNo();
    double mem;
    double nodeNoTot, eltNoTot, memTot;

    mem = eltNo;                  // Total number of nodes in the graph
    mem = mem * 4. * sizeof(int); // Each node stores four integers
    mem = mem / (1024. * 1024.);  // Change to MB

    MPI_Reduce(&nodeNo, &nodeNoTot, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(&eltNo , &eltNoTot , 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(&mem   , &memTot   , 1, MPI_DOUBLE, MPI_SUM, 0, comm);


    if(rank == 0)
    {
        printf("\n");
        printf("    Number of nodes:\n");
        printf("         Rank0.........%.3E\n", nodeNo);
        printf("         Averg.........%.3E\n", nodeNoTot / npes);
        printf("         Total.........%.3E\n", nodeNoTot);

        printf("    Number of edges:\n");
        printf("         Rank0.........%.3E\n", eltNo - nodeNo);
        printf("         Averg.........%.3E\n", (eltNoTot - nodeNoTot) / npes);
        printf("         Total.........%.3E\n", eltNoTot - nodeNoTot);

        printf("    Memory usage.[MB]:\n");
        printf("         Rank0.........%.3E\n", mem);
        printf("         Averg.........%.3E\n", memTot / npes);
        printf("         Total.........%.3E\n", memTot);
    }

    Printf(rank, "\n\n");
}


//
// Matrix allocation
//
void MatrixAllocation(const Schr* schr, Mat& mtxH, Mat& mtxS, MPI_Comm comm)
{
int rank, npes;
double nnz, rowNo, nnzTot, rowNoTot, mem, memTot; // In order to avoid overflow, 'double' type is used.
double sprFrac, sprFracAvg, sprRow, sprRowAvg;
int offdiag, diag;
double offdiagTmp, diagTmp, offdiagTot, diagTot;

    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &rank);

    Printf(rank, "----MATRIX-ALLOCATION----=====================================================\n");

    rowNo = schr->Iupper() - schr->Ilower() + 1;
    nnz = schr->Alloc(mtxH, mtxS, offdiag, diag, comm);

    // Two matrices H and S
    // It is assumed, that to store one double entry in the matrix,
    // two additional integers are required in CSR format
    mem = 2 * nnz * (sizeof(double) + 2 * sizeof(int)) / (1024 * 1024);


    MPI_Reduce(&nnz  , &nnzTot  , 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(&rowNo, &rowNoTot, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(&mem  , &memTot  , 1, MPI_DOUBLE, MPI_SUM, 0, comm);

    offdiagTmp = offdiag;
    diagTmp = diag;
    MPI_Reduce(&diagTmp   , &diagTot   , 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(&offdiagTmp, &offdiagTot, 1, MPI_DOUBLE, MPI_SUM, 0, comm);


    sprFrac = nnz / (rowNo * rowNoTot);
    sprFracAvg = nnzTot / (rowNoTot * rowNoTot);

    sprRow = nnz / rowNo;
    sprRowAvg = nnzTot / rowNoTot;

    if(rank == 0)
    {
        printf("    Dimension of matrix H or S:\n");
        printf("              Rank0..............%.3E (%.0f)\n", rowNo   , rowNo   );
        printf("              Averg..............%.3E       \n", rowNoTot / npes   );
        printf("              Total..............%.3E (%.0f)\n", rowNoTot, rowNoTot);

        printf("    Nonzero entries in H or S:\n");
        printf("              Rank0..............%.3E (%.0f)\n", nnz   , nnz   );
        printf("              Averg..............%.3E       \n", nnzTot / npes );
        printf("              Total..............%.3E (%.0f)\n", nnzTot, nnzTot);

        printf("    Sparseness of H or S (fraction):\n");
        printf("              Rank0..............%.3E\n", sprFrac   );
        printf("              Averg..............%.3E\n", sprFracAvg);

        printf("    Sparseness of H or S (non-zero elemnts in a row):\n");
        printf("              Rank0..............%.1lf\n", sprRow   );
        printf("              Averg..............%.1lf\n", sprRowAvg);

/*
        printf("    Off-diagonal entries of H or S:\n");
        printf("              Rank0..............%.3E (%.0f)\n", offdiagTmp, offdiagTmp);
        printf("              Averg..............%.3E       \n", offdiagTot / npes     );
        printf("              Total..............%.3E (%.0f)\n", offdiagTot, offdiagTot);

        printf("    Diagonal entries of H or S:\n");
        printf("              Rank0..............%.3E (%.0f)\n", diagTmp, diagTmp);
        printf("              Averg..............%.3E       \n", diagTot / npes  );
        printf("              Total..............%.3E (%.0f)\n", diagTot, diagTot);
*/
        printf("    Memory usage for H and S [MiB]:\n");
        printf("              Rank0..............%.3E\n", mem);
        printf("              Averg..............%.3E\n", memTot / npes);
        printf("              Total..............%.3E\n", memTot);
    }


    //Printf(rank, "done\n\n");
    Printf(rank, "\n\n");
}

//
// Assembling
//
void Assembling( Schr* schr, Mat& mtxH, Mat& mtxS, MPI_Comm comm )
{
int rank;

    MPI_Comm_rank(comm, &rank);

    Printf(rank, "----ASSEMBLING----============================================================\n");

    time_t beg, end;
    if(rank == 0)
        ::time(&beg);

    schr->Assem(mtxH, mtxS, comm);

    Printf(rank, "done\n");

    if(rank == 0)
    {
        ::time(&end);
        printf("Assembling time = %d [s]\n\n", static_cast<int>(end - beg));
        fflush(stdout);
    }

    Printf(rank, "\n\n");
}

//
// Solving
//
void Solving( Mat& mtxH, Mat& mtxS, int ilower, int iupper, MPI_Comm comm, const ParamDb& db )
{
int rank;

    MPI_Comm_rank(comm, &rank);

    Printf(rank, "----SOLVING----===============================================================\n");

    EigSol eigSol(db);
    eigSol.Run(mtxH, mtxS, ilower, iupper, comm);

    if(rank == 0)
    {
        // const double alpha = Db_GetDbl("General_Alpha");

        fprintf(stdout, "%s", "\n\n");
        fprintf(stdout, "%s", "             EIGENVALUES            ERROR\n");
        fprintf(stdout, "%s", "====================================================\n");

        for(int i = 0; i < eigSol.EigValNo(); i++)
            fprintf(stdout, "%5d %20.10E %20.10E\n", i, eigSol.EigVal(i), eigSol.EigErr(i));
    }


    // Printf(rank, "\ndone\n\n");
}

//
// Read in the file "path" into the buffer "buf"
//
int ReadInFile( const std::string& path, std::string& buf )
{
    std::ifstream in( path.c_str() );
    if( !in )
    {
        std::cout << "Cannot open file for reading. Path = " << path << std::endl;
        return 1;
    }

    std::stringstream strStream;
    strStream << in.rdbuf();  // read the file
    buf = strStream.str();    // buf holds the content of the file


    return 0;
}

//
// Bradcasts the file to other PEs and reads the parameters into "aa" object of type "T"
//
template<typename T>
void BcastFile(const std::string& file, MPI_Comm comm, T& aa)
{
int rank, len;
char *buf;

    MPI_Comm_rank(comm, &rank);

    if(rank == 0)
        len = file.size();

    MPI_Bcast(&len, 1, MPI_INT, 0, comm);

    buf = new char[len + 1];
    if(rank == 0)
        strcpy(buf, file.c_str());

    MPI_Bcast(buf, len, MPI_CHAR, 0, comm);
    // printf("RANK=%d,'%s'\n", rank, buf);

    aa.Read(buf, len);

    delete [] buf;

}

