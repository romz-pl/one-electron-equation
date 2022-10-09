#include <cassert>
#include <mpi.h>
#include <petscpcmg.h>
#include "eigsol.h"
#include "util.h"


//
// Constructor.
//
EigSol::EigSol(const ParamDb& db) : m_db(db)
{
}

//
// Destructor
//
EigSol::~EigSol(void)
{
    // hypre_TFree(m_eigenvalues);
    //free(m_eigenvalues);
}

//
// Calculates a few the lowest eigenvalues for generalized eigenvalue problem H x = \lambda S x
//
int EigSol::Run(Mat& mtxH, Mat& mtxS, int ilower, int iupper, MPI_Comm comm)
{
    PetscErrorCode ierr;
    EPS eps;
    ST st;
    KSP ksp;
    PC pc;

    ConfigEPS(eps, mtxH, mtxS);

    ierr = EPSGetST(eps, &st);
    CHKERRQ(ierr);
    ConfigST(st);

    ierr = STGetKSP(st, &ksp);
    CHKERRQ(ierr);
    ConfigKSP(ksp);

    ierr = KSPGetPC(ksp, &pc);
    CHKERRQ(ierr);
    ConfigPC(pc);

    // KSPSetUp(ksp);

    ierr = EPSView(eps, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "\n\n");

    // Solve the eigensystem
    ierr = EPSSolve(eps);
    CHKERRQ(ierr);

    WriteInfo(eps);


    SaveEigen(eps);

    ierr = EPSDestroy(&eps);
    CHKERRQ(ierr);

    return 0;
}


int EigSol::ConfigEPS(EPS& eps, Mat mtxH, Mat mtxS)
{
    PetscErrorCode ierr;
    const PetscInt nev     = m_db.GetInt("Slepc_eps_nev");
    const PetscInt ncv     = m_db.GetInt("Slepc_eps_ncv");
    const PetscInt mpd     = m_db.GetInt("Slepc_eps_mpd");
    const PetscInt maxit   = m_db.GetInt("Slepc_eps_max_it");
    const PetscReal tol    = m_db.GetDbl("Slepc_eps_tol");
    const std::string type = m_db.GetStr("Slepc_eps_type");


    ierr = EPSCreate(PETSC_COMM_WORLD, &eps); // Create eigensolver context
    CHKERRQ(ierr);

    ierr = EPSSetOperators(eps, mtxH, mtxS);
    CHKERRQ(ierr);

    ierr = EPSSetProblemType(eps, EPS_GHEP);
    CHKERRQ(ierr);

    ierr = EPSSetDimensions(eps, nev, ncv, mpd);
    CHKERRQ(ierr);

    ierr = EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL);
    CHKERRQ(ierr);

    ierr = EPSSetTolerances(eps, tol, maxit);
    CHKERRQ(ierr);

    ierr = EPSSetType(eps, type.c_str());
    CHKERRQ(ierr);

    ierr = EPSMonitorSet(eps, EPSMonitorFirst, PETSC_NULL, PETSC_NULL);
    CHKERRQ(ierr);

    ierr = EPSSetUp(eps);
    CHKERRQ(ierr);

    //   ierr = EPSView(eps, PETSC_VIEWER_STDOUT_WORLD);  CHKERRQ(ierr);

    return 0;

}

int EigSol::ConfigST(ST& st)
{
PetscErrorCode ierr;

    ierr = STSetMatStructure(st, SAME_NONZERO_PATTERN);
    CHKERRQ(ierr);

    return 0;
}

int EigSol::ConfigKSP(KSP& ksp)
{
    PetscErrorCode ierr;
    const PetscReal rtol   = m_db.GetDbl("Slepc_st_ksp_rtol");
    const std::string type = m_db.GetStr("Slepc_st_ksp_type");


    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    CHKERRQ(ierr);

    ierr = KSPSetInitialGuessKnoll(ksp, PETSC_TRUE);
    CHKERRQ(ierr);

    ierr = KSPSetTolerances(ksp, rtol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    CHKERRQ(ierr);

    // KSPGMRES
    // KSPMINRES
    // KSPCG
    ierr = KSPSetType(ksp, type.c_str());
    CHKERRQ(ierr);


    ierr = KSPSetUp(ksp);
    CHKERRQ(ierr);

    return 0;

}

int EigSol::ConfigPC(PC& pc)
{
    PetscErrorCode ierr;
    const std::string pc_type = m_db.GetStr("Slepc_st_pc_type");


    ierr = PCSetType(pc, pc_type.c_str());
    CHKERRQ(ierr);

    if(pc_type == "asm")
        ConfigASM(pc);


    // if(pc_type == "gamg")
    //	ConfigGAMG(pc);

    ierr = PCSetUp(pc);
    CHKERRQ(ierr);

    return 0;
}


int EigSol::ConfigGAMG(PC& pc)
{
    const std::string typeStr      = m_db.GetStr("Slepc_gamg_type");
    const std::string cycleTypeStr = m_db.GetStr("Slepc_gamg_cycle_type");
    const int smoothUp             = m_db.GetInt("Slepc_gamg_smooth_up");
    const int smoothDown           = m_db.GetInt("Slepc_gamg_smooth_down");
    const int levels               = m_db.GetInt("Slepc_gamg_levels");

    PetscErrorCode ierr;
    PCMGType type;
    PCMGCycleType cycleType;


    ierr = PCMGSetLevels(pc, levels, PETSC_NULL);
    CHKERRQ(ierr);


    if(typeStr == "MULTIPLICATIVE")
    {
        type = PC_MG_MULTIPLICATIVE;
    }
    else if (typeStr == "ADDITIVE")
    {
        type = PC_MG_ADDITIVE;
    }
    else if (typeStr == "FULL")
    {
        type = PC_MG_FULL;
    }
    else if (typeStr == "KASKADE")
    {
        type = PC_MG_KASKADE;
    }
    else
    {
        PetscPrintf(PETSC_COMM_WORLD, "Unknown PCMGType: '%s'\n", typeStr.c_str());
        SlepcFinalize();
        exit(1);
    }
    ierr = PCMGSetType(pc, type);
    CHKERRQ(ierr);



    if(cycleTypeStr == "V")
    {
        cycleType = PC_MG_CYCLE_V;
    }
    else if (cycleTypeStr == "W")
    {
        cycleType = PC_MG_CYCLE_W;
    }
    else
    {
        PetscPrintf(PETSC_COMM_WORLD, "Unknown PCMGCycleType: '%s'\n", cycleTypeStr.c_str());
        SlepcFinalize();
        exit(1);
    }
    ierr = PCMGSetCycleType(pc, cycleType);
    CHKERRQ(ierr);


    ierr = PCMGSetNumberSmoothUp(pc, smoothUp);
    CHKERRQ(ierr);

    ierr = PCMGSetNumberSmoothDown(pc, smoothDown);
    CHKERRQ(ierr);

    return 0;
}


int EigSol::ConfigASM(PC& pc)
{
    const int overlap            = m_db.GetInt("Slepc_asm_overlap");
    const int local_subdomains   = m_db.GetInt("Slpec_asm_local_subdomains");
    const std::string typeStr    = m_db.GetStr("Slepc_asm_type");

    PetscErrorCode ierr;
    PCASMType type;


    if(typeStr == "BASIC")
    {
        type = PC_ASM_BASIC;
    }
    else if (typeStr == "RESTRICT")
    {
        type = PC_ASM_RESTRICT;
    }
    else if (typeStr == "INTERPOLATE")
    {
        type = PC_ASM_INTERPOLATE;
    }
    else if (typeStr == "NONE")
    {
        type = PC_ASM_NONE;
    }
    else
    {
        PetscPrintf(PETSC_COMM_WORLD, "Unknown PCASMType: '%s'\n", typeStr.c_str());
        SlepcFinalize();
        exit(1);
    }

    ierr = PCASMSetType(pc, type);
    CHKERRQ(ierr);

    ierr = PCASMSetOverlap(pc, overlap);
    CHKERRQ(ierr);

    ierr = PCASMSetLocalSubdomains(pc, local_subdomains, PETSC_NULL, PETSC_NULL);
    CHKERRQ(ierr);

    return 0;
}



int EigSol::WriteInfo(EPS eps)
{
    PetscErrorCode ierr;
    const EPSType type;
    PetscInt its, lits, nev, maxit;
    PetscReal tol;

    // Optional: Get some information from the solver and display it
    ierr = EPSGetIterationNumber(eps, &its);
    CHKERRQ(ierr);

    ierr = EPSGetOperationCounters(eps, PETSC_NULL, PETSC_NULL, &lits);
    CHKERRQ(ierr);

    ierr = EPSGetType(eps, &type);
    CHKERRQ(ierr);

    ierr = EPSGetDimensions(eps, &nev, PETSC_NULL, PETSC_NULL);
    CHKERRQ(ierr);

    ierr = EPSGetTolerances(eps, &tol, &maxit);
    CHKERRQ(ierr);


    ierr = PetscPrintf(PETSC_COMM_WORLD, "\n"); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------\n"); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "    Number of iterations of the method.............%D\n", its);   CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "    Number of linear iterations of the method......%D\n", lits);  CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "    Solution method................................%s\n", type);  CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "    Number of requested eigenvalues................%D\n", nev);   CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "    Stopping condition: eps_tol....................%E\n", tol); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "    Stopping condition: eps_max_it.................%D\n", maxit); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------\n"); CHKERRQ(ierr);

    // Display solution and clean up
    // ierr = EPSPrintSolution(eps, PETSC_NULL); CHKERRQ(ierr);

    return 0;
}


int EigSol::SaveEigen(EPS eps)
{
    PetscInt nconv;
    PetscScalar kr, ki;
    // PetscReal err;

    EPSGetConverged(eps, &nconv);
    m_eigenvalues.resize(nconv);
    m_err.resize(nconv);

    for(int j = 0; j < nconv; j++)
    {
        EPSGetEigenpair(eps, j, &kr, &ki, PETSC_NULL, PETSC_NULL);
        m_eigenvalues[j] = kr;

        // This function does not work on IBM Blue Gene/P NOTOS
        // EPSComputeRelativeError(eps, j, &err);
        // m_err[j] = err;
        m_err[j] = 0;
    }

    return 0;
}


//
// Returns number of calculated eigenvalues
//
int EigSol::EigValNo() const
{
    // return Db_GetInt("Lobpcg_BlockSize");
    return static_cast< int >( m_eigenvalues.size() );
}


//
// Returns i-th eigenvalue
//
double EigSol::EigVal(int i) const
{
    return m_eigenvalues[i];
}

//
// Returns error of i-th eigenvalue
//
double EigSol::EigErr(int i) const
{
    return m_err[i];
}

