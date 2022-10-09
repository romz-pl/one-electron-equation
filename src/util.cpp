#include <cassert>
#include <cstdio>
#include <iostream>
#include <mpi.h>
#include "util.h"


//
// Prinnts text "txt" on stdout for processs with rank=0
//
void Printf(int rank, const char* txt)
{
    if(rank == 0)
    {
        printf("%s", txt);
        fflush(stdout);
    }
}

//
// Prints Intro into file "out"
//
void Intro(FILE* out)
{
    assert(out);

    fprintf(out,
        "===============================================================================\n"
        " RRRRR  SSSSS  H   H  RRRRR    Zbigniew ROMANOWSKI                             \n"
        " R   R  S      H   H  R   R                                                    \n"
        " RRRRR  SSSSS  HHHHH  RRRRR    romz@wp.pl                                      \n"
        " R  R       S  H   H  R  R                                                     \n"
        " R   R  SSSSS  H   H  R   R    version: 5.x                                    \n"
        "===============================================================================\n\n" );
    fflush(out);
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------


//
// "Throws" run-time error.
//
void Throw_runtime_error( const std::string& str )
{
    std::cout << "RUNTIME ERROR: " << str << std::endl;
    MPI_Finalize();
    exit(1);

}

//
// "Throws" invalid-argument error.
//
void Throw_invalid_argument( const std::string& str )
{
    std::cout << "INVALID ARGUMENT: " <<  str << std::endl;
    MPI_Finalize();
    exit(1);
}




