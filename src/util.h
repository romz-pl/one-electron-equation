#ifndef RSCHR_UTIL__H
#define RSCHR_UTIL__H

//
// AUTHOR: Zbigniew Romanowski [romz@wp.pl]
//

//
// A few utility functions
//


#include <cstdio>
#include <string>


void Printf(int rank, const char* txt);

void Intro(FILE* out);

void Throw_runtime_error( const std::string& str );
void Throw_invalid_argument( const std::string& str );

#endif

