#include <stdio.h>
#include "recmtx.h"

//
// Construcor. Creates empty matrix.
//
RecMtx::RecMtx() : m_colNo(0), m_rowNo(0), m_array(nullptr)
{
}

//
// Constructor. Creates matrix of dimension rowNo and colNo
//
RecMtx::RecMtx(int rowNo, int colNo) : m_array(nullptr)
{
    Malloc(rowNo, colNo);
}

//
// Destructor
//
RecMtx::~RecMtx(void)
{
    delete [] m_array;
}

//
// Allocates matrix of dimension rowNo, colNo
//
void RecMtx::Malloc(int rowNo, int colNo)
{
    delete [] m_array;

    assert(rowNo > 0);
    assert(colNo > 0);

    m_rowNo = rowNo;
    m_colNo = colNo;
    m_array = new double[rowNo * colNo];
    assert(m_array);
}

//
// Returns pointer to column "col".
//
double* RecMtx::Col(int col)
{
    assert(col < m_colNo);
    return m_array + col * m_rowNo;
}

//
// Returns const pointer to column "col".
//
const double* RecMtx::Col(int col) const
{
    assert(col < m_colNo);
    return m_array + col * m_rowNo;
}

//
// Sets all enrties in matrix to zero.
//
void RecMtx::Zero(void)
{
    for(int i = 0; i < m_rowNo * m_colNo; i++)
        m_array[i] = 0;

    // memset(m_array, m_rowNo * m_colNo * sizeof(double), 0);
}

//
// Multiplies column "col" by scalar "alpha"
//
void RecMtx::MulCol(int col, double alpha)
{
double *p = Col(col);

    for(int i = 0; i < m_rowNo; i++)
        p[i] *= alpha;
}


//
// Writes matrix into file "path".
// Format: each line contains (row, col, value)
//
void RecMtx::Write(const char* path) const
{
FILE* out;
int row, col;

    out = fopen(path, "wt");
    for(row = 0; row < RowNo(); row++)
    {
        for(col = 0; col < ColNo(); col++)
            fprintf(out, "(%d %d %E\n)", row, col, (*this)(row, col));
    }
    fprintf(out, "\n");
    fclose(out);
}



