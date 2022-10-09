#ifndef RSCHR_RECMTX__H
#define RSCHR_RECMTX__H

//
// AUTHOR: Zbigniew Romanowski [romz@wp.pl]
//

//
// 1. Rectangular matrix compatible with BLAS, LAPACK.
//
// 2. Column-major order, like in FORTRAN 77
//


#include <cassert>


class RecMtx
{
public:
    RecMtx();
    RecMtx(int rowNo, int colNo);
    ~RecMtx();

    void Malloc(int rowNo, int colNo);

    double  operator()(int row, int col) const { return m_array[Elt(row, col)]; }
    double& operator()(int row, int col)       { return m_array[Elt(row, col)]; }

    int ColNo(void) const { return m_colNo; }
    int RowNo(void) const { return m_rowNo; }

    double*       Col(int col);
    const double* Col(int col) const;

    void MulCol(int col, double alpha);
    void Zero(void);

    void Write(const char* path) const;

private:
    int Elt(int row, int col) const;

private:
    // Number of columns
    int m_colNo;

    // Number of rows
    int m_rowNo;

    // Array of data
    double* m_array;

};

//
// Returns the index of the array.
// Column-major order, like in FORTRAN 77.
//
inline
int RecMtx::Elt(int row, int col) const
{
    assert(row < m_rowNo);
    assert(col < m_colNo);

    return col * m_rowNo + row;
}
#endif

