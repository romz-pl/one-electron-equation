#ifndef RSCHR_GRAPH__H
#define RSCHR_GRAPH__H

//
// AUTHOR: Zbigniew Romanowski [romz@wp.pl]
//

//
// Graph represented as an adjacency list
//
// 1. Each node of the graph consist of four integers: (dof, ix, iy, iz).
//
// 2. Member "dof" is degree of freedom (DOF) of three dimensional basis
//    function, i.e. tensor product of B-Splines.
//
// 3. Triple (ix,iy,iz) defines identyfier related to tensor product of B-Splines.
//
// 4. The graph can be read from and wrritiin to file.
//


#include <vector>

class Graph
{
public:
    Graph() = default;

    int NodeNo() const;

    int EltNo() const;
    int EltNo(int nodeId) const;

    const int* GetIdx(int nodeId, int ngb) const;
    int        GetDof(int nodeId, int ngb) const;

    void Sort();

    void Clear();
    void Resize(int s);
    void Copy(int nodeId, const int* a, int eltNo);

#if 0
    void ReadTxt (const char* path);
    void ReadBin (const char* path);
    void WriteTxt(const char* path);
#endif

private:
    static bool CmpFun( const std::vector<int>& a, const std::vector<int>& b );


private:
    // Array of neigbours. Adjacency structure
    std::vector< std::vector< int > > m_node;
};



#endif
