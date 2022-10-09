#include <cassert>
#include <algorithm>
#include "graph.h"


//
// Returns number of elements for "nodeId".
// Number of neigbours is one less.
// "nodeId" start from 0. It is local numbering for distributed graph.
//
int Graph::EltNo(int nodeId) const
{
    assert(nodeId >= 0 && nodeId < NodeNo());

    // Each element occupies 4 integers
    return static_cast<int>(m_node[nodeId].size() / 4);
}

//
// Returns total number of elements in graph.
//
int Graph::EltNo() const
{
    int sum = 0;

    for(int i = 0; i < NodeNo(); ++i)
        sum += EltNo(i);
    return sum;
}

//
// Returns number of nodes in graph.
//
int Graph::NodeNo() const
{
    return static_cast< int >( m_node.size() );
}

//
// For node "nodeId" and neighbour "ngb",
// it returns pointer to array int[3] holding indexes (ix,iy,iz) of B-splines
// creating tensor product.
//
const int* Graph::GetIdx(int nodeId, int ngb) const
{
    // Each node consists of four integers, hence multiplication by 4.
    // The triple (ix,iy,iz) starts from element "1", hence addition of 1.
    const int k = 4 * ngb + 1;

    assert(nodeId >= 0 && nodeId < NodeNo());
    assert(k >= 0 && k < (int)m_node[nodeId].size());

    return &(m_node[nodeId][k]);
}

//
// Returns DOF of B-spline.
// B-spline is defined by "nodeId" and neighbours "ngb".
//
int Graph::GetDof(int nodeId, int ngb) const
{
    // Each node consists of four integers, hence multiplication by 4.
    // The DOF starts from element "0".
    const int k = 4 * ngb;

    assert(nodeId >= 0 && nodeId < NodeNo());
    assert(k >= 0 && k < (int)m_node[nodeId].size());

    return m_node[nodeId][k];
}


//
// Comparison function
// First elements of vectors are compared.
// First elements of vectors are DOF-s of the nodes.
// Hence, sorting is by DOF of nodes of the graph.
//
bool Graph::CmpFun(const std::vector<int>& a, const std::vector<int>& b)
{
    return (a[0] < b[0]);
}

//
// Sorts adjacency graph.
// Sorting is by DOF of nodes of the graph
//
void Graph::Sort(void)
{
    std::sort(m_node.begin(), m_node.end(), CmpFun);
}

//
// Clears graph
//
void Graph::Clear(void)
{
    m_node.clear();
}

//
// Resizes the number of nodes in graph
//
void Graph::Resize(int s)
{
    assert(s > 0);
    m_node.resize(s);
}

//
// Copies array "a" to node "nodeId".
// Array "a" has "eltNo" members.
// "eltNo" must be divisible by four.
//
void Graph::Copy(int nodeId, const int* a, int eltNo)
{
    assert(nodeId < (int)m_node.size());
    assert(a);
    assert(eltNo % 4 == 0);

    m_node[nodeId].resize(eltNo);
    for(int i = 0; i < eltNo; ++i)
        m_node[nodeId][i] = a[i];
}


#if 0
//
// Reads graph from file (TEXT mode)
//
void Graph::ReadTxt(const char* path)
{
FILE* in;
int nodeNo, eltNo, dof;

    in = fopen(path, "rt");
    if(!in)
    {
        printf("Cannot open file '%s' for reading.\n", path);
        exit(1);
    }

    fscanf(in, "%d %d %d %d %d %d\n", &m_nx, &m_ny, &m_nz, &m_deg, &m_ilower, &m_iupper);
    nodeNo = m_iupper - m_ilower + 1;
    m_node.resize(nodeNo);

    for(int i = 0; i < nodeNo; ++i)
    {
        fscanf(in, "%d ", &eltNo);
        // printf("%d: ", eltNo);
        m_node[i].resize(4 * eltNo);
        for(int j = 0; j < 4 * eltNo; ++j)
        {
            fscanf(in, "%d", &dof);
            m_node[i][j] = dof;
            // printf("%d ", dof);
        }
        // printf("\n");
    }
    fclose(in);
}

//
// Reads graph from file (BINARY mode)
//
void Graph::ReadBin(const char* path)
{
FILE* in;
int nodeNo, eltNo;
int buf[6];

    in = fopen(path, "rb");
    if(!in)
    {
        printf("Cannot open file '%s' for reading.\n", path);
        exit(1);
    }

    fread(buf, sizeof(int), 6, in);
    m_nx     = buf[0];
    m_ny     = buf[1];
    m_nz     = buf[2];
    m_deg    = buf[3];
    m_ilower = buf[4];
    m_iupper = buf[5];


    // printf("%d %d\n", m_ilower, m_iupper);
    nodeNo = m_iupper - m_ilower + 1;
    m_node.resize(nodeNo);

    for(int i = 0; i < nodeNo; ++i)
    {
        fread(&eltNo, sizeof(int), 1, in);
        m_node[i].resize(4 * eltNo);
        fread(&(m_node[i][0]), sizeof(int), 4 * eltNo, in);
        // printf(".");
    }
    fclose(in);
}

//
// Writes graph into file (TEXT mode)
//
void Graph::WriteTxt(const char* path)
{
FILE *out;

    out = fopen(path, "wt");
    if(!out)
    {
        printf("Cannot open file '%s' for writing.\n", path);
        exit(1);
    }

    fprintf(out, "%d %d %d %d %d\n", m_nx, m_ny, m_nz, m_ilower, m_iupper);

    for(size_t i = 0; i < m_node.size(); ++i)
    {
        fprintf(out, "%d ", static_cast<int>(m_node[i].size() / 4));
        for(size_t j = 0; j < m_node[i].size(); j++)
        {
            fprintf(out, "%d ", m_node[i][j]);
        }
        fprintf(out, "\n");
    }
    fprintf(out, "\n");
    fclose(out);
}

#endif
