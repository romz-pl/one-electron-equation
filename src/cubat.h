#ifndef RSCHR_CUBAT__H
#define RSCHR_CUBAT__H


//
// AUTHOR: Zbigniew Romanowski [romz@wp.pl]
//

//
// Cubature: quadrature over cube [-1,1] x [-1,1] x [-1,1]
//
// This is a tensor product of one-dimensional Gauss quadratures.
//


#include <vector>

class CubatNode
{
    friend class Cubat;
public:
    CubatNode() : m_r(0), m_s(0), m_t(0), m_w(0) { }
    CubatNode(double r, double s, double t, double w) : m_r(r), m_s(s), m_t(t), m_w(w) { }
    ~CubatNode() = default;

private:
    // Coordinates (r, s, t) of node
    double m_r, m_s, m_t;

    // Weight of node
    double m_w;

};


class Cubat
{
public:
    explicit Cubat(int order);
    ~Cubat() = default;

    int Q() const { return m_node.size(); }

    double R(int i) const { return m_node[i].m_r; }
    double S(int i) const { return m_node[i].m_s; }
    double T(int i) const { return m_node[i].m_t; }
    double W(int i) const { return m_node[i].m_w; }

    int Order() const { return m_order; }

private:
    void InitGaussProd(int order);

    static void GauLeg(double x1, double x2, std::vector<double>& x, std::vector<double>& w);

private:
    // Array of cubature nodes
    std::vector<CubatNode> m_node;

    // Order of cubature
    const int m_order;
};

#endif

