#include <cassert>
#include <vector>
#include <cmath>
#include "cubat.h"



//
// Constructor
//
Cubat::Cubat(int order) : m_order(order)
{
    assert(order > 0);
    InitGaussProd(order);
}

//
// Initializes cubatures as a tensor product of one-dimensional Gauss quadratures
//
void Cubat::InitGaussProd(int order)
{
const int size = order * order * order;
std::vector<double> w1D(order), x1D(order);
int k;

    GauLeg(-1., 1., x1D, w1D);

    m_node.resize(size);

    k = 0;
    for(int s = 0; s < order; s++)
    {
        for(int t = 0; t < order; t++)
        {
            for(int u = 0; u < order; u++)
            {
                m_node[k] = CubatNode(x1D[s], x1D[t], x1D[u], w1D[s] * w1D[t] * w1D[u]);
                k++;
            }
        }
    }
}

//
// Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
// arrays x[1..n] and w[1..n] of length n, containing the abscissas and weights of the Gauss-
// Legendre n-point quadrature formula.
//
// From "Numerical Recipes in C"
//
void Cubat::GauLeg(double x1, double x2, std::vector<double>& x, std::vector<double>& w)
{
// const double pi = 3.1415926535897932384626433832795;
const double pi = 4 * atan(1.);
const double EPS = 3.0e-15;
int m, j, i, n;
double z1, z, xm, xl, pp, p3, p2, p1; // High precision is a good idea for this routine.

    assert(x.size() == w.size());
    n = static_cast<int>(x.size());

    m = (n + 1) / 2; // The roots are symmetric in the interval, so we only have to find half of them.
    xm = 0.5 * (x2 + x1);
    xl = 0.5 * (x2 - x1);
    for(i = 0; i < m; i++) // Loop over the desired roots.
    {
        z = cos(pi * (i + 1 - 0.25) / (n + 0.5));
        // Starting with the above approximation to the ith root, we enter the main loop of
        // refinement by Newton's method.
        do
        {
            p1 = 1.0;
            p2 = 0.0;
            // Loop up the recurrence relation to get the Legendre polynomial evaluated at z.
            for(j = 1; j <= n; j++)
            {
                p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
            }
            // p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
            // by a standard relation involving also p2, the polynomial of one lower order.
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp; // Newton's method.
        }
        while(fabs(z - z1) > EPS);

        x[i] = xm - xl * z;			// Scale the root to the desired interval,
        x[n - 1 - i] = xm + xl * z; // and put in its symmetric counterpart.
        w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp); // Compute the weight...
        w[n - 1 - i] = w[i];						 // and its symmetric counterpart.
    }
}

