#include <cassert>
#include "bspline.h"


//
// Returns value of one dimensional B-spline of degree "n" at point "x".
// Evaluated recursively.
//
double BSpline::Value(double x, int n)
{

    assert(n >= 0);
    if(n == 0)
    {
        if(0. <= x && x < 1.)
            return 1;
        else
            return 0;

        // return static_cast<double>(0. <= x && x <= 1.);
    }

    const double a = x / n;
    const double b = (n + 1 - x) / n;

    return a * Value(x, n - 1) + b * Value(x - 1, n - 1);
}

//
// Returns value of first derivative of one dimensional B-spline of degree "n" at point "x"
//
double BSpline::ValueD(double x, int n)
{
    return Value(x, n - 1) - Value(x - 1, n - 1);
}

//
// Helper function used for evaluation of overlap and
// kinetic integrals analytically
//
double BSpline::S(int k, int n)
{
    return BSpline::Value(n + 1 + k, 2 * n + 1);
}


//
//-----------------------------------------------------------------------------
//

//
// Returns value of one dimensional B-spline of degree "n" at point "x"
//
// IT DOES NOT WORK!
//
double BSpline::ValueFast(double x, int n)
{
double v;
    assert(n >= 0);

//	if(x < 0. || x > n + 1)
//		return 0.;

    switch(n)
    {
        case 0: v = Core0(x); break;
        case 1: v = Core1(x); break;
        case 2: v = Core2(x); break;
        case 3: v = Core3(x); break;
        case 4: v = Core4(x); break;
        case 5: v = Core5(x); break;

        default: v = Value(x, n); break;
    }
    return v;
}

double BSpline::Core0(double x)
{
//	if(0. <= x && x < 1.)
//		return 1;
//
//	return 0;

    return static_cast<double>(0. <= x && x < 1.);
}

double BSpline::Core1(double x)
{
//	if(0. <= x && x < 1.)
//		return x;
//
//	if(1. <= x && x < 2.)
//		return 2 - x;
//
//	return 0;

    // This avoids "if" statement
    const bool b1 = (0. <= x && x < 1.);
    const bool b2 = (1. <= x && x < 2.);

    return b1 * x + b2 * (2 - x);
}


double BSpline::Core2(double x)
{
static const double c[3][3] = {
        {   0.,  0., 0.5},
        { -1.5,  3., -1.},
        {9./2., -3., 0.5}
    };
const int p = static_cast<int>(x);

    assert(x >= 0. && x <= 3.);
    return c[p][0] + x * (c[p][1] + x * c[p][2]);
}


/*
double BSpline::Core3(double x)
{
    assert(x >= 0.);
    assert(x <= 4.);

    if(0. <= x && x < 1.)
        return x * x * x / 6.;

    if(1. <= x && x < 2.)
    {
        static const double c[4] = {2./3., -2, 2, -0.5};
        return c[0] + x*c[1] + x*x*c[2] + x*x*x*c[3];
    }

    if(2. <= x && x < 3.)
    {
        static const double c[4] = {-22./3., 10, -4, 0.5};
        return c[0] + x*c[1] + x*x*c[2] + x*x*x*c[3];
    }

    if(3. <= x && x < 4.)
    {
        x -= 4;
        return -x*x*x / 6.;
    }

    return 0;
}
*/

//
// Returns B-spline of degree 3 at "x"
//
double BSpline::Core3(double x)
{
static const double c[4][4] = {
        {     0.,   0.,   0.,  1./6.},
        {  2./3.,  -2.,   2.,   -0.5},
        {-22./3.,  10.,  -4.,    0.5},
        { 32./3.,  -8.,   2., -1./6.}
    };
const int p = static_cast<int>(x);
const double *r = c[p];

    assert(x >= 0. && x <= 4.);
    return r[0] + x * (r[1] + x * (r[2] + x * r[3]));
}

//
// Returns B-spline of degree 4 at "x"
//
double BSpline::Core4(double x)
{
static const double c[5][5] = {
        {       0.,       0.,        0.,      0.,  1./24.},
        {  -5./24.,  20./24.,  -30./24.,  20./24, -4./24.},
        { 155./24., -300./24., 210./24., -60./24., 6./24.},
        {-655./24., 780./24., -330./24.,  60./24, -4./24.},
        { 625./24., -500./24., 150./24., -20./24., 1./24.}
    };
const int p = static_cast<int>(x);
const double *r = c[p];

    assert(x >= 0. && x <= 5.);
    return r[0] + x * (r[1] + x * (r[2] + x * (r[3] + x * r[4])));
}

//
// Returns B-spline of degree 5 at "x"
//
double BSpline::Core5(double x)
{
static const double c[6][6] = {
        {        0.,          0.,         0.,         0.,       0.,  1./120.},
        {   6./120.,   -30./120.,   60./120.,  -60./120., 30./120., -5./120.},
        { -237./60.,    585./60.,  -570./60.,   270./60., -60./60.,   5./60.},
        { 2193./60.,  -3465./60.,  2130./60.,  -630./60.,  90./60.,  -5./60.},
        {-1829./20.,     409./4.,    -89./2.,     19./2.,      -1.,   1./24.},
        {7776./120., -6480./120., 2160./120., -360./120., 30./120., -1./120.}
    };
const int p = static_cast<int>(x);
const double *r = c[p];

    assert(x >= 0. && x <= 6.);
    return r[0] + x * (r[1] + x * (r[2] + x * (r[3] + x * (r[4] + x * r[5]))));
}



