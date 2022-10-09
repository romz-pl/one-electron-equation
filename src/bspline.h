#ifndef RSCHR_BSPLINE__H
#define RSCHR_BSPLINE__H

//
// AUTHOR: Zbigniew Romanowski [romz@wp.pl]
//

//
// Three dimmensional B-spline.
//
// This class evaluates:
//    1. Values of B-spline of arbitrary degree.
//    2. Values of the first derivative of B-spline of arbitrary degree.
//    3. Values needed for evaluation of overlap and kinetic integrals
//


class BSpline
{
public:
    BSpline(void);
    ~BSpline(void);

    static double S(int k, int n);

    static double Value(double x, int n);
    static double ValueD(double x, int n);

private:
    static double ValueFast(double x, int n);



    static double Core0(double x);
    static double Core1(double x);
    static double Core2(double x);
    static double Core3(double x);
    static double Core4(double x);
    static double Core5(double x);
};

#endif

