#ifndef RSCHR_POTEN__H
#define RSCHR_POTEN__H

//
// AUTHOR: Zbigniew Romanowski [romz@wp.pl]
//

//
// 1. Interaction potential.
//
// 2. In the current version of program, interaction potential is represented
//    as a sum of "basic potentials" centered at "potential centers".
//
// 3. "Basisc potentials" are:
//        a) Gauss potential, class PotGauss.
//        b) Soft Coulomb potential, class PotSoftCoul.
//        c) Potential of three-dimensional harmonic oscilator, class PotHarm.
//        d) Hydrogen atom with HGH pseudopotential, class PotHGHhydro.
//
// 4. "Basisc potentials" are deriveed form class Fun3D, which represents
//     function from R^3 -> R
//
// 5. In order to extend functionality, new class must be created.
//



#include <vector>
#include <cstdio>
#include <cmath>



//
// Virtual basis class representing potential
//
class Fun3D
{
public:
    Fun3D(void) { }
    virtual ~Fun3D(void) { }
    virtual double Get(double x, double y, double z) const = 0;
    virtual void Info(FILE* out) const = 0;
};

//
// Gauss potential
//
class PotGauss : public Fun3D
{
public:
    PotGauss(double d, double a) : m_d(d), m_a(a) { }
    virtual ~PotGauss(void) { }

    virtual double Get(double x, double y, double z) const
    {
        const double r2 = x*x + y*y + z*z;
        return m_d * exp(m_a * r2);
    }
    virtual void Info(FILE* out) const
    {
        fprintf(out, "    Gauss(d,a)=(%5.2lf,%5.2lf)", m_d, m_a);
    }
private:
    double m_d, m_a;
};

//
// Soft Coulomb potential
//
class PotSoftCoul : public Fun3D
{
public:
    PotSoftCoul(double d, double a) : m_d(d), m_a(a) { }
    virtual ~PotSoftCoul(void) { }

    virtual double Get(double x, double y, double z) const
    {
        const double r2 = x*x + y*y + z*z;
        return m_d / sqrt(r2 + m_a);
    }
    virtual void Info(FILE* out) const
    {
        fprintf(out, "    Soft-Coulomb(d,a)=(%5.2lf,%5.2lf)", m_d, m_a);
    }
private:
    double m_d, m_a;
};


//
// potentiala d/(r^2 + a)
//
class PotA : public Fun3D
{
public:
    PotA(double d, double a) : m_d(d), m_a(a) { }
    virtual ~PotA(void) { }

    virtual double Get(double x, double y, double z) const
    {
        const double r2 = x*x + y*y + z*z;
        return m_d / (r2 + m_a);
    }
    virtual void Info(FILE* out) const
    {
        fprintf(out, "    Soft-Squared(d,a)=(%5.2lf,%5.2lf)", m_d, m_a);
    }
private:
    double m_d, m_a;
};

//
// potentiala d/(r^3 + a)
//
class PotB : public Fun3D
{
public:
    PotB(double d, double a) : m_d(d), m_a(a) { }
    virtual ~PotB(void) { }

    virtual double Get(double x, double y, double z) const
    {
        const double r2 = x*x + y*y + z*z;
        const double r3 = r2 * sqrt(r2);
        return m_d / (r3 + m_a);
    }
    virtual void Info(FILE* out) const
    {
        fprintf(out, "    Soft-Cubed(d,a)=(%5.2lf,%5.2lf)", m_d, m_a);
    }
private:
    double m_d, m_a;
};


//
// potentiala d/(r^4 + a)
//
class PotD : public Fun3D
{
public:
    PotD(double d, double a) : m_d(d), m_a(a) { }
    virtual ~PotD(void) { }

    virtual double Get(double x, double y, double z) const
    {
        const double r2 = x*x + y*y + z*z;
        const double r4 = r2 * r2;
        return m_d / (r4 + m_a);
    }
    virtual void Info(FILE* out) const
    {
        fprintf(out, "    Soft-r^4(d,a)=(%5.2lf,%5.2lf)", m_d, m_a);
    }
private:
    double m_d, m_a;
};

//
// potentiala s * exp(a*r^2)
//
class PotE : public Fun3D
{
public:
    PotE(double d, double a) : m_d(d), m_a(a) { }
    virtual ~PotE(void) { }

    virtual double Get(double x, double y, double z) const
    {
        const double r2 = x*x + y*y + z*z;
        return m_d * exp(m_a * r2);
    }
    virtual void Info(FILE* out) const
    {
        fprintf(out, "    Gauss(d,a)=(%5.2lf,%5.2lf)", m_d, m_a);
    }
private:
    double m_d, m_a;
};

//
// potentiala d/(r^5 + a)
//
class PotF : public Fun3D
{
public:
    PotF(double d, double a) : m_d(d), m_a(a) { }
    virtual ~PotF(void) { }

    virtual double Get(double x, double y, double z) const
    {
        const double r2 = x*x + y*y + z*z;
        const double r5 = r2 * r2 * sqrt(r2);
        return m_d / (r5 + m_a);
    }
    virtual void Info(FILE* out) const
    {
        fprintf(out, "    Soft-r^5(d,a)=(%5.2lf,%5.2lf)", m_d, m_a);
    }
private:
    double m_d, m_a;
};


//
//
// Potential of three-dimensional harmonic oscilator
//
class PotHarm : public Fun3D
{
public:
    PotHarm(void) { }
    virtual ~PotHarm(void) { }

    virtual double Get(double x, double y, double z) const
    {
        const double k = 0.5;
        return k * (x*x + y*y + z*z);
    }
    virtual void Info(FILE* out) const
    {
        fprintf(out, "%s", "    Harmonic-oscylator()");
    }
};


//
// Hydrogen atom with HGH pseudopotential
//
class PotHGHhydro : public Fun3D
{
public:
    PotHGHhydro(void) { }
    virtual ~PotHGHhydro(void) { }

    virtual double Get(double x, double y, double z) const
    {
        // Parameters from:
        //     Phys. Rev B, vol. 58, pp. 3641 (1998)
        const double zion = 1, rloc = 0.2, c1 = -4.180237, c2 = 0.725075, c3 = 0, c4 = 0;
        static const double sqrt2 = sqrt(2.);

        const double r = ::sqrt(x*x + y*y + z*z);
        const double w = r / rloc;
        const double w2 = w * w;
        const double w4 = w2 * w2;
        const double w6 = w4 * w2;
        const double pol = c1 + c2 * w2 + c3 * w4 + c4 * w6;
        const double tt = ::erf(w / sqrt2) / r;

        // return -(zion / r) * ::erf(w / sqrt2) + exp(-0.5 * w2) * pol;
        return -zion * tt + exp(-0.5 * w2) * pol;
    }
    virtual void Info(FILE* out) const
    {
        fprintf(out, "%s", "    HGH-pp-hydrogen()");
    }
};

//
// Potential centre
//
class PotCntr
{
    friend class Poten;
public:
    PotCntr(void) : m_x(0.), m_y(0), m_z(0.) { }
    PotCntr(double x, double y, double z) : m_x(x), m_y(y), m_z(z) { }
private:
    double m_x, m_y, m_z;
};



//
// Interaction potential as a sum of "Basic potentials".
//
class Poten
{
public:
    Poten( ) = default;
    ~Poten( ) = default;

    double Get(double x, double y, double z) const;

    int Read(const char* buf, int len);
    void Info(FILE* out) const;


private:
    void Clear();
    int ReadLine( const char* line );

private:
    // Array of potential centers
    std::vector< PotCntr > m_cntr;

    // Array of "basic potentials"
    std::vector< Fun3D* > m_fun;
};

#endif

