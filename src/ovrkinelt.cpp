#include "ovrkinelt.h"
#include "bspline.h"

OvrKinElt::OvrKinElt(const ParamDb& db)
{
    Precalc(db);
}

//
// Calculates one-dimensional overlap and kinetic orbitals ANALITYCALY
//
void OvrKinElt::Precalc(const ParamDb& db)
{
const int deg = db.GetInt("General_BsplDeg");
const double h = db.GetDbl("Domain_h");

double a, b, c;

    m_ovr.resize(deg + 1);
    m_der.resize(deg + 1);

    for(int k = 0; k < static_cast<int>(m_ovr.size()); k++) // For each overlap integral
        m_ovr[k] = h * BSpline::S(k, deg);


    for(int k = 0; k < static_cast<int>(m_der.size()); k++) // For each derivative integral
    {
        a = BSpline::S(k, deg - 1);
        b = BSpline::S(k - 1, deg - 1);
        c = BSpline::S(k + 1, deg - 1);
        m_der[k] = (2 * a - b - c) / h;
    }

//	for(int k = 0; k < m_der.size(); k++)
//		fprintf(logGlob, "Ovr[%d] = %lf    Der[%d] = %lf\n", k, m_ovr[k], k, m_der[k]);
}

