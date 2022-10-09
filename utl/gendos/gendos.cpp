#include "gendos.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


GenDOS::GenDOS(int pntNo, double sigma) :
	m_pntNo(pntNo), m_sigma(sigma)
{
}

void GenDOS::Read(const char* path)
{
FILE* in;
double v, err;
int n;

	in = fopen(path, "rt");
	if(!in)
	{
		printf("Cannot open file '%s'\n", path);
		exit(1);
	}
	
	m_eig.reserve(10000);
	while(1)
	{
		if( fscanf(in, "%d %lf %lf\n", &n, &v, &err) != 3)
			break;
			
		// printf("%E\n", v);
		m_eig.push_back(v);
	}
	
	m_eMin = m_eig.front();
	m_eMax = m_eig.back();
	
	fclose(in);
}

void GenDOS::Write(const char* path)
{
FILE* out;
double e, v;
double dE = (m_eMax - m_eMin) / m_pntNo;

	out = fopen(path, "wt");
	if(!out)
	{
		printf("Cannot open file '%s'\n", path);
		exit(1);
	}

	for(int i = 0; i <= m_pntNo; i++)
	{
		e = m_eMin + i * dE;
		v = GetVal(e);
		
		fprintf(out, "%E %E\n", e, v);
	}
	fclose(out);

}


double GenDOS::GetVal(double e) const
{
double ret = 0;

	for(size_t i = 0; i < m_eig.size(); i++)
		ret += GetGauss(e - m_eig[i]);
		
	return ret;
}

double GenDOS::GetGauss(double x) const
{
const double s = 2.0 * m_sigma * m_sigma;
const double w = x * x / s;

	return ( exp(-w) / sqrt(M_PI * s) );
	
}
