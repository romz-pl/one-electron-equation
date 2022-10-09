#pragma once

#include <vector>

class GenDOS
{
public:
	GenDOS(int pntNo, double sigma);
	
	void Read(const char* path);
	
	void Write(const char* path);
	
private:
	double GetVal(double e) const;
	double GetGauss(double x) const;
	
private:

	std::vector<double> m_eig;
	
	const int m_pntNo;
	
	const double m_sigma; 
	
	double m_eMin;
	
	double m_eMax; 
	
	
};
