#include <algorithm>
#include <cassert>
#include <string>
#include "poten.h"
#include "util.h"



//
// Reads potential definition from file "path"
//
int Poten::Read(const char* buf, int len)
{
std::string line;
int i = -1;

    Clear();
    while(true)
    {
        line.clear();
        while(1)
        {
            i++;
            if(i >= len)
                return 0;

            const char c = buf[i];
            if(c == '\n' || c == '\r')
                break;

            line += c;
        }

        // Skip empty lines
        if( line.empty() )
            continue;

        // Skip line with only whitespaces
        if( std::all_of( line.begin(), line.end(), isspace ) )
            continue;

        // Skip comments
        if( line[0] == '#' )
            continue;


        if( ReadLine(line.c_str()) != 0)
            return 1;

    }

    assert(m_cntr.size() == m_fun.size());

    return 0;
}

int Poten::ReadLine(const char* line)
{
double x, y, z;
char flag = ' ';
Fun3D *f;

    sscanf(line, "%c", &flag);
    // printf("flag = %c\n", flag);

    switch(flag)
    {
        case 'G':
        {
            double d, a;
            if(sscanf(line, "%c %lf %lf %lf %lf %lf", &flag, &x, &y, &z, &d, &a) != 6) goto LABELERR;
            f = new PotGauss(d, a);
        }
        break;


        case 'C':
        {
            double d, a;
            if(sscanf(line, "%c %lf %lf %lf %lf %lf", &flag, &x, &y, &z, &d, &a) != 6) goto LABELERR;
            f = new PotSoftCoul(d, a);
        }
        break;


        case 'A':
        {
            double d, a;
            if(sscanf(line, "%c %lf %lf %lf %lf %lf", &flag, &x, &y, &z, &d, &a) != 6) goto LABELERR;
            f = new PotA(d, a);
        }
        break;


        case 'B':
        {
            double d, a;
            if(sscanf(line, "%c %lf %lf %lf %lf %lf", &flag, &x, &y, &z, &d, &a) != 6) goto LABELERR;
            f = new PotB(d, a);
        }
        break;


        case 'D':
        {
            double d, a;
            if(sscanf(line, "%c %lf %lf %lf %lf %lf", &flag, &x, &y, &z, &d, &a) != 6) goto LABELERR;
            f = new PotD(d, a);
        }
        break;


        case 'E':
        {
            double d, a;
            if(sscanf(line, "%c %lf %lf %lf %lf %lf", &flag, &x, &y, &z, &d, &a) != 6) goto LABELERR;
            f = new PotE(d, a);
        }
        break;


        case 'F':
        {
            double d, a;
            if(sscanf(line, "%c %lf %lf %lf %lf %lf", &flag, &x, &y, &z, &d, &a) != 6) goto LABELERR;
            f = new PotF(d, a);
        }
        break;


        case 'H':
        {
            if(sscanf(line, "%c %lf %lf %lf", &flag, &x, &y, &z) != 4) goto LABELERR;
            f = new PotHarm();
        }
        break;

        case 'P':
        {
            if(sscanf(line, "%c %lf %lf %lf", &flag, &x, &y, &z) != 4) goto LABELERR;
            f = new PotHGHhydro();
        }
        break;

        default:
            printf("ERROR-INPUT-PARAMETER: Bad potential");
            goto LABELERR;
        break;
    }

    // printf("flag = %c  x = %lf   y = %lf  z = %lf\n", flag, x, y, z);
    m_cntr.push_back(PotCntr(x, y, z));
    m_fun.push_back(f);

    return 0;


LABELERR:
    printf("POT.DEF: Invalid format at line '%s'.", line);
    return 1;
}

//
// Frees memory
//
void Poten::Clear( )
{
    for(size_t i = 0; i < m_fun.size(); i++)
        delete m_fun[i];

    m_fun.clear();
    m_cntr.clear();
}


//
// Writes info about interaction potential
//
void Poten::Info(FILE* out) const
{

    fprintf(out, "NUMBER OF 'ATOMS' = %d\n", static_cast<int>(m_fun.size()));
    for(size_t i = 0; i < m_fun.size(); i++)
    {
        m_fun[i]->Info(out);
        fprintf(out, " at (x,y,z) = (%5.2lf,%5.2lf,%5.2lf)\n", m_cntr[i].m_x, m_cntr[i].m_y, m_cntr[i].m_z);
    }
    printf("\n\n");
    fflush(out);

}


//
// Returns value of potential at (x,y,z)
//
double Poten::Get(double x, double y, double z) const
{
double val = 0.;

    assert(m_cntr.size() == m_fun.size());

    for(size_t i = 0; i < m_fun.size(); i++)
    {

        const PotCntr& c = m_cntr[i];
        // shifted center (xs, ys, zs)
        const double xs = x - c.m_x;
        const double ys = y - c.m_y;
        const double zs = z - c.m_z;

        val += m_fun[i]->Get(xs, ys, zs);

    }
    return val;
}

