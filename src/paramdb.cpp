#include <cassert>
#include <algorithm>
#include <iostream>
#include <sstream>
#include "paramdb.h"
#include "util.h"


//
// Returns value of parameter "param".
//
std::string ParamDb::GetStr( const std::string& param ) const
{
    std::map < std::string, std::string >::const_iterator iter = find( param );

    if( iter != end() )
        return iter->second;


    const std::string buf = "Parameter not found in input file. Parm = " + param;
    Throw_invalid_argument( buf );
    return std::string();
}

//
// Returns value of parameter "param". Parameter must be of type "double".
//
double ParamDb::GetDbl( const std::string& param ) const
{
    return std::stod( GetStr( param ) );
}

//
// Returns value of parameter "param". Parameter must be of type "int".
//
int ParamDb::GetInt( const std::string& param ) const
{
    return std::stoi( GetStr( param ) );
}

//
// Returns value of parameter "param". Parameter must be of type "char".
//
char ParamDb::GetChr( const std::string& param ) const
{
    return GetStr( param )[0];
}

//
// Returns "true" if the paramater "param" is equal to "Y"
//
bool ParamDb::GetBool( const std::string& param ) const
{
    return ( GetStr( param ) == "Y" );
}



//
// Inserts parameter "param" with vale "val"
//
void ParamDb::Insert( const std::string& param, const std::string& val )
{
    if( find( param ) != end() )
    {
        const std::string buf = "Duplication of parameter in input file. Param = " + param;
        Throw_invalid_argument( buf );
    }

    insert( std::pair < std::string, std::string >( param, val ) );
}


//
// Reads parameters from the buffer "buf".
// The length of the buffer is "len".
//
void ParamDb::Read( const char* buf, int len )
{
int i = -1;
std::string line, param, val;


    while(1)
    {
        line.clear();

        while(1)
        {
            i++;
            if(i >= len)
                return;

            const char c = buf[i];
            if(c == '\n' || c == '\r')
                break;

            line += c;

        }
        // printf("%s\n", line.c_str());


        // Skip empty lines
        if( line.empty() )
            continue;

        // Skip line with only whitespaces
        if( std::all_of( line.begin(), line.end(), isspace ) )
            continue;

        // Skip comments
        if( line[0] == '#' )
            continue;

        param.clear();
        val.clear();
        std::stringstream( line ) >> param >> val;

        if( param.empty() || val.empty() )
        {
            const std::string buf = "Error during reading input file. Line = " +  line;
            Throw_invalid_argument( buf );
        }

        Insert( param, val );
    }

}


