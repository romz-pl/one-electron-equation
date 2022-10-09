#ifndef RSCHR_PARAMDB__H
#define RSCHR_PARAMDB__H

//
// AUTHOR: Zbigniew Romanowski [romz@wp.pl]
//

//
// 1. In-memory database of input parameters.
//
// 2. The databe is represented as hash array of pait (key, value),
//    where both "key" and "value" are strings.
//
// 3. Conversion of values to proper type is done by access functions.
//
// 4. The database is represented as singleton.
//
// 5. For convenience auxiliary functions are provided, which are recomended to use.
//

#include <map>
#include <string>

class ParamDb : private std::map< std::string, std::string >
{
public:
    ParamDb( ) = default;
    ~ParamDb( ) = default;


public:
    void Read( const char* buf, int len );


    std::string GetStr ( const std::string& param ) const;
    double      GetDbl ( const std::string& param ) const;
    int         GetInt ( const std::string& param ) const;
    char        GetChr ( const std::string& param ) const;
    bool        GetBool( const std::string& param ) const;


private:
    void Insert( const std::string& param, const std::string& val );

};

#endif

