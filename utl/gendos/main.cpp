#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "gendos.h"



int main(int argc, char* argv[])
{

	if(argc != 3)
	{
		printf("gendox.x <filename> <sigma>\n");
		exit(1);
	}
	
	std::string in(argv[1]);
	std::string out = in + ".dos";

	GenDOS gen(10000, atof(argv[2]));

	gen.Read( in.c_str() );
	gen.Write( out.c_str() );

	return 0;
}

