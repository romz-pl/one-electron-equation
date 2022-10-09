#include <stdio.h>


int main(void)
{
const int Nx = 8, Ny = 8, Nz = 8;
const double dx = 1.7, dy = 1.7, dz = 1.7;
const double a = -2.0, d = -20.0;
int ix, iy, iz;
double x, y, z;

	x = 0;
	for(ix = 0; ix < Nx; ix++)
	{
		y = 0;
		for(iy = 0; iy < Ny; iy++)
		{
			z = 0;
			for(iz = 0; iz < Nz; iz++)
			{
				printf("E  %lf  %lf  %lf  %lf  %lf\n", x, y, z, d, a);
				z += dz;
			}
			printf("#\n");
			y += dy;
		}
		x += dx;
	}
}
