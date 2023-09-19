/*
	gcc -O 00_base.c ofd_datalib.c -o 00
	./00
*/


#include "ofd_datalib.h"
#include <math.h>
#include <stdio.h>

int main(void)
{
	const int nxarray = 5;
	const int nyarray = 50;
	const int ns = 5;
	const int nt = 5;
	const double freq = 14e9;
	const double lambda = 3e8 / freq;
	const double dx = 0.1*lambda;
	const double dy = 0.1*lambda;
	const double er = 2.3;
	const double tand = 0.001;
	const double dt = 1.6e-3;
	const double h = 20e-3;
	const double z_div = 16e-3;
	const double dh = 16e-3;
	const double margin = 0.25*lambda;
	const double l[] = {2.13e-3, 2.04e-3, 1.91e-3, 1.75e-3, 1.57e-3};
	const double m[] = {2e-3, 1.5e-3, 1e-3, 0.5e-3, 0.2e-3};

	const double eps0 = 8.854e-12;
	const double pi = 4 * atan(1);

	// initialize

	ofd_init();

	// title

	ofd_title("00_base");

	// mesh

	const double x0 = -nxarray * ns * dx;
	const double x1 = -x0;
	const int xmesh_margin = (margin) / (0.01*lambda);
	const int xmesh_base = (x1 - x0) / (0.005*lambda);
	ofd_xsection(4, x0 - margin, x0, x1, x1 + margin);
	ofd_xdivision(3, xmesh_margin, xmesh_base, xmesh_margin);

	const double y0 = -nyarray / 2.0 * dy;
	const double y1 = -y0;
	const int ymesh_margin =  (margin) / (0.01*lambda);
	const int ymesh_base = (y1 - y0) / (0.005*lambda);
	ofd_ysection(4, y0 - margin, y0, y1, y1 + margin);
	ofd_ydivision(3, ymesh_margin, ymesh_base, ymesh_margin);

	const double z0 = 0;
	const double z1 = z0 + h;
	const double z2 = z1 + dt;
	const double z3 = z2 + 0.5*lambda;
	const int zmesh_d = (0.5*lambda) / (0.005*lambda);
	const int zmesh_base = (dt) / (0.005*lambda);
	const int zmesh_h = h / (0.005*lambda);
	const int zmesh_margin = (margin) / (0.01*lambda);
	ofd_zsection(5, z0 - margin, z0, z1, z2, z3);
	ofd_zdivision(4, zmesh_margin, zmesh_h, zmesh_base, zmesh_d);

	// material

	const double sigma = er * tand * (2 * pi * freq) * eps0;
	ofd_material(er, sigma, 1, 0, "");

	// geometry

	// substrate
	ofd_geometry(2, 1, x0, x1, y0, y1, z1, z2);


	for (int i = 0; i < nxarray; i++) {
	for (int j = 0; j < nyarray; j++) {
		const double x = x0 + (2*i * ns * dx);
		const double y = y0 + (j + 0.5) * dy;
		const double xg = x0 + (2*i * ns * dx) + (ns * dx);
		// reconfigurable stracture
		ofd_geometry(2, 1, x, x + (ns * dx), y0, y1, z0, z_div);
		ofd_geometry(2, 1, xg, xg + (ns * dx), y0, y1, z_div, z1);
		// metalic wall
		ofd_geometry(1, 1, x, x, y0, y1, z0, z1);
		ofd_geometry(1, 1, xg, xg, y0, y1, z0, z1);
		ofd_geometry(1, 1, x1, x1, y0, y1, z0, z1);

		for (int p = 0; p < ns; p++) {
			const double xa = x + (p + 0.5) * dx;
			ofd_geometry(1, 1, xa - l[p] /2, xa + l[p]/2, y -l[p]/2, y +l[p]/2, z2, z2);
		}
		for (int g = 0; g < nt; g++) {
			const double xb = xg + (g + 0.5) * dx;
			//ofd_geometry(1, 1, xb, x1, -y/2, +y/2, z2, z2);
			ofd_geometry(0, 1, xb - m[g] /2, xb + m[g]/2, y -m[g]/2, y +m[g]/2, z2, z2);
		}
	}
	}

	// ground
	ofd_geometry(1, 1, x0, x1, y0, y1, z0, z0);

	// plane wave

	ofd_planewave(0, 90, 1);

	// ABC

	ofd_pml(5, 2, 1e-5);

	// frequency

	ofd_frequency1(freq, freq, 0);
	ofd_frequency2(freq, freq, 0);

	// solver

	ofd_solver(30000, 100, 1e-3);

	// iteration

	ofd_plotiter(1);

	// far-1d

	ofd_plotfar1d('Y', 360, 0);

	// output options

	ofd_far1dstyle(0);
	ofd_far1dcomponent(1, 0, 0);
	ofd_far1ddb(1);

	// near-2d
	ofd_plotnear2d("Ey",'Y',0);


	// output option
	ofd_near2ddim(1,0);
	// output

	ofd_outdata("00_base.ofd");

	return 0;
}