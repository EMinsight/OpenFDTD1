#include "ofd_datalib.h"
#include <math.h>
#include <stdio.h>

int main(void)
{
	const int nxarray = 5;
	const int nyarray = 25;
	const int ns = 5;
	const double freq = 14e9;
	const double lambda = 3e8 / freq;
	const double dx = 0.1*lambda;
	const double dy = 0.1*lambda;
	const double er = 2.3;
	const double tand = 0.001;
	const double dt = 1.6e-3;
	const double h = 20e-3;
	const double dh = 16e-3;
	const double margin = 0.25*lambda;
	const double l[] = {2.13e-3, 2.04e-3, 1.91e-3, 1.75e-3, 1.57e-3};

	const double eps0 = 8.854e-12;
	const double pi = 4 * atan(1);

	// initialize

	ofd_init();

	// title

	ofd_title("00_base");

	// mesh

	const double x0 = -nxarray / 2.0 * ns * dx;
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

	// ground
	ofd_geometry(1, 1, x0, x1, y0, y1, z0, z0);

	// substrate
	ofd_geometry(2, 1, x0, x1, y0, y1, z1, z2);
	ofd_geometry(2, 1, x0, x1, y0, y1, z1 - dh, z1);

	for (int i = 0; i < nxarray; i++) {
	for (int j = 0; j < nyarray; j++) {
		const double x = x0 + (i * ns * dx);
		const double y = y0 + (j + 0.5) * dy;
		for (int n = 0; n < ns; n++) {
			const double xa = x + (n + 0.5) * dx;
			ofd_geometry(1, 1, xa - l[n] /2, xa + l[n]/2, y -l[n]/2, y +l[n]/2, z2, z2);
		}
	}
	}

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