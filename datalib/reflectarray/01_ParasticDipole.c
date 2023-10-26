/*
	gcc -O 01_ParasticDipole.c ofd_datalib.c -o 01
	./01
*/


#include "ofd_datalib.h"
#include <math.h>
#include <stdio.h>

int main(void)
{
	const int nxarray = 1;
	const int nyarray = 1;
	const int ns = 4;           // Count of Upper cell
	const int nt = 4;           // Count of Lower cell
	const double freq = 28e9;
	const double lambda = 3e8 / freq;
	const double dx = 0.5*lambda;   // UnitCell size
	const double dy = 0.5*lambda;   // UnitCell size
	const double er = 2.5;
	const double tand = 0.0015;
	const double dt = 1.6e-3;   // Thicness of BaseDielectric
	const double h = 10e-3;     // Hight of Space of InsertDielectric 
	const double z_div = 5e-3;// Diverce of InsertDielectric
	const double margin = 1*lambda;
	const double l[] = {0.67*dx, 0.52*dx, 0.43*dx, 0.31*dx};
	const double m[] = {0.58*dx, 0.55*dx, 0.52*dx, 0.43*dx};
    const double r = 0.7;       // ratio of lm/lp
    const double wm = 0.2*dx;   // Width of main Dipole
    const double wp = 0.1*dx;   // Width od main Dipole
    const double dist = 0.5e-3; // Distance between MainDipole and ParasticDipole 

	const double eps0 = 8.854e-12;
	const double pi = 4 * atan(1);

	// initialize

	ofd_init();

	// title

	ofd_title("01_ParasticDipole");

	// mesh

	const double x0 = -nxarray * (ns+nt) /2 * dx;
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
	// ofd_geometry(2, 1, x0, x1, y0, y1, z1, z2);


	for (int i = 0; i < nxarray; i++) {
	for (int j = 0; j < nyarray; j++) {
		const double x = x0 + (i * (ns+nt) * dx);
		const double y = y0 + (j + 0.5) * dy;
		const double xg = x0 + (i * (ns+nt) * dx) + (ns * dx);

		// Upper and Lower Dielectric
		ofd_geometry(2, 1, xg, xg + (nt * dx), y0, y1, z0, z_div);
		ofd_geometry(2, 1, x, x + (ns * dx), y0, y1, z_div, z1);

		for (int p = 0; p < ns; p++) {  // Upper Type
			const double xa = x + (p + 0.5) * dx;   // Center of UC
            // Main Dipole
			ofd_geometry(1, 1, xa - l[p]/2, xa + l[p]/2, y -wm/2, y +wm/2, z2, z2);
			// Parastic Dipole
            ofd_geometry(1, 1, xa - l[p]*r/2, xa + l[p]*r/2, y-(wm/2)-dist-wp, y-(wm/2)-dist, z2, z2);
            ofd_geometry(1, 1, xa - l[p]*r/2, xa + l[p]*r/2, y+(wm/2)+dist, y+(wm/2)+dist+wp, z2, z2);
            // metalic wall
			ofd_geometry(1, 1, x + p * dx, x + p * dx, y0, y1, z0, z1);
		}
		for (int g = 0; g < nt; g++) {  // Lower Type
			const double xb = xg + (g + 0.5) * dx;  // Center of UC
			// Main Dipole
            ofd_geometry(1, 1, xb - m[g]/2, xb + m[g]/2, y -wm/2, y +wm/2, z2, z2);
            // Parastic Dipole
            ofd_geometry(1, 1, xb - m[g]*r/2, xb + m[g]*r/2, y-(wm/2)-dist-wp, y-(wm/2)-dist, z2, z2);
			ofd_geometry(1, 1, xb - m[g]*r/2, xb + m[g]*r/2, y+(wm/2)+dist, y+(wm/2)+dist+wp, z2, z2);
            // metalic wall
			ofd_geometry(1, 1, xg + g * dx, xg + g * dx, y0, y1, z0, z1);
			ofd_geometry(1, 1, x1, x1, y0, y1, z0, z1);
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

	ofd_outdata("01_ParasticDipole.ofd");

	return 0;
}