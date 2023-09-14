#ifndef _OFD_DATALIB_H_
#define _OFD_DATALIB_H_
extern void ofd_init(void);
extern void ofd_section_size(int);
extern void ofd_material_size(int);
extern void ofd_geometry_size(int);
extern void ofd_feed_size(int);
extern void ofd_load_size(int);
extern void ofd_point_size(int);
extern void ofd_far1d_size(int);
extern void ofd_near1d_size(int);
extern void ofd_near2d_size(int);
extern void ofd_title(const char []);
extern void ofd_xsection1(double);
extern void ofd_ysection1(double);
extern void ofd_zsection1(double);
extern void ofd_xdivision1(int);
extern void ofd_ydivision1(int);
extern void ofd_zdivision1(int);
extern void ofd_xsection(int, ...);
extern void ofd_ysection(int, ...);
extern void ofd_zsection(int, ...);
extern void ofd_xdivision(int, ...);
extern void ofd_ydivision(int, ...);
extern void ofd_zdivision(int, ...);
extern void ofd_material(double, double, double, double, const char []);
extern void ofd_material_dispersion(double, double, double, double, const char []);
extern void ofd_geometry(int, int, double, double, double, double, double, double);
extern void ofd_geometry_array(int, int, const double []);
extern void ofd_geometry_pillar(int, int, const double []);
extern void ofd_geometry_name(const char []);
extern void ofd_feed(char, double, double, double, double, double, double);
extern void ofd_planewave(double, double, int);
extern void ofd_load(char, double, double, double, char, double);
extern void ofd_point(char, double, double, double, const char []);
extern void ofd_rfeed(double);
extern void ofd_pulsewidth(double);
extern void ofd_timestep(double);
extern void ofd_frequency1(double, double, int);
extern void ofd_frequency2(double, double, int);
extern void ofd_solver(int, int, double);
extern void ofd_pml(int, double, double);
extern void ofd_pbc(int, int, int);
extern void ofd_matchingloss(void);
extern void ofd_plotiter(int);
extern void ofd_plotfeed(int);
extern void ofd_plotpoint(int);
extern void ofd_plotsmith(void);
extern void ofd_plotzin(int, double, double, int);
extern void ofd_plotyin(int, double, double, int);
extern void ofd_plotref(int, double, double, int);
extern void ofd_plotspara(int, double, double, int);
extern void ofd_plotcoupling(int, double, double, int);
extern void ofd_plotfar0d(double, double, int, double, double, int);
extern void ofd_freqdiv(int);
extern void ofd_plotfar1d(char, int, double);
extern void ofd_far1dcomponent(int, int, int);
extern void ofd_far1dstyle(int);
extern void ofd_far1ddb(int);
extern void ofd_far1dnorm(void);
extern void ofd_far1dscale(double, double, int);
extern void ofd_plotfar2d(int, int);
extern void ofd_far2dcomponent(int, int, int, int, int, int, int);
extern void ofd_far2ddb(int);
extern void ofd_far2dscale(double, double);
extern void ofd_far2dobj(double);
extern void ofd_plotnear1d(const char [], char, double, double);
extern void ofd_near1ddb(int);
extern void ofd_near1dscale(double, double, int);
extern void ofd_near1dnoinc(void);
extern void ofd_plotnear2d(const char [], char, double);
extern void ofd_near2ddim(int, int);
extern void ofd_near2dframe(int);
extern void ofd_near2ddb(int);
extern void ofd_near2dscale(double, double);
extern void ofd_near2dcontour(int);
extern void ofd_near2dobj(int);
extern void ofd_near2dnoinc(void);
extern void ofd_near2dzoom(double, double, double, double);
extern void ofd_window2d(int, int, int, int);
extern void ofd_window3d(int, int, int, int, double, double);
extern void ofd_outdata(const char []);

// nearest integer
#define NINT(l, d) (int)((fabs(l) / fabs(d)) + 0.5)
#endif  // _OFD_DATALIB_H_
