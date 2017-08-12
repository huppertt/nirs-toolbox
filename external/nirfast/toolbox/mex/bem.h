/* BEM.h is a compilation of all the functions and macros used
 * by the BEM (boundary element method) code located in the C files
 * main_build_matrix.c and build_source.c.
 * 
 * BY: Senate Taka
 * 
 * Created: Sept-08-2008
 * Modified: Hamid Ghadyani
 *				March 2010
 *           	  Most of the input arguments are not duplicated within mex file anymore.
 *           	  Fixed seg. fault crashes under Windows
 *				Apr 2010: 
 *				  Cleaned up usage of macros to access matrices.
 *				  Removed unnecessary copying of elements of integral results.
 *				May 2010:
 *				  Added OpenMP to significantly reduce the runtime
 */

#ifndef BEM_H
#define BEM_H

#include <math.h>
#include <mex.h>
#include <assert.h>
#define PI 3.141592653589793115997963468544185161590576171875

/* macros used for multiplying and dividing complex numbers.
 *  form complex numbers z = a+bi and w = c+di,
 *  multiplication: z*w = (a*c - b*d) + (a*d + b*c)*i
 *  stores the real and img parts in real & img.
 *  we use a do while loop with condition of 0 so that semicolons can be
 *  neatly hidden in the do part of the statement. NOTE, use macro by calling
 * mult_complex(...); (note semicolon at the end of macro call).
 */
#define mult_complex(a, b, c, d, real, img) \
do{ 	(real) = ( ((a)*(c)) - ((b)*(d)) ); \
	     (img) = ( ((a)*(d)) + ((b)*(c)) ); \
}\
while(0)

/* macro for cis(x), where x is a complex number a+bi
 * cis(yi) for real y = cos(y)+isin(y) (in other words, euler's formula)
 * _a_ is a pseudo gensym() variable name (the gensym in lisp).
 */
#define cis(a, b, real, img) \
do {	double _a_ = exp((a)); \
	(real) = _a_ * (cos((b))); \
	(img) = _a_ * (sin((b))); \
}\
while(0)


/* functions called by mexFunction */
void get_gausspts_3d(double *eta1, double *eta2, double *w, int ngp);
void get_area_normals(double *elem_list, double *nod_list, long num_elem,
			long num_nodes, double *area_elem_list, double *normal_x_list,
			double *normal_y_list, double *normal_z_list);
void singular_integrand_polar(double* int_val_real, double* int_val_img,
					double omegar, double omegai, double diff_coeff,
					double area, double *gp1, double *gp2, double *wgt, long ngpts);
static void compute_integral_gausspts(double *AAr, double *AAi, double *BBr, double *BBi,
		long size, long ngpts, 
		double *eta1, double *eta2, double *w, long ith_nod,
		long *jelem, double *xelem, double *yelem, double *zelem,
		double xx, double yy, double zz, double area, double norm_x, 
		double norm_y, double norm_z, double omegar, double omegai, double diff_coeff);

#endif /* bem_h */
