#include "config.h"
#include "nan.h"
#include "minell.h"
#include "nn.h"

#include <mex.h>

#if !defined(_DBL_ARRAY_)
#define _DBL_ARRAY
typedef struct {
  double *pr;
  int n;
  int m;
  int s;
} dbl_array;
#endif

void load_dbl_array (const mxArray *pr, dbl_array *arr) 
{ 
  arr->pr = mxGetPr(pr);
  arr->m  = mxGetM(pr); 
  arr->n  = mxGetN(pr);
  arr->s  = arr->n * arr->m;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    dbl_array   x_in, y_in, z_in, x_out, y_out;
    point       *in_points = NULL;
    point       *out_points = NULL;
    delaunay    *d = NULL;
    void        *interpolator = NULL;
    double      *x_out_pr, *y_out_pr, *z_out_pr;

    int         ii;

    /* Declarations */
    if (nrhs != 5) { 
        mexErrMsgTxt("Five input arguments required.\n"); 
    } else if (nlhs > 3) {
        mexErrMsgTxt("Too many output arguments.\n"); 
    } 

    /* We only handle doubles */
    if (!mxIsDouble(prhs[0])) {
        mexErrMsgTxt("Input image should be double.\n");
    }

    /* load array and size into struct */
    load_dbl_array(prhs[0], &x_in);
    load_dbl_array(prhs[1], &y_in);
    load_dbl_array(prhs[2], &z_in);
    load_dbl_array(prhs[3], &x_out);
    load_dbl_array(prhs[4], &y_out);

    /* Check intput sizes */
    if (x_in.s != y_in.s || x_in.s != z_in.s)
        mexErrMsgTxt("x, y, z should be equal sized.\n");

    /* Check output sizes */
    if (x_out.s != y_out.s)
        mexErrMsgTxt("output coordinates should be equal sized.\n");

    in_points = malloc(x_in.s * sizeof(point));

    /* Convert inputs to array of points. */
    for (ii = 0; ii < x_in.s; ii++) {
        in_points[ii].x = x_in.pr[ii];
        in_points[ii].y = y_in.pr[ii];
        in_points[ii].z = z_in.pr[ii];
    }

    out_points = malloc(x_out.s * sizeof(point));
    
    /* Convert outputs to array of points. */
    for (ii = 0; ii < x_out.s; ii++) {
        out_points[ii].x = x_out.pr[ii];
        out_points[ii].y = y_out.pr[ii];
        out_points[ii].z = NaN;
    }

    /* Create delaunay triangulation*/
    d = delaunay_build(x_in.s, in_points, 0, NULL, 0, NULL);

    
    /* Create interpolator */    
    interpolator = nnpi_create(d);
    nnpi_setwmin(interpolator, -DBL_MAX);

    /* Serially interpolate */
    for (ii = 0; ii < x_out.s; ii++) {
        nnpi_interpolate_point(interpolator, &out_points[ii]);
    }

    /* Convert back to mxArray format */

    /* Allocate everything */

/*--------------------------------------------------
*     if (nlhs == 3) {
*         plhs[0] = mxCreateDoubleMatrix(x_out.m, x_out.n, mxREAL);
*         plhs[1] = mxCreateDoubleMatrix(y_out.m, y_out.n, mxREAL);
*         plhs[2] = mxCreateDoubleMatrix(z_out.m, z_out.n, mxREAL);
*         x_out_pr = mxGetPr(plhs[0]);
*         y_out_pr = mxGetPr(plhs[1]); 
*         z_out_pr = mxGetPr(plhs[2]);
*         for (ii = 0; ii < x_out.s; ii++) {
*             x_out_pr[ii] = out_points[ii].x;
*             y_out_pr[ii] = out_points[ii].y;
*             z_out_pr[ii] = out_points[ii].z;
*         }
*     } else 
*--------------------------------------------------*/
    if (nlhs == 1) {
        plhs[0] = mxCreateDoubleMatrix(x_out.m, x_out.n, mxREAL);
        z_out_pr = mxGetPr(plhs[0]);

        for (ii = 0; ii < x_out.s; ii++) {
            z_out_pr[ii] = out_points[ii].z;
        }
    } else {
        mexErrMsgTxt("Number of output arguments should be 1 or 3.\n");
    }




    /* Free everything */
    nnpi_destroy(interpolator);
    delaunay_destroy(d);
    free(in_points);
    free(out_points);
}
