#include <mex.h>
#include "anisoDouble.h"

#include "defs.h"

#define NARGS 6  /* Number of input arguments */

#define ROWS rows
#define COLS cols
#define PAD  (int)pad_amt[0]

/* Looping anisotropic gaussian filtering. */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, 
	const mxArray *prhs[])
{
    /* Declarations*/
    char *msg;
    mxArray *gauss1, *gauss2, *t;
    double *out, *out_2, *tmp, *angle, *g1, *g2;
    double *img, *img_2, *pad_amt, *pad_pr;
    int i, j, ni, nj, crd;
    int cc = 0; 
    int rows, cols, rows_2, cols_2, g1Dim, g2Dim;
    double *avg, *sig_u, *sig_v;

    if (nrhs != NARGS) { 
	sprintf(msg, "%d input arguments required.\n", NARGS);
	mexErrMsgTxt( msg); 
    } else if (nlhs > 2) {
	mexErrMsgTxt("Too many output arguments.\n"); 
    } 

    /* We only handle doubles */
    if (!mxIsDouble(prhs[0]) 
	    || !mxIsDouble(prhs[1]) 
	    || !mxIsDouble(prhs[2])
	    || !mxIsDouble(prhs[3])
	    || !mxIsDouble(prhs[4])) {
	mexErrMsgTxt("Input image and stdevs should be double.\n");
    }
   
    img   = mxGetPr(prhs[0]);
    cols  = mxGetN(prhs[0]);
    rows  = mxGetM(prhs[0]);

    img_2   = mxGetPr(prhs[1]);
    cols_2  = mxGetN(prhs[1]);
    rows_2  = mxGetM(prhs[1]);

    sig_u = mxGetPr(prhs[2]);
    sig_v = mxGetPr(prhs[3]);
    angle = mxGetPr(prhs[4]);

    /* Pad amount */
    pad_amt = mxGetPr(prhs[5]);

    if (cols != cols_2 || cols != cols_2) {
	mexErrMsgTxt("Input images should be the same size.\n");
	return;
    }

    if (cols != mxGetN(prhs[2]) + (int)pad_amt[0] * 2 
	   || cols != mxGetN(prhs[3]) + (int)pad_amt[0] * 2
	   || cols != mxGetN(prhs[4]) + (int)pad_amt[0] * 2
	   || rows != mxGetM(prhs[2]) + (int)pad_amt[0] * 2
	   || rows != mxGetM(prhs[3]) + (int)pad_amt[0] * 2
	   || rows != mxGetM(prhs[4]) + (int)pad_amt[0] * 2) {
	mexErrMsgTxt("Std devs. should be same size as image - (2*pad)\n");
	return;
    }

    if (DEBUG) {
	printf("r: %d, c: %d\n", rows, cols);
    }

    /* Create output */ 
    plhs[0] = mxCreateDoubleMatrix(rows-(int)pad_amt[0] * 2, 
	    cols-(int)pad_amt[0] * 2, mxREAL);
    out = mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(rows-(int)pad_amt[0] * 2, 
	    cols-(int)pad_amt[0] * 2, mxREAL);
    out_2 = mxGetPr(plhs[1]);

    /* Loop over inputs, filtering */
    for (i = (int)PAD; i < ROWS-(int)PAD; i++) {
	for (j = (int)PAD; j < COLS-(int)PAD; j++) {
	/*--------------------------------------------------
	*     Calculate coords of angle/sigs.
	*--------------------------------------------------*/
	    ni = i-(int)PAD;
	    nj = j-(int)PAD;
	    crd = coord(ni, nj, COLS - (int)PAD * 2, ROWS - (int)PAD * 2);
	    if (DEBUG)
		printf("i: %d, j: %d, ni: %d, nj: %d, crd: %d\n", 
		    i, j, ni, nj, crd);
	    anisoGauss(img, img_2, ROWS, COLS, i, j, sig_u[crd], sig_v[crd], 
		    0.5, angle[crd], &out[crd], &out_2[crd]);

	}
    }

    return;
}
