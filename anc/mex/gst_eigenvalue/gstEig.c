/* Calculate the eigenvalues of a given set of multiplied gradients. */
/* These results aren't the same as matlabs..... */
 
#include <mex.h>
#include <math.h>

#include "defs.h"
#include "common.h"

#define NARGS 3  /* Number of input arguments */

#define ROWS rows
#define COLS cols

/* Looping anisotropic gaussian filtering. */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, 
	const mxArray *prhs[])
{
    /* Declarations*/

    char 	*msg;
    double 	*gx2, *gy2, *gxy, *lam_1, *lam_2;
    double	pre, post;   
    int 	i, j, rows, cols;

    if (nrhs != NARGS) { 
	sprintf(msg, "%d input arguments required.\n", NARGS);
	mexErrMsgTxt( msg); 
    } else if (nlhs > 2) {
	mexErrMsgTxt("Too many output arguments.\n"); 
    } 

    /* We only handle doubles */
    if (!mxIsDouble(prhs[0]) 
	    || !mxIsDouble(prhs[1]) 
	    || !mxIsDouble(prhs[2])) {
	mexErrMsgTxt("Input image and stdevs should be double.\n");
    }
   
    gxy   = mxGetPr(prhs[0]);
    cols  = mxGetN(prhs[0]);
    rows  = mxGetM(prhs[0]);
    
    gx2   = mxGetPr(prhs[1]);
    gy2   = mxGetPr(prhs[2]);

    if (cols != mxGetN(prhs[1]) || cols != mxGetN(prhs[2]) 
	   || rows != mxGetM(prhs[1]) || rows != mxGetM(prhs[2])) { 
	mexErrMsgTxt("Input vectors should be of equal size.\n");
    }


    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    lam_1 = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    lam_2 = mxGetPr(plhs[1]);

    for (i = 0; i < ROWS * COLS; i++) { 
	pre  = 0.5 * (gx2[i] + gy2[i]);
	post = 0.5 * sqrt(pow(gx2[i] - gy2[i], 2.0) + 4.0 * pow(gxy[i], 2.0));
	lam_1[i] = pre + post;
	lam_2[i] = pre - post;
    }
}
