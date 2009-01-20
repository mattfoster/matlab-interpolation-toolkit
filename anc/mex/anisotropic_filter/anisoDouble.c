#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defs.h"

/* Create 1/2 a gaussian function, using the given dimension and sigma */
void createGauss(int dim, double sigma, double g[])
{
    int 	i;
    double 	e;

    if (DEBUG)
	printf("Creating gaussian (size: %d, sig: %2.4f)\n", dim, sigma);

    e = sqrt(2*PI) * sigma; 

    for (i = 0; i < dim; i++) {
	g[i] = (double) exp(-0.5 * pow((double) i, 2) / pow(sigma, 2)) / e;
    }
}	

/* Return x filtered value for a given point. */
double g(double img[], int x, int y, int cols, int rows, double *g1, int g1Dim)
{
    int 	k;
    double 	out;

    if ( x-g1Dim < 0 || x > cols || y < 0 || y > rows) {
	if (DEBUG) printf("Index outside array.\n");
	return 0.0;
    }
   
    out = g1[0]*img[coord(x, y, cols, rows)];
    
    for (k = 1; k < g1Dim; k++) { /* Loop over filter */
	out += 
	    g1[k] * ((double)img[coord(x-k, y, cols, rows)] + 
		     (double)img[coord(x+k, y, cols, rows)]);
    }

    return out;
}

void anisoGauss(double img1[], double img2[], int rows, int cols, 
	int x, int y, double sig_u, double sig_v, double a,
	double angle, double *out1, double *out2) 
{

    double 	den, t_phi, phi, sig_x, sig_phi, den2;
    double	*w, *w1, *g1, *g2; 
    int 	i, k, xmfk, xpfk;
    int		g1Dim, g2Dim;

    /* Compute required parameters: */
    if (0 == angle) {
	t_phi 	= 0;
	phi 	= 0;
	sig_phi = sig_v;
	sig_x   = sig_u;
    } else {
	den   	= pow(sig_v, 2) * pow(cos(angle), 2) + 
	    pow(sig_u, 2) * pow(sin(angle), 2);

	den2 	= (pow(sig_u, 2) - pow(sig_v, 2)) * 
	    cos(angle)*sin(angle);

	if (fabs(den2) < TOL)
	    t_phi = 0;
	else
	    t_phi = den/den2;

	phi	= atan(t_phi);

	/* Filter standard deviations */
	if (fabs(phi) < TOL)
	    sig_phi = sig_u;
	else
	    sig_phi = fabs(1/sin(phi) * sqrt(den));

	sig_x	= fabs((sig_u * sig_v) / sqrt(den));

	sig_x 	= MIN(sig_x, SIG_CAP);
	sig_phi = MIN(sig_phi, SIG_CAP);
	
	if (DEBUG) {
	    printf("i: %d, j: %d\n", x, y);

	    printf("sig_u: %2.4f, sig_v: %2.4f, sig_x: %2.4f\n", 
		    sig_u, sig_v, sig_x);
		
	    printf("sig_phi: %2.4f, phi: %2.4f, t_phi: %2.4f\n", 
		    sig_phi, phi, t_phi);
	}
    }

    /* Create gaussians (use 3 sigma + 1 rule)*/
    g1Dim = 3 * sig_x + 1;
    g1 = (double *) calloc(g1Dim, sizeof(double));
    createGauss(g1Dim, sig_x, g1);

    g2Dim = 3 * sig_phi + 1;
    g2 = (double *) calloc(g2Dim, sizeof(double));
    createGauss(g2Dim, sig_phi, g2);   

    if ( x-g1Dim < 0 || y-g2Dim < 0 ) {
	/*--------------------------------------------------
	* printf("Error: x and/or y coordinates are too small.\n");
	*--------------------------------------------------*/
	return;
    }

    /* Allocate memory for premultiplied factors */
    w  = (double *) calloc(g2Dim, sizeof (double));
    w1 = (double *) calloc(g2Dim, sizeof (double));

    /* Premultiply the interpolation factors */
    for (i = 0; i < g2Dim; i++) {
	w[i]  = g2[i] * a;
	w1[i] = g2[i] * (1-a);
    }   

/*--------------------------------------------------
*     Only filter in x...
*--------------------------------------------------*/
/*--------------------------------------------------
*     out[0] =  g(img, x, y, cols, rows, g1, g1Dim);
*--------------------------------------------------*/

    out1[0] = g2[0]*g(img1, x, y, cols, rows, g1, g1Dim);
    out2[0] = g2[0]*g(img2, x, y, cols, rows, g1, g1Dim);

    /* Filter in u */
    for (k = 1; k < g2Dim; k++) {

	if (fabs(t_phi) < TOL) {
	    xpfk = x;
	    xmfk = x;
	    if (DEBUG)
		printf("Set interpolation factor to 0\n");
	} else {
	    xpfk = (int)floor(x + (double)k/t_phi);
	    xmfk = (int)floor(x - (double)k/t_phi);
	}
	
	/*--------------------------------------------------
	* out[0] += 
	*     w[k] * (g(img, xmfk,   y-k, cols, rows, g1, g1Dim) +
	* 	    g(img, xpfk,   y+k, cols, rows, g1, g1Dim))+
	*     w1[k]* (g(img, xmfk-1, y-k, cols, rows, g1, g1Dim) +
	* 	    g(img, xpfk+1, y+k, cols, rows, g1, g1Dim));
	*--------------------------------------------------*/
 
	out1[0] += 
	    w[k] * (g(img1, xmfk,   y-k, cols, rows, g1, g1Dim) +
		    g(img1, xpfk,   y+k, cols, rows, g1, g1Dim))+
	    w1[k]* (g(img1, xmfk,   y-k, cols, rows, g1, g1Dim) +
		    g(img1, xpfk,   y+k, cols, rows, g1, g1Dim));
	out2[0] += 
	    w[k] * (g(img2, xmfk,   y-k, cols, rows, g1, g1Dim) +
		    g(img2, xpfk,   y+k, cols, rows, g1, g1Dim))+
	    w1[k]* (g(img2, xmfk,   y-k, cols, rows, g1, g1Dim) +
		    g(img2, xpfk,   y+k, cols, rows, g1, g1Dim));
    }

    if (DEBUG) printf("out: %2.4f\n", out1[0]);
    if (DEBUG) printf("out: %2.4f\n", out2[0]);

    free(g1);
    free(g2);
    free(w);
    free(w1);
}
