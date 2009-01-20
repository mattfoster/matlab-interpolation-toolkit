/* anisoDouble.c */
void createGauss(int dim, double sigma, double g[]);
double g(double img[], int x, int y, int cols, int rows, double *g1, int g1Dim);
void anisoGauss(double img1[], double img2[], int rows, int cols, int x, int y, double sig_u, double sig_v, double a, double angle, double *out1, double *out2);
