/* Functions to implement fast anisotropic gaussian filtering */

/* Return the index for the given coordinates (assuming matlab 
 * style indexing
 * */ 
int coord(int i, int j, int cols, int rows)
{
    int c;
    c = j * rows + i;
/*--------------------------------------------------
*     c = i * cols + j;
*--------------------------------------------------*/

/*--------------------------------------------------
*     return c < rows * cols ? c : rows * cols -1;
*--------------------------------------------------*/
    return c;
}
