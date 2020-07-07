/*******************************************************************
 * NAME : void initmod(void (* odeparms)(int *, double *))
 *
 * DESCRIPTION :     initialises the parameter values, passed from the R- code
 *
 * INPUTS :
 *       PARAMETERS:
 *           void  * odeparms              fills a double array with double
 *                                         precision values, to copy the parameter
 *                                         values into the global variable
 */
#include "parameters.h"

void initmod(void (* odeparms)(int *, double *))
{
  int N = 48;
  odeparms(&N, p.value);
}