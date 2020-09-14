#ifndef TH_GENERIC_FILE
#define TH_GENERIC_FILE "generic/THGamma.h"
#else

/* Level 1 */
// TH_API void THBlas_(swap)(long n, real *x, long incx, real *y, long incy);

TH_API void THTensor_(polygamma)(long n, THTensor *input, THTensor *output);
TH_API void THTensor_(lbeta)(THTensor *a, THTensor *b, THTensor *output);


#endif
