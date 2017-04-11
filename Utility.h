#ifndef _UTILITY_HEADER_
#define _UTILITY_HEADER_

#include<vector>

// return largest positive root of quadratic equation with coefficients a, b, c
// if both roots negative or complex, return a really big number
double quad_solve( double a, double b, double c );

// returns the index of e in vector A: A must be in ascending order
int bin_search( std::vector<double> A, double e );

#endif
