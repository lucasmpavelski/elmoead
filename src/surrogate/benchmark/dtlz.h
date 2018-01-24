#ifndef _DTLZ_H
#define _DTLZ_H

#include <cmath>
#include <cstdlib>

namespace dtlz
{

void DTLZ1(const double *x, double* f, const unsigned no_vars, const unsigned no_objs);
void DTLZ2(const double *x, double* f, const unsigned no_vars, const unsigned no_objs);
void DTLZ3(const double *x, double* f, const unsigned no_vars, const unsigned no_objs);
void DTLZ4(const double *x, double* f, const unsigned no_vars, const unsigned no_objs);
void DTLZ5(const double *x, double* f, const unsigned no_vars, const unsigned no_objs);
void DTLZ6(const double *x, double* f, const unsigned no_vars, const unsigned no_objs);
void DTLZ7(const double *x, double* f, const unsigned no_vars, const unsigned no_objs);

}

#endif /* _DTLZ_H */
