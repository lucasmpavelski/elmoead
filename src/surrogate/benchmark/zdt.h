#ifndef _ZDT_H
#define _ZDT_H

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

void ZDT1(const double *x, const unsigned int n, double *f);
void ZDT2(const double *x, const unsigned int n, double *f);
void ZDT3(const double *x, const unsigned int n, double *f);
void ZDT4(const double *x, const unsigned int n, double *f);
void ZDT6(const double *x, const unsigned int n, double *f);

#ifdef __cplusplus
}
#endif

#endif /* _ZDT_H */
