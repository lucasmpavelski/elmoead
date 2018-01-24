#include "zdt.h"

#define PI 3.1415926535897932384626433832795
//#include <stdio.h>
//static int count_eval = 0;

void ZDT1(const double *x, const unsigned int n, double *f)
{

    unsigned int i = 0;
    double f1 = 0;
    double g = 0;
    double h = 0;


    f1 = x[0];

    for (i = 1; i < n; i++)
    {
        g += x[i];
    }
    g = 1 + 9 * g / (n-1);
    h = 1 - sqrt(f1 / g);

    f[0] = f1;
    f[1] = g * h;

}

void ZDT2(const double *x, const unsigned int n, double *f)
{
    unsigned int i = 0;
    double f1 = 0;
    double g = 0;
    double h = 0;


    f1 = x[0];

    for (i = 1; i < n; i++)
    {
        g += x[i];
    }
    g = 1.0 + 9.0 * g / (n-1.0);
    h = 1.0 - pow(f1 / g, 2.0);

    f[0] = f1;
    f[1] = g * h;
    //printf("%d\n", ++count_eval); fflush(stdout);

}

void ZDT3(const double *x, const unsigned int n, double *f)
{
    unsigned int i = 0;
    double f1 = 0;
    double g = 0;
    double h = 0;


    f1 = x[0];

    for (i = 1; i < n; i++)
    {
        g += x[i];
    }
    g = 1 + 9 * g / (n-1);
    h = 1 - sqrt(f1 / g) - (f1 / g) * sin(10 * PI * f1);

    f[0] = f1;
    f[1] = g * h;

}

void ZDT4(const double *x, const unsigned int n, double *f)
{
    unsigned int i = 0;
    double f1 = 0;
    double g = 0;
    double h = 0;

    f1 = x[0];

    for (i = 1; i < n; i++)
    {
        double x_i = x[i];
        g += x_i * x_i - 10.0 * cos(4.0 * PI * x_i);
    }
    g = 1.0 + 10.0 * (n - 1) + g;
    h = 1.0 - sqrt(f1 / g);

    f[0] = f1;
    f[1] = g * h;
}

void ZDT6(const double *x, const unsigned int n, double *f)
{
    unsigned int i = 0;
    double f1 = 0;
    double g = 0;
    double h = 0;

    f1 = 1 - exp(-4.0 * x[0]) * pow(sin(6.0 * PI * x[0]), 6.0);

    for (i = 1; i < n; i++)
    {
        g += x[i];
    }
    g = 1.0 + 9.0 * pow(g / (n-1.0), 0.25);
    h = 1.0 - pow(f1 / g, 2.0);

    f[0] = f1;
    f[1] = g * h;
}
