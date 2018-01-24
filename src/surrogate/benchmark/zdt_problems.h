#ifndef ZDT_PROBLEMS_H
#define	ZDT_PROBLEMS_H

#include <cmath>
#include <cstdio>
#include <algorithm>
using std::fill_n;

#include "../mo.h"
#include "zdt.h"


struct ZDT1Problem : public MOProblem
{
    ZDT1Problem(const int no_vars) :
        MOProblem("ZDT1", no_vars, 2,
                  dvector(no_vars, 0.0), dvector(no_vars, 1.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        double g = x[1];
        for (int i = 2; i < no_vars; i++)
            g += x[i];
        g = 1.0 + 9.0 * g / (no_vars - 1.0);
        double h = 1 - sqrt(x[0] / g);
        f[0] = x[0];
        f[1] = g * h;
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

struct ZDT2Problem : public MOProblem
{
    ZDT2Problem(const int no_vars) :
        MOProblem("ZDT2", no_vars, 2,
                   dvector(no_vars, 0.0), dvector(no_vars, 1.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        double g = x[1];
        for (int i = 2; i < no_vars; ++i)
            g += x[i];
        g = 1.0 + 9.0 * g / (no_vars - 1.0);
        const double h = 1.0 - sqr(x[0] / g);
        f[0] = x[0];
        f[1] = g * h;
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

struct ZDT3Problem : public MOProblem
{
    ZDT3Problem(const int no_vars) :
        MOProblem("ZDT3", no_vars, 2,
                   dvector(no_vars, 0.0), dvector(no_vars, 1.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        double g = x[1];
        for (int i = 2; i < no_vars; i++)
            g += x[i];
        g = 1.0 + 9.0 * g / (no_vars - 1.0);
        const double h = 1.0 - sqrt(x[0] / g) - (x[0] / g) *
                sin(10 * M_PI * x[0]);
        f[0] = x[0];
        f[1] = g * h;
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

struct ZDT4Problem : public MOProblem
{
    ZDT4Problem(const int no_vars) :
        MOProblem("ZDT4", no_vars, 2, make_lb(no_vars), make_ub(no_vars))
    {
    };

    static dvector make_lb(const int no_vars) {
        dvector v(no_vars, -5.0);
        v[0] = 0.0;
        return v;
    };

    static dvector make_ub(const int no_vars) {
        dvector v(no_vars, 5.0);
        v[0] = 1.0;
        return v;
    };

    inline void evaluate(const double x[], double f[]) final {
        double g = 0.0;
        for (int i = 1; i < no_vars; ++i)
            g += x[i] * x[i] - 10.0 * cos(4.0 * M_PI * x[i]);
        g = 1.0 + 10.0 * (no_vars - 1) + g;
        const double h = 1.0 - sqrt(x[0] / g);
        f[0] = x[0];
        f[1] = g * h;
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

struct ZDT6Problem : public MOProblem
{
    ZDT6Problem(const int no_vars) :
        MOProblem("ZDT6", no_vars, 2,
                   dvector(no_vars, 0.0), dvector(no_vars, 1.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        double g = x[1];
        for (int i = 2; i < no_vars; ++i)
            g += x[i];
        g = 1.0 + 9.0 * pow(g / (no_vars - 1.0), 0.25);
        const double f1 = 1.0 - exp(-4.0 * x[0]) *
                intpow(sin(6.0 * M_PI * x[0]), 6);
        const double h = 1.0 - pow(f1 / g, 2.0);
        f[0] = f1;
        f[1] = g * h;
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};
/*
#ifdef __cplusplus
extern "C" {
#endif

void initKNO1Problem(MOProblem* prob);
void deinitKNO1Problem(MOProblem* prob);

// ZDTs

MOProblem* newZDTProblem(const char* name, const int no_vars);

void freeZDTProblem(MOProblem* prob);

#ifdef __cplusplus
}
#endif
*/
#endif /* ZDT_PROBLEMS_H */

