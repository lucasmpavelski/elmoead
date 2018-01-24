#ifndef DTLZ_PROBLEMS_H
#define	DTLZ_PROBLEMS_H

#include <cmath>
#include <cstdio>
#include <algorithm>
using std::fill_n;

#include "../mo.h"
#include "dtlz.h"

struct DTLZ1Problem : public MOProblem
{
    DTLZ1Problem(const int no_vars, const int no_objs) :
        MOProblem("DTLZ1", no_vars, no_objs,
                  dvector(no_vars, 0.0), dvector(no_vars, 1.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        dtlz::DTLZ1(x, f, noVars(), noObjs());
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

struct DTLZ2Problem : public MOProblem
{
    DTLZ2Problem(const int no_vars, const int no_objs) :
        MOProblem("DTLZ2", no_vars, no_objs,
                  dvector(no_vars, 0.0), dvector(no_vars, 1.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        dtlz::DTLZ2(x, f, noVars(), noObjs());
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

struct DTLZ3Problem : public MOProblem
{
    DTLZ3Problem(const int no_vars, const int no_objs) :
        MOProblem("DTLZ3", no_vars, no_objs,
                  dvector(no_vars, 0.0), dvector(no_vars, 1.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        dtlz::DTLZ3(x, f, noVars(), noObjs());
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

struct DTLZ4Problem : public MOProblem
{
    DTLZ4Problem(const int no_vars, const int no_objs) :
        MOProblem("DTLZ4", no_vars, no_objs,
                  dvector(no_vars, 0.0), dvector(no_vars, 1.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        dtlz::DTLZ4(x, f, noVars(), noObjs());
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

struct DTLZ5Problem : public MOProblem
{
    DTLZ5Problem(const int no_vars, const int no_objs) :
        MOProblem("DTLZ5", no_vars, no_objs,
                  dvector(no_vars, 0.0), dvector(no_vars, 1.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        dtlz::DTLZ5(x, f, noVars(), noObjs());
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

struct DTLZ6Problem : public MOProblem
{
    DTLZ6Problem(const int no_vars, const int no_objs) :
        MOProblem("DTLZ6", no_vars, no_objs,
                  dvector(no_vars, 0.0), dvector(no_vars, 1.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        dtlz::DTLZ6(x, f, noVars(), noObjs());
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

struct DTLZ7Problem : public MOProblem
{
    DTLZ7Problem(const int no_vars, const int no_objs) :
        MOProblem("DTLZ7", no_vars, no_objs,
                  dvector(no_vars, 0.0), dvector(no_vars, 1.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        dtlz::DTLZ7(x, f, noVars(), noObjs());
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

#endif /* DTLZ_PROBLEMS_H */

