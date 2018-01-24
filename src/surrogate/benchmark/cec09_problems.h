#ifndef CEC09_PROBLEMS_H
#define CEC09_PROBLEMS_H

#include "../mo.h"
#include "../aux.h"
#include "cec09.hpp"

static MOProblem::dvector make2dUFBound(const int no_vars, const double first,
                                        const double remaining)
{
    MOProblem::dvector r(no_vars, remaining);
    r[0] = first;
    return r;
}

static MOProblem::dvector make3dUFBound(const int no_vars, const double first,
                                        const double second, const double remaining)
{
    MOProblem::dvector r(no_vars, remaining);
    r[0] = first;
    r[1] = second;
    return r;
}

struct UF1Problem : public MOProblem
{
    UF1Problem(const int no_vars) :
        MOProblem("UF1", no_vars, 2,
                  make2dUFBound(no_vars, 0.0, -1.0), dvector(no_vars, 1.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        CEC09::UF1(x, f, no_vars);
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

struct UF2Problem : public MOProblem
{
    UF2Problem(const int no_vars) :
        MOProblem("UF2", no_vars, 2,
                  make2dUFBound(no_vars, 0.0, -1.0), dvector(no_vars, 1.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        CEC09::UF2(x, f, no_vars);
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

struct UF3Problem : public MOProblem
{
    UF3Problem(const int no_vars) :
        MOProblem("UF3", no_vars, 2,
                  dvector(no_vars, 0.0), dvector(no_vars, 1.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        CEC09::UF3(x, f, no_vars);
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

struct UF4Problem : public MOProblem
{
    UF4Problem(const int no_vars) :
        MOProblem("UF4", no_vars, 2,
                  make2dUFBound(no_vars, 0.0, -2.0),
                  make2dUFBound(no_vars, 1.0,  2.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        CEC09::UF4(x, f, no_vars);
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

struct UF5Problem : public MOProblem
{
    UF5Problem(const int no_vars) :
        MOProblem("UF5", no_vars, 2,
                  make2dUFBound(no_vars, 0.0, -1.0), dvector(no_vars, 1.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        CEC09::UF5(x, f, no_vars);
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

struct UF6Problem : public MOProblem
{
    UF6Problem(const int no_vars) :
        MOProblem("UF6", no_vars, 2,
                  make2dUFBound(no_vars, 0.0, -1.0), dvector(no_vars, 1.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        CEC09::UF6(x, f, no_vars);
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

struct UF7Problem : public MOProblem
{
    UF7Problem(const int no_vars) :
        MOProblem("UF7", no_vars, 2,
                  make2dUFBound(no_vars, 0.0, -1.0), dvector(no_vars, 1.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        CEC09::UF1(x, f, no_vars);
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

struct UF8Problem : public MOProblem
{
    UF8Problem(const int no_vars) :
        MOProblem("UF8", no_vars, 2,
                  make3dUFBound(no_vars, 0.0, 0.0, -2.0),
                  make3dUFBound(no_vars, 1.0, 1.0,  2.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        CEC09::UF1(x, f, no_vars);
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

struct UF9Problem : public MOProblem
{
    UF9Problem(const int no_vars) :
        MOProblem("UF9", no_vars, 2,
                  make3dUFBound(no_vars, 0.0, 0.0, -2.0),
                  make3dUFBound(no_vars, 1.0, 1.0,  2.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        CEC09::UF9(x, f, no_vars);
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

struct UF10Problem : public MOProblem
{
    UF10Problem(const int no_vars) :
        MOProblem("UF10", no_vars, 2,
                  make3dUFBound(no_vars, 0.0, 0.0, -2.0),
                  make3dUFBound(no_vars, 1.0, 1.0,  2.0))
    {};

    inline void evaluate(const double x[], double f[]) final {
        CEC09::UF10(x, f, no_vars);
    };

    inline void validate(double x[]) final { truncateToBounds(x); };
};

#endif /* CEC09_PROBLEMS_H */
