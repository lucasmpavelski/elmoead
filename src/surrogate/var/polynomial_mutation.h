#ifndef GENETIC_OPERATORS_INC
#define GENETIC_OPERATORS_INC

#include <cmath>
#include <functional>

#include "../aux.h"

#ifdef __cplusplus
extern "C" {
#endif

void polynomialMutation(const double in[], const double lb[], const double ub[],
                        const int no_vars, const double mut_rate,
                        const double eta_m, double out[]);

void SBXover(const double lb[], const double ub[], const int no_vars,
             const double xover_rate, const double eta_c,
             const double x1[], const double x2[],
             double offs1[], double offs2[]);

#ifdef __cplusplus
}
#endif


class SBXoverOperator
{
    const double xover_rate;
    const double eta_c;

public:
    SBXoverOperator(const double xover_rate = 1.0, const double eta_c = 20.0) :
        xover_rate(xover_rate), eta_c(eta_c)
    {}

    inline double* apply(const double lb[], const double ub[], const int no_vars,
                         const double x1[], const double x2[],
                         double offs1[], double offs2[]) const {
        SBXover(lb, ub, no_vars, xover_rate, eta_c, x1, x2, offs1, offs2);
        return offs1;
    }

    inline double* operator()(const double lb[], const double ub[], const int no_vars,
                              const double x1[], const double x2[],
                              double offs1[], double offs2[]) const {
        return apply(lb, ub, no_vars, x1, x2, offs1, offs2);
    }

    friend std::ostream& operator<<(std::ostream& os, const SBXoverOperator& pmo) {
        os << "SBXoverOperator:\n";
        os << "  xover_rate: " << pmo.xover_rate << "\n";
        os << "  eta_c: " << pmo.eta_c << endl;
        return os;
    }
};

class PMOperator
{
    const double mut_rate;
    const double eta_m;

public:
    explicit PMOperator(const double mut_rate, const double eta_m = 20.0) :
        mut_rate(mut_rate), eta_m(eta_m)
    {}

    explicit PMOperator(const int no_vars, const double eta_m = 20.0) :
        PMOperator{1.0 / no_vars, eta_m}
    {}

    inline double* apply(const double lb[], const double ub[], const int no_vars,
                         const double in[], double out[]) const {
        polynomialMutation(in, lb, ub, no_vars, mut_rate, eta_m, out);
        return out;
    }

    inline double* apply(const double lb[], const double ub[], const int no_vars,
                         const double in[]) const {
        double* out = new double[no_vars];
        return apply(lb, ub, no_vars, in, out);
    }

    inline double* operator()(const double lb[], const double ub[], const int no_vars,
                              const double in[], double out[]) const {
        return apply(lb, ub, no_vars, in, out);
    }

    inline double* operator()(const double lb[], const double ub[], const int no_vars,
                              const double in[]) const {
        return apply(lb, ub, no_vars, in);
    }

    inline double* mutate(const double lb[], const double ub[], const int no_vars,
                          double inout[]) const {
        return apply(lb, ub, no_vars, inout, inout);
    }

    friend std::ostream& operator<<(std::ostream& os, const PMOperator& pmo) {
        os << "PMOperator:\n";
        os << "  mut_rate: " << pmo.mut_rate << "\n";
        os << "  eta_m: " << pmo.eta_m << endl;
        return os;
    }
};

int tournament(const int competitors[], const int no_comp, const int t_size,
               const std::function<bool(const int&,const int&)>& compare);

#endif /* GENETIC_OPERATORS_INC */
