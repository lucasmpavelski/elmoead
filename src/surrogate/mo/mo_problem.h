#ifndef MO_PROBLEM_H
#define MO_PROBLEM_H

#include <algorithm>
#include <functional>

#include "../aux.h"

struct MOProblem
{
public:
    using dvector = std::vector<double>;
    using bounds = std::vector<std::pair<double,double>>;

    MOProblem(const std::string& name, const size_t no_vars, const size_t no_objs,
              const dvector& lower_bounds, const dvector& upper_bounds) :
        no_vars(no_vars),
        no_objs(no_objs),
        lower_bounds(lower_bounds),
        upper_bounds(upper_bounds),
        name(name)
    {}

    const char* getName() const { return name.c_str(); }
    size_t noVars() const { return no_vars; }
    size_t noObjs() const { return no_objs; }

    bool isValid(const double x[]) const {
        for (int i = 0; i < no_vars; ++i)
            if (x[i] < lower_bounds[i] || x[i] > upper_bounds[i])
                return false;
        return true;
    }

    virtual void evaluate(const double x[], double f[]) = 0;
    virtual void validate(double x[]) = 0;

    inline void truncateToBounds(double x[]) const {
        for (int i = 0; i < no_vars; ++i)
            x[i] = truncateVar(x[i], lower_bounds[i], upper_bounds[i]);
    }

    inline void reflectToBounds(double x[]) const {
        for (int i = 0; i < no_vars; ++i)
            x[i] = reflectVar(x[i], lower_bounds[i], upper_bounds[i]);
    }

    static inline double reflectVar(double x, const double& lb,
                                    const double& ub) {
        while (x < lb)
            x = ub - std::fabs(x - lb);
        while (x > ub)
            x = lb + std::fabs(x - ub);
        return x;
    }

    static inline double truncateVar(double x, const double& lb,
                                     const double& ub) {
        if (x < lb)
            return lb;
        if (x > ub)
            return ub;
        return x;
    }

    bounds getBounds() const {
        bounds b(noVars());
        for (int i = 0; i < noVars(); i++)
            b[i] = {lower_bounds[i], upper_bounds[i]};
        return b;
    }

    const size_t no_vars;
    const size_t no_objs;
    const dvector lower_bounds;
    const dvector upper_bounds;
    const std::string name;

    friend std::ostream& operator<<(std::ostream& os, const MOProblem& p) {
        return os << "name: " << p.getName() << '\n'
                  << "no_vars: " << p.noVars() << '\n'
                  << "no_objs: " << p.noObjs() << '\n';
    }

};

struct TestProblem : public MOProblem
{
public:
    TestProblem(const int no_vars, const int no_objs, const double lb = 0.0,
                const double ub = 1.0) :
        MOProblem("TestProblem", no_vars, no_objs,
                  dvector(no_vars, lb), dvector(no_vars, ub))
    {};

    inline void evaluate(const double x[], double f[]) override {}

    inline void validate(double x[]) override {}
};

class GenericProblem : public MOProblem
{
public:
    typedef std::function<void(double const*,double*)> EvaluateFunc;
    typedef std::function<void(double*)> ValidateFunc;

    GenericProblem(const std::string& name,
                   const int& no_vars,
                   const int& no_objs,
                   const dvector& lower_bounds,
                   const dvector& upper_bounds,
                   EvaluateFunc eval, ValidateFunc val) :
        MOProblem(name, no_vars, no_objs, lower_bounds, upper_bounds),
        eval(eval),
        val(val)
    {}

    inline void evaluate(const double x[], double f[]) final override {
        eval(x, f);
    }

    inline void validate(double x[]) final override {
        val(x);
    }

private:
    EvaluateFunc eval;
    ValidateFunc val;
};

#endif /* MO_PROBLEM_H */
