#pragma once

#include "wfg.h"
#include "../mo.h"

template<int n>
struct WFGProblem : public MOProblem {
    const size_t k;
    size_t evc;

    WFGProblem(const size_t no_vars, const size_t no_objs)
        : MOProblem("WFG" + std::to_string(n), no_vars, no_objs,
                    dvector(no_vars, 0.0), dvector(no_vars, 1.0))
        , k(initK(no_vars)) {}

    inline void evaluate(const double x[], double f[]) final {
        wfg_eval(x, noVars(), k, noObjs(), getName(), f);
    }

    inline void validate(double x[]) final { truncateToBounds(x); }

    static size_t initK(size_t no_vars) {
        size_t k = std::ceil(no_vars / 6.0);
        return (no_vars - k) % 2 == 0 ? k : k + 1;
    }
};

using WFG1Problem = WFGProblem<1>;
using WFG2Problem = WFGProblem<2>;
using WFG3Problem = WFGProblem<3>;
using WFG4Problem = WFGProblem<4>;
using WFG5Problem = WFGProblem<5>;
using WFG6Problem = WFGProblem<6>;
using WFG7Problem = WFGProblem<7>;
using WFG8Problem = WFGProblem<8>;
using WFG9Problem = WFGProblem<9>;



//struct WFG1Problem : public WFGProblem {
//    WFG1Problem(const size_t no_vars, const size_t no_objs)
//        : WFGProblem("WGF1", no_vars, no_objs)
//    {
//    }
//};
//
//struct WFG2Problem : public WFGProblem {
//    WFG2Problem(const size_t no_vars, const size_t no_objs)
//        : WFGProblem("WGF2", no_vars, no_objs)
//    {
//    }
//};
//
//struct WFG3Problem : public WFGProblem {
//    WFG3Problem(const size_t no_vars, const size_t no_objs)
//        : WFGProblem("WGF3", no_vars, no_objs)
//    {
//    }
//};
//
//struct WFG4Problem : public WFGProblem {
//    WFG4Problem(const size_t no_vars, const size_t no_objs)
//        : WFGProblem("WGF4", no_vars, no_objs)
//    {
//    }
//};
//
//struct WFG5Problem : public WFGProblem {
//    WFG5Problem(const size_t no_vars, const size_t no_objs)
//        : WFGProblem("WGF5", no_vars, no_objs)
//    {
//    }
//};
//
//struct WFG6Problem : public WFGProblem {
//    WFG6Problem(const size_t no_vars, const size_t no_objs)
//        : WFGProblem("WGF6", no_vars, no_objs)
//    {
//    }
//};
//
//struct WFG7Problem : public WFGProblem {
//    WFG7Problem(const size_t no_vars, const size_t no_objs)
//        : WFGProblem("WGF7", no_vars, no_objs)
//    {
//    }
//};
//
//struct WFG8Problem : public WFGProblem {
//    WFG8Problem(const size_t no_vars, const size_t no_objs)
//        : WFGProblem("WGF8", no_vars, no_objs)
//    {
//    }
//};
//
//struct WFG9Problem : public WFGProblem {
//    WFG9Problem(const size_t no_vars, const size_t no_objs)
//        : WFGProblem("WGF9", no_vars, no_objs)
//    {
//    }
//};
