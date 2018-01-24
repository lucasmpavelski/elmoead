#pragma once

#include "../surrogate/var/de_operators.h"

TEST(DifferentialEvolution, optimize) {
    DifferentialEvolution de(
        DEOperator(DEType::RAND_1_BIN, 1.0, 0.5),
        PMOperator(0.5));
    svectord lb = {-5, -5};
    svectord ub = { 5,  5};
    GenericProblem prob("test", 2, 1, lb, ub, [](const double* x, double* y) {
        y[0] = x[0] * x[0] + x[1] * x[1];
    }, [](double* x) {
        x[0] = MOProblem::truncateVar(x[0], -5, 5);
        x[1] = MOProblem::truncateVar(x[1], -5, 5);
    });
    Population pop(10, prob, LATIN_HYP_SAMPLING);
    pop.evaluateAll();
    de.optimize(10000, prob, pop, NullCallback);
    /*[](size_t no_evals, const Population& pop) {
        if (no_evals % 10 == 0)
            cout << "vars:\n" << pop.vars() << endl;
    });*/
    auto near = [](double a, double b) { return a - b < 0.01; };
    for (auto& s : pop)
        ASSERT_TRUE(near(s.vars[0], 0) && near(s.vars[1], 0));
}
