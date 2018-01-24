#include <cstdio>
#include <cmath>

#include "helpers.h"
#include "../surrogate/moead/stm.h"
#include "../surrogate/aux.h"

class STMTest : public testing::Test {
protected:
    MOProblem* prob;
    static const int no_vars = 2;
    static const int no_objs = 2;

    virtual void SetUp()
    {
        prob = new TestProblem(no_vars, no_objs);
    }

    virtual void TearDown()
    {
        delete prob;
    }
};

TEST_F(STMTest, ComputeProbPreferences)
{
    const int no_objs = 2;
    const int no_problems = 3;
    const int no_solutions = 3;

    TestProblem prob(1, no_objs);
    WeightSet<double> w(no_objs, no_solutions - 1, WeightSet<>::Type::UNIFORM, 0.0);
    AggrProblem<double> aggr_prob(INVERTED_TCHEBYCHEFF, no_objs);

    const double ideal[no_objs] = { 0.0, 0.0 };
    const double nadir[no_objs] = { 1.0, 1.0 };
    aggr_prob.updateReferencePoints(ideal);
    aggr_prob.updateReferencePoints(nadir);

    double pop_objs[no_solutions][no_objs] = {
        { 0.1, 0.9 },
        { 0.4, 0.6 },
        { 0.9, 0.1 },
    };
    Population pop(no_solutions, prob);
    for (int i = 0; i < no_solutions; ++i)
        pop[i].setObjs(pop_objs[i]);

    ematrixi out;
    problemPreferences(pop, w, aggr_prob, out);

    ematrixi expected(no_problems, no_solutions);
    expected << 0, 1, 2,
        1, 0, 2,
        2, 1, 0;

    EXPECT_TRUE(out == expected);
}

TEST_F(STMTest, ComputeSolsPreferences)
{
    const int no_objs = 2;
    const int no_problems = 3;
    const int no_solutions = 3;

    TestProblem prob(1, no_objs);
    WeightSet<double> w(no_objs, no_solutions - 1, WeightSet<>::Type::UNIFORM, 0.0);
    AggrProblem<double> aggr_prob(INVERTED_TCHEBYCHEFF, no_objs);

    const double ideal[no_objs] = { 0.0, 0.0 };
    const double nadir[no_objs] = { 1.0, 1.0 };
    aggr_prob.updateReferencePoints(ideal);
    aggr_prob.updateReferencePoints(nadir);

    double pop_objs[no_solutions][no_objs] = {
        { 0.9, 0.1 },
        { 0.4, 0.6 },
        { 0.1, 0.9 },
    };
    Population pop(no_solutions, prob);
    for (int i = 0; i < no_solutions; ++i)
        pop[i].setObjs(pop_objs[i]);

    ematrixi out;
    solutionPreferences(pop, w, aggr_prob, out);

    ematrixi expected(no_solutions, no_problems);
    expected << 2, 1, 0,
        1, 0, 2,
        0, 1, 2;

    EXPECT_TRUE(out == expected);

    /*   const int no_objs = 2;
    const int no_problems = 3;
    const int no_solutions = 3;
    Weights<double> w(no_objs, no_solutions-1, 0.0);
    const double ideal[no_objs] = {0.0, 0.0};
    const double nadir[no_objs] = {1.0, 1.0};
    double pop_objs[no_solutions][no_objs] = {
        {0.9, 0.1},
        {0.4, 0.6},
        {0.1, 0.9},
    };
    Solution* pop = newPopulationForProblem(no_solutions, prob);
    for (int i = 0; i < no_solutions; ++i)
        memcpy(pop[i].objs, pop_objs[i], no_objs * sizeof(double));
    int** p = computeSolsPreferences(pop, no_solutions, w, no_problems, ideal,
                                     nadir, NULL, NULL);

    const int gp[no_problems][no_solutions] = {
        {0, 1, 2},
        {1, 2, 0},
        {2, 1, 0},
    };
    for (int i = 0; i < no_problems; ++i)
    {
        EXPECT_TRUE(allEq(p[i], gp[i], no_solutions));
        //printIntVector(stdout, p[i], no_problems);
        //printf("\n");
    }

    freeMatrix(p);
    freePopulation(pop);*/
}
/*
TEST_F(STMTest, STM) {
    const int pop_size = 10;
    const int no_problems = 5;

    const int spsi_p[no_problems][pop_size] = {
        {1, 3, 4, 2, 5, 8, 7, 6, 9,10},
        {1, 4, 3, 2, 5, 8, 7, 6, 9,10},
        {2, 1, 5, 8, 4, 7, 3, 6, 9,10},
        {2, 8, 9,10, 1, 5, 7, 4, 6, 3},
        {9, 2,10, 8, 1, 5, 7, 4, 6, 3},
    };
    const int spsi_x[pop_size][no_problems] = {
        {1, 2, 3, 4, 5},
        {4, 5, 3, 2, 1},
        {1, 2, 3, 4, 5},
        {1, 2, 3, 4, 5},
        {2, 3, 1, 4, 5},
        {3, 4, 2, 5, 1},
        {3, 4, 2, 5, 1},
        {4, 5, 3, 2, 1},
        {5, 4, 3, 2, 1},
        {5, 4, 3, 2, 1},
    };
    int** psi_p = newMatrix<int>(no_problems, pop_size);
    int** psi_x = newMatrix<int>(pop_size, no_problems);

    for (int i = 0; i < no_problems; ++i)
    {
        for (int j = 0; j < pop_size; ++j)
        {
            psi_p[i][j] = spsi_p[i][j] - 1;
            psi_x[j][i] = spsi_x[j][i] - 1;
        }
    }

    int f_p[no_problems];
    int pc[no_problems];
    bool f_x[pop_size];
    int selected[pop_size];
    stm(pop_size, no_problems, psi_p, psi_x, pc, f_p, f_x,
        nullptr, nullptr, TCHEBYCHEFF, nullptr,
        selected);

    const int g_selected[pop_size] = {0, 3, 6, 1, 2, 6, 6, 6, 4, 6};
    EXPECT_TRUE(allEq(selected, g_selected, pop_size));

    freeMatrix(psi_p);
    freeMatrix(psi_x);
}
*/
