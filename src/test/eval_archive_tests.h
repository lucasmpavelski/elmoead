#include <cstdio>
#include <cmath>

#include "helpers.h"
#include "../surrogate/mo.h"
#include "../surrogate/aux.h"
/*
class EvalArchiveTest : public testing::Test
{
protected:
    MOProblem* prob;
    static const int no_vars = 2;
    static const int no_objs = 2;

    EvalArchive ea;
    Solution* sols[2];
    static const int no_sols = 2;
    Solution* sol_similar;
    Solution* sol_different;

    static void evaluate(void* problem, const double ins[],
                                   double outs[]) {
        outs[0] = ins[0] + ins[1];
        outs[1] = ins[0] - ins[1];
    };

    virtual void SetUp() {
        double lb[2] = {0, 1};
        double ub[2] = {1, 2};

        prob = newMOProblem("test", no_vars, no_objs, lb, ub, NULL, evaluate,
                            validateProbTruncating, NULL);

        sols[0] = newSolutionAlloc(prob, NULL);
        sols[1] = newSolutionAlloc(prob, NULL);
        sol_similar = newSolutionAlloc(prob, NULL);
        sol_different = newSolutionAlloc(prob, NULL);

        sols[0]->vars[0]       = 0.0; sols[0]->vars[1]       = 1.0;
        sols[1]->vars[0]       = 1.0; sols[1]->vars[1]       = 0.0;
        sol_similar->vars[0]   = 1.0; sol_similar->vars[1]   = 0.8;
        sol_different->vars[0] = 0.9; sol_different->vars[1] = 1.0;
        double eps = 0.9;

        newEvalArchive(prob, 4, eps, &ea);
    };

    void addSolutions() {
        for (int i = 0; i < no_sols; ++i)
            addSampleToArchive(&ea, sols[i]);
    };

    virtual void TearDown() {
        freeMOProblem(prob);
        for (int i = 0; i < no_sols; ++i)
            freeSolution(sols[i]);
        freeSolution(sol_similar);
        freeSolution(sol_different);
        freeEvalArchiveContents(&ea);
    };
};

TEST_F(EvalArchiveTest, Constructor) {
    EXPECT_EQ(ea.act_size, 0);
    for (int i = 0; i < no_sols; ++i)
        EXPECT_FALSE(archiveContains(&ea, sols[i]));
}

TEST_F(EvalArchiveTest, ContainSimilar) {
    addSolutions();
    EXPECT_TRUE(archiveContains(&ea, sol_similar));
}

TEST_F(EvalArchiveTest, NotContainDifferent) {
    addSolutions();
    EXPECT_FALSE(archiveContains(&ea, sol_different));
}

/*
TEST_F(PopulationTest, PrintObjs) {
    double v = 0.0;
    for (int i = 0; i < pop_size; ++i)
        for (int j = 0; j < no_objs; ++j)
            pop[i].objs[j] = v++;
    FILE* f = tmpfile();
    printPopulationObjs(f, pop, pop_size);
    char* f_str = fileToString(f);
    fclose(f);
    EXPECT_TRUE(isSubString(f_str,
                            "0.000000 1.000000\n"
                            "2.000000 3.000000\n"
                            "4.000000 5.000000"));
}

TEST_F(PopulationTest, Read) {
    FILE* f = tmpfile();
    fprintf(f,
            "0.000000 -1.000000;6.  7.00000001   \n"
            "  \n" // blank lines are ok
            "   2 3.00  ;8.000000  -9.000000\n");
    rewind(f);
    int ind_read = readPopulation(f, pop_size, pop);
    fclose(f);
    double v[3][2] = {{ 0.0, -1.0},
                      { 2.0, 3.0}};
    double o[3][2] = {{ 6.0, 7.00000001},
                      { 8.0, -9.0}};
    EXPECT_EQ(ind_read, 2); // only 2 solutions read
    for (int i = 0; i < ind_read; ++i)
    {
        EXPECT_TRUE(allDoubleEq(v[i], pop[i].vars, no_vars));
        EXPECT_TRUE(allDoubleEq(o[i], pop[i].objs, no_objs));
    }
}

TEST_F(SolutionTest, Eval) {
    sol->vars[0] = 0;
    sol->vars[1] = 1;
    evaluateSolution(sol);
    EXPECT_DOUBLE_EQ(sol->objs[0], +1.0);
    EXPECT_DOUBLE_EQ(sol->objs[1], -1.0);
}

TEST_F(SolutionTest, Validate) {
    sol->vars[0] = +10.0;
    sol->vars[1] = -10.0;
    validateSolution(sol);
    EXPECT_DOUBLE_EQ(sol->vars[0], 1.0);
    EXPECT_DOUBLE_EQ(sol->vars[1], 1.0);
}
*/
