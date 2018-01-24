#include "helpers.h"
#include "../surrogate/mo.h"
#include "../surrogate/benchmark.h"
#include "../surrogate/aux.h"

class ZDTTest : public testing::Test
{
protected:
    double* vars;
    double* objs;
    MOProblem* prob;

    virtual void SetUp() {
    };

    void setProblem(const char* nm) {
        prob = getBenchmark(nm, 2);
        //sol.problem = *prob;
    };

    virtual void TearDown() {
        delete prob;
    };
};
/*
void t1() {
    const int no_vars = 30;
    prob_ptr prob = getBenchmark("ZDT1", no_vars);
    Solution2 sol(prob);
    Solution::VarsT v[30]{1.0};
    sol.setVars(v);

    for (int i = 0; i < 1000000; ++i)
    {
     //   for (int j = 0; j < no_vars; ++j)
     //       sol.vars[j] = RNG::realUniform(0.0, 1.0);
      //  sol.validate();
        sol.evaluate();
    }
};

void t2() {
    const int no_vars = 30;
    MOProblem* oprob = newBenchmarkProblem("ZDT1", no_vars);
    Solution sol(*oprob);
    Solution::VarsT v[30]{1.0};
    sol.setVars(v);

    for (int i = 0; i < 1000000; ++i)
    {
    //    for (int j = 0; j < no_vars; ++j)
    //        sol.vars[j] = RNG::realUniform(0.0, 1.0);
     //   sol.validate();
        sol.evaluate();
    }

    freeMOProblem(oprob);
};
*/

TEST_F(ZDTTest, ZDT1) {
    setProblem("ZDT1");
    AllocSolution sol(*prob);

    sol.vars[0] = 0; sol.vars[1] = 0;
    sol.evaluate();
    EXPECT_DOUBLE_EQ(sol.objs[0], 0.0);
    EXPECT_DOUBLE_EQ(sol.objs[1], 1.0);

    sol.vars[0] = .5; sol.vars[1] = .5;
    sol.evaluate();
    EXPECT_DOUBLE_EQ(sol.objs[0], .5);
    EXPECT_NEAR(sol.objs[1], 3.841687605, 0.000001);

    sol.vars[0] = 1; sol.vars[1] = 1;
    sol.evaluate();
    EXPECT_DOUBLE_EQ(sol.objs[0], 1.0);
    EXPECT_NEAR(sol.objs[1], 6.83772234, 0.000001);
}

TEST_F(ZDTTest, ZDT2) {
    setProblem("ZDT2");
    AllocSolution sol(*prob);

    sol.vars[0] = 0; sol.vars[1] = 0;
    sol.evaluate();
    EXPECT_DOUBLE_EQ(sol.objs[0], 0.0);
    EXPECT_DOUBLE_EQ(sol.objs[1], 1.0);

    sol.vars[0] = .5; sol.vars[1] = .5;
    sol.evaluate();
    EXPECT_DOUBLE_EQ(sol.objs[0], .5);
    EXPECT_NEAR(sol.objs[1], 5.454545455, 0.000001);

    sol.vars[0] = 1; sol.vars[1] = 1;
    sol.evaluate();
    EXPECT_DOUBLE_EQ(sol.objs[0], 1.0);
    EXPECT_DOUBLE_EQ(sol.objs[1], 9.9);
}

TEST_F(ZDTTest, ZDT3) {
    setProblem("ZDT3");
    AllocSolution sol(*prob);

    sol.vars[0] = 0; sol.vars[1] = 0;
    sol.evaluate();
    EXPECT_DOUBLE_EQ(sol.objs[0], 0.0);
    EXPECT_DOUBLE_EQ(sol.objs[1], 1.0);

    sol.vars[0] = .5; sol.vars[1] = .5;
    sol.evaluate();
    EXPECT_DOUBLE_EQ(sol.objs[0], .5);
    EXPECT_NEAR(sol.objs[1], 3.841687605, 0.000001);

    sol.vars[0] = 1; sol.vars[1] = 1;
    sol.evaluate();
    EXPECT_DOUBLE_EQ(sol.objs[0], 1.0);
    EXPECT_NEAR(sol.objs[1], 6.83772234, 0.000001);
}

TEST_F(ZDTTest, ZDT4) {
    setProblem("ZDT4");
    AllocSolution sol(*prob);

    sol.vars[0] = 0; sol.vars[1] = 0;
    sol.evaluate();
    EXPECT_DOUBLE_EQ(sol.objs[0], 0.0);
    EXPECT_DOUBLE_EQ(sol.objs[1], 1.0);

    sol.vars[0] = .5; sol.vars[1] = .5;
    sol.evaluate();
    EXPECT_DOUBLE_EQ(sol.objs[0], .5);
    EXPECT_NEAR(sol.objs[1], 0.459430585, 0.000001);

    sol.vars[0] = 1; sol.vars[1] = 1;
    sol.evaluate();
    EXPECT_DOUBLE_EQ(sol.objs[0], 1.0);
    EXPECT_NEAR(sol.objs[1], 0.585786438, 0.000001);
}

TEST_F(ZDTTest, ZDT6) {
    setProblem("ZDT6");
    AllocSolution sol(*prob);

    sol.vars[0] = 0; sol.vars[1] = 0;
    sol.evaluate();
    EXPECT_DOUBLE_EQ(sol.objs[0], 1.0);
    EXPECT_DOUBLE_EQ(sol.objs[1], 0.0);

    sol.vars[0] = .5; sol.vars[1] = .5;
    sol.evaluate();
    EXPECT_DOUBLE_EQ(sol.objs[0], 1.0);
    EXPECT_NEAR(sol.objs[1], 8.451355308, 0.000001);

    sol.vars[0] = 1; sol.vars[1] = 1;
    sol.evaluate();
    EXPECT_DOUBLE_EQ(sol.objs[0], 1.0);
    EXPECT_DOUBLE_EQ(sol.objs[1], 9.9);
}

