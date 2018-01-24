#include <cstdio>

#include "helpers.h"
#include "../surrogate/mo.h"
#include "../surrogate/aux.h"


TEST(SolutionTests, EmptyConstructor) {
    Solution s;
    EXPECT_EQ(s.vars, nullptr);
    EXPECT_EQ(s.objs, nullptr);
    EXPECT_EQ(s.noVars(), 0);
    EXPECT_EQ(s.noObjs(), 0);
}

TEST(SolutionTests, SharedVarsConstructor) {
    double vars[2];
    double objs[3];
    TestProblem prob(2, 3);
    Solution s(&prob, vars, objs);
    EXPECT_EQ(s.vars, vars);
    EXPECT_EQ(s.objs, objs);
    EXPECT_EQ(s.noVars(), prob.noVars());
    EXPECT_EQ(s.noObjs(), prob.noObjs());
}

TEST(SolutionTests, AllocatedSolution) {
    TestProblem prob(2, 3);
    AllocSolution s(prob);
    EXPECT_NE(s.vars, nullptr);
    EXPECT_NE(s.objs, nullptr);
    EXPECT_EQ(s.noVars(), prob.noVars());
    EXPECT_EQ(s.noObjs(), prob.noObjs());
}

TEST(SolutionTests, CopyConstructor) {
    double vars[2];
    double objs[3];
    TestProblem prob(2, 3);
    Solution s2(&prob, vars, objs);
    Solution s(s2);
    EXPECT_EQ(s.vars, s2.vars);
    EXPECT_EQ(s.objs, s2.objs);
    EXPECT_EQ(s.problem, s2.problem);
}

TEST(SolutionTests, Assignment) {
    TestProblem prob(2, 3);
    AllocSolution s(prob);
    AllocSolution s2(prob);
    s.setVars({0, 1});
    s.setObjs({2, 3, 4});
    s2 = s;
    EXPECT_NE(s.vars, s2.vars);
    EXPECT_NE(s.objs, s2.objs);
    EXPECT_TRUE(std::equal(s.vars, s.vars + 2, s2.vars));
    EXPECT_TRUE(std::equal(s.objs, s.objs + 3, s2.objs));
    EXPECT_EQ(s.noVars(), prob.no_vars);
    EXPECT_EQ(s.noObjs(), prob.no_objs);
}

TEST(SolutionTests, Equality) {
    TestProblem prob(2, 3);
    AllocSolution s(prob);
    s.setVars({0, 1});
    s.setObjs({2, 3, 4});
    AllocSolution s2(prob);
    s2 = s;
    EXPECT_TRUE(s == s2);
    s.setVars({0.1, 1.1});
    EXPECT_TRUE(s.equals(s2, 0.11));
    EXPECT_FALSE(s.equals(s2, 0.09));
}

TEST(SolutionTests, Dominate) {
    TestProblem prob(1, 2);
    AllocSolution s1(prob);
    AllocSolution s2(prob);
    s1.setObjs({0, 0});
    s2.setObjs({0, 0});
    EXPECT_FALSE(s1 < s2);
    EXPECT_TRUE(s1 <= s2);
    s2.setObjs({0, 1});
    EXPECT_TRUE(s1 < s2);
    EXPECT_TRUE(s1 <= s2);
}

class SolutionTest : public testing::Test
{
protected:
    MOProblem* prob;
    const int no_vars = 2;
    const int no_objs = 2;
    std::vector<double> lb;
    std::vector<double> ub;
    double vars[2];
    double objs[2];
    Solution* sol;

    static void evaluate(void* problem, const double ins[],
                                   double outs[]) {
    };

    virtual void SetUp() {
        lb = {0, 1};
        ub = {1, 2};
        auto ef = [](const double ins[], double outs[]) {
            outs[0] = ins[0] + ins[1];
            outs[1] = ins[0] - ins[1];
        };
        auto vf = [](double x[]) {
            if (x[0] < 0) x[0] = 0;
            if (x[0] > 1) x[0] = 1;
            if (x[1] < 1) x[1] = 1;
            if (x[1] > 2) x[1] = 2;
        };
        prob = new GenericProblem("test", no_vars, no_objs, lb, ub, ef, vf);
        sol = new Solution(prob, vars, objs);
    }

    virtual void TearDown() {
        delete prob;
        delete sol;
    }
};
/*
TEST_F(SolutionTest, MoveConstructor) {
    double vars[2];
    double objs[3];
    Solution s2(prob, vars, objs);
    Solution s(std::move(s2));
    EXPECT_EQ(s.vars, vars);
    EXPECT_EQ(s.objs, objs);
    EXPECT_EQ(s.noVars(), prob->no_vars);
    EXPECT_EQ(s.noObjs(), prob->no_objs);
}
*/
TEST_F(SolutionTest, AssignmentOperator) {
    AllocSolution s2(*prob);
    AllocSolution s(*prob);
    for (int i = 0; i < prob->noVars(); ++i)
        s2.vars[i] = i;
    for (int i = 0; i < prob->noObjs(); ++i)
        s2.objs[i] = i;
    s = s2;
    EXPECT_NE(s.vars, s2.vars);
    EXPECT_NE(s.objs, s2.objs);
    EXPECT_TRUE(allEq(s.vars, s2.vars, prob->noVars()));
    EXPECT_TRUE(allEq(s.objs, s2.objs, prob->noVars()));
    EXPECT_EQ(s.noVars(), prob->no_vars);
    EXPECT_EQ(s.noObjs(), prob->no_objs);
}
/*
TEST_F(SolutionTest, MoveAssignmentOperator) {
    double vars[2];
    double objs[3];
    Solution s2(prob, vars, objs);
    Solution s = std::move(s2);
    EXPECT_EQ(s.vars, vars);
    EXPECT_EQ(s.objs, objs);
    EXPECT_EQ(s.noVars(), prob->no_vars);
    EXPECT_EQ(s.noObjs(), prob->no_objs);
}
*/
TEST_F(SolutionTest, Print) {
    sol->setVars({1, 2});
    sol->setObjs({3, 4});
    std::stringstream ss;
    ss << *sol;
    std::string r = ss.str();
    std::size_t found = r.find("1");
    EXPECT_TRUE(found != std::string::npos);
    found  = r.find("2", found + 1);
    EXPECT_TRUE(found != std::string::npos);
    found  = r.find(";", found + 1);
    EXPECT_TRUE(found != std::string::npos);
    found  = r.find("3", found + 1);
    EXPECT_TRUE(found != std::string::npos);
    found  = r.find("4", found + 1);
    EXPECT_TRUE(found != std::string::npos);
}

TEST_F(SolutionTest, Read) {
    std::stringstream ss;
    ss << "   5.0   -6  ;-7.0 8.  ";
    ss >> *sol;
    EXPECT_DOUBLE_EQ(sol->vars[0],  5.0);
    EXPECT_DOUBLE_EQ(sol->vars[1], -6.0);
    EXPECT_DOUBLE_EQ(sol->objs[0], -7.0);
    EXPECT_DOUBLE_EQ(sol->objs[1],  8.0);
}

TEST_F(SolutionTest, Eval) {
    sol->vars[0] = 0;
    sol->vars[1] = 1;
    sol->evaluate();
    EXPECT_DOUBLE_EQ(sol->objs[0], +1.0);
    EXPECT_DOUBLE_EQ(sol->objs[1], -1.0);
}

TEST_F(SolutionTest, Validate) {
    sol->vars[0] = +10.0;
    sol->vars[1] = -10.0;
    sol->validate();
    EXPECT_DOUBLE_EQ(sol->vars[0], 1.0);
    EXPECT_DOUBLE_EQ(sol->vars[1], 1.0);
}
