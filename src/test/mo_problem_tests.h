#include "helpers.h"
#include "../surrogate/mo.h"
#include "../surrogate/aux.h"

TEST(MOProblemTests, ValidationByTruncating) {
    class ProbEx : public TestProblem
    {
    public:
        using TestProblem::TestProblem;
        void validate(double* x) final override {
            truncateToBounds(x);
        };
    };

    ProbEx mop(5, 0, -180.0, 180.0);
    double x[10] = {0.0, -190.0, 190.0, -180-360-10, 180+360+10};
    mop.validate(x);
    double vx[10] = {0.0, -180.0, +180.0, -180.0, +180.0};
    EXPECT_TRUE(allDoubleEq(x, vx, 5));
}

TEST(MOProblemTests, ValidationByReflecting) {
    class ProbEx : public TestProblem
    {
    public:
        using TestProblem::TestProblem;
        void validate(double* x) final override {
            reflectToBounds(x);
        };
    };

    ProbEx mop(5, 0, -180.0, 180.0);
    double x[10] = {0.0, -190.0, 190.0, -180-360-10, 180+360+10};
    mop.validate(x);
    double vx[10] = {0.0, 170.0, -170.0, 170.0, -170.0};
    EXPECT_TRUE(allDoubleEq(x, vx, 5));
}

/*
class MOProblemTest : public testing::Test
{
protected:
    MOProblem* prob;
    int no_vars;
    int no_objs;
    double lb[2];
    double ub[2];
    void* info;

    static void evaluate(void* problem, const double ins[], double outs[]) {
        outs[0] = ins[0] + ins[1];
        outs[1] = ins[0] - ins[1];
    };

    static void validate(void* problem, double inputs[]) {
        if (inputs[0] <= 0) inputs[0] = 1;
        if (inputs[0] >= 5) inputs[0] = 4;
        if (inputs[1] <= 5) inputs[1] = 6;
        if (inputs[1] >= 9) inputs[1] = 5;
    };

    virtual void SetUp() {
        no_vars = no_objs = 2;
        info = (void*)'T';
        lb[0] = 0; lb[1] = 1;
        ub[0] = 1; ub[1] = 2;
        prob = TestProblem(no_vars, no_objs, lb, ub, info, evaluate,
                            validate, NULL);
    };

    virtual void TearDown() {
        freeMOProblem(prob);
    };
};

TEST_F(MOProblemTest, Constructor) {
    EXPECT_STRCASEEQ(prob->name, "test");
    EXPECT_EQ(prob->no_vars, no_vars);
    EXPECT_EQ(prob->no_objs, no_objs);
    EXPECT_TRUE(allDoubleEq(prob->lower_bounds, lb, no_vars));
    EXPECT_TRUE(allDoubleEq(prob->upper_bounds, ub, no_objs));
}

TEST_F(MOProblemTest, Print) {
    FILE* f = fopen("/tmp/mo_print", "w+");
    printMOProblem(f, prob);
    char* f_str = fileToString(f);
    fclose(f);
    EXPECT_TRUE(isSubString(f_str, "name: \"test\""));
    EXPECT_TRUE(isSubString(f_str, "no_vars: 2"));
    EXPECT_TRUE(isSubString(f_str, "no_objs: 2"));
    EXPECT_TRUE(isSubString(f_str, "lower_bounds: [0.000000 1.000000]"));
    EXPECT_TRUE(isSubString(f_str, "upper_bounds: [1.000000 2.000000]"));
    free(f_str);
}

TEST_F(MOProblemTest, Eval) {
    double x[2] = {1.0, 2.0};
    double y[2];
    evaluate(prob, x, y);
    double eval_y[2];
    prob->evaluate(prob, x, eval_y);
    EXPECT_TRUE(allDoubleEq(y, eval_y, no_objs));
}

TEST_F(MOProblemTest, Validate) {
    double x[2] = {10.0, -10.0};
    validate(prob, x);
    double val_x[2] = {10.0, -10.0};
    prob->validate(prob, val_x);
    EXPECT_TRUE(allDoubleEq(x, val_x, no_objs));
}

// TODO: test pre-defined validateions (truncate + reflect)
*/
