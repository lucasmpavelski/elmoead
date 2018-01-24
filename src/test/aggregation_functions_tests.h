#include "../surrogate/obj_dec.h"

class AggregationFunctionsTest : public testing::Test {
protected:
    double w[3];
    double x[3];
    double zero[3];
    double ideal[3];
    double nadir[3];

    virtual void SetUp() {
        w[0]     = 1.0; w[1]     = 2.0; w[2]     = 3.0;
        x[0]     = 5.0; x[1]     = 7.0; x[2]     = 4.5;
        zero[0]  = 0.0; zero[1]  = 0.0; zero[2]  = 0.0;
        ideal[0] = 2.0; ideal[1] = 3.0; ideal[2] = 4.5;
        nadir[0] =10.0; nadir[1] = 7.0; nadir[2] = 5.0;
    }
};

TEST_F(AggregationFunctionsTest, WeightedSum) {
    WeightSet<double> weights(w, 1, 3);
    //EXPECT_DOUBLE_EQ(weights.eval(WEIGHTED_SUM, 0, x, ideal, nadir), 32.5);

    EXPECT_DOUBLE_EQ(weightedSum(x, 3, w), 32.5);
    EXPECT_DOUBLE_EQ(aggregationFunction(WEIGHTED_SUM, x, 3, w, NULL, NULL),
                     32.5);
}

TEST_F(AggregationFunctionsTest, Tchebycheff) {
    WeightSet<double> weights(w, 1, 3);
    //EXPECT_DOUBLE_EQ(weights.eval(TCHEBYCHEFF, 0, x, zero, nadir), 14.0);

    EXPECT_DOUBLE_EQ(tchebycheff(x, 3, w, zero), 14.0);
    EXPECT_DOUBLE_EQ(aggregationFunction(TCHEBYCHEFF, x, 3, w, zero, NULL),
                     14.0);
}

TEST_F(AggregationFunctionsTest, AugmentedTchebycheff) {
    WeightSet<double> weights(w, 1, 3);
    //EXPECT_DOUBLE_EQ(weights.eval(AUGMENTED_TCHEBYCHEFF, 0, x, zero, nadir), 14.325);

    /*
    EXPECT_DOUBLE_EQ(augmentedTchebycheff(x, 3, w, zero), 14.325);
    EXPECT_DOUBLE_EQ(aggregationFunction(AUGMENTED_TCHEBYCHEFF, x, 3, w, zero,
                                         NULL), 14.325);
    */

    EXPECT_DOUBLE_EQ(augmentedTchebycheff(x, 3, w, zero), 5.325);
    EXPECT_DOUBLE_EQ(aggregationFunction(AUGMENTED_TCHEBYCHEFF, x, 3, w, zero,
                                         NULL), 5.325);
}

TEST(AggregationFunctionTests, NormalizedTchebycheff) {
    const int no_objs = 3.0;
    const double ideal[] = {0, 1, 2};
    const double nadir[] = {2, 4, 6};
    const double w[] = {1, 2, 3};
    const double x[] = {1, 2, 4};
    AggrFunc<NORMALIZED_TCHEBYCHEFF> af;
    EXPECT_DOUBLE_EQ(af.apply(x, 3, w, ideal, nadir), 1.5);
}

TEST_F(AggregationFunctionsTest, PenaltyBoundaryIntersection) {
  AggrFunc<PENALTY_BOUNDARY_INTERSECTION> af(5.0);
  EXPECT_DOUBLE_EQ(af.apply(x, 3, w, ideal, nadir), 48.029445475701678);
}
