#include "helpers.h"
#include "../surrogate/aux.h"
#include "../surrogate/aux/bounds.h"

TEST(RNGTests, generation) {
    RNG::seed(time(NULL));
    const int nrolls = 10000; // number of experiments
    const int nstars = 95;    // maximum number of stars to distribute

    int pi[10] = {}, pr[10] = {};
    int count = 0;  // count number of trues

    for (int i = 0; i < nrolls; ++i)
    {
        int ri = RNG::intUniform(2, 7);
        pi[ri]++;
        int rr = RNG::realUniform(2.0, 8.0);
        pr[int(rr)]++;
        if (RNG::flipCoin())
            ++count;
    }

    std::cout << "uniform_int_distribution [2,7]:" << std::endl;
    for (int i = 0; i < 10; ++i)
        std::cout << i << ": " << std::string(pi[i]*nstars/nrolls,'*') << std::endl;

    std::cout << "uniform_real_distribution [2,8):" << std::endl;
    for (int i = 0; i < 10; ++i)
        std::cout << i << ": " << std::string(pr[i]*nstars/nrolls,'*') << std::endl;

    std::cout << "bernoulli_distribution (0.5):" << std::endl;
    std::cout << "true:  " << count << std::endl;
    std::cout << "false: " << nrolls-count << std::endl;

}

TEST(AuxTests, Intpow) {
    EXPECT_DOUBLE_EQ(intpow(5.0, 2), 25.0);
    EXPECT_DOUBLE_EQ(intpow(3.0, 3), 27.0);
    EXPECT_DOUBLE_EQ(intpow(3.141592, 2), 9.8696002944640018);
    EXPECT_NEAR(intpow(0.2, 8), 0.00000256, 0.0000001);
}

TEST(AuxTests, Binomial) {
    EXPECT_EQ(binomial(0, 0), 1);
    EXPECT_EQ(binomial(0, 1), 1);
    EXPECT_EQ(binomial(1, 0), 1);
    EXPECT_EQ(binomial(2, 1), 2);
    EXPECT_EQ(binomial(2, 2), 1);
    EXPECT_EQ(binomial(3, 2), 3);
    EXPECT_EQ(binomial(5, 3), 10);
}

TEST(AuxTests, ranks) {
    double ins[6][3] = {
        {0, 1, 2},
        {0, 2, 1},
        {1, 0, 2},
        {1, 2, 0},
        {2, 0, 1},
        {2, 1, 0}
    };
    int gos[6][3] = {
        {0, 1, 2},
        {0, 2, 1},
        {1, 0, 2},
        {2, 0, 1},
        {1, 2, 0},
        {2, 1, 0}
    };
    for (int i = 0; i < 6; ++i)
    {
        int* o = ranks(ins[i], 3, NULL);

        EXPECT_TRUE(allEq(o, gos[i], 3));
        free(o);
    }
}

TEST(AuxTests, VectorSum) {
    double v[3] = {-1, 4, 3};
    EXPECT_DOUBLE_EQ(sum(v, 3), 6.0);
}

TEST(AuxTests, VectorMean) {
    double v[8] = {2,4,4,4,5,5,7,9};
    EXPECT_DOUBLE_EQ(mean(v, 8), 5.0);
}

TEST(AuxTests, VectorNorm) {
    double v[3] = {0, 4, 3};
    EXPECT_DOUBLE_EQ(norm(v, 3), 5.0);
}

TEST(AuxTests, VectorProjection) {
    double v1[2] = {4, 3};
    double v2[2] = {5, 0};
    double gp[2] = {4, 0};
    double* ep = projection(v1, v2, 2, NULL);
    EXPECT_TRUE(allDoubleEq(ep, gp, 2));
    free(ep);
}

TEST(AuxTests, MatrixColsMean) {
    double m1[3] = {1, 2,  3};
    double m2[3] = {3, 2, -1};
    double* m[2] = {m1, m2};
    double cm[3];
    colwiseMean(m, 2, 3, cm);
    EXPECT_DOUBLE_EQ(cm[0], 2.0);
    EXPECT_DOUBLE_EQ(cm[1], 2.0);
    EXPECT_DOUBLE_EQ(cm[2], 1.0);
}

TEST(AuxTests, VectorStd) {
    double v[8] = {2,4,4,4,5,5,7,9};
    EXPECT_DOUBLE_EQ(stdWithMean(v, 8, 5.0), 2.0);
    EXPECT_DOUBLE_EQ(stddev(v, 8), 2.0);
}

TEST(AuxTests, MatrixColsStd) {
    double m1[3] = {1, 2,  3};
    double m2[3] = {3, 2, -1};
    double* m[2] = {m1, m2};
    double cm[3];
    double cstd[3];
    colwiseMeanStd(m, 2, 3, cm, cstd);
    EXPECT_DOUBLE_EQ(cstd[0], 1.0);
    EXPECT_DOUBLE_EQ(cstd[1], 0.0);
    EXPECT_DOUBLE_EQ(cstd[2], 2.0);
}

TEST(AuxTests, Serialization) {
    std::vector<double> v = {10e-6, M_PI, 1.1};
    std::stringstream ss;
    for (auto& el : v)
    {
        ss << serialize(el);
        double r;
        ss >> read_serialized(r);
        EXPECT_DOUBLE_EQ(el, r);
    }
}

TEST(AuxTests, SerializationRvalue) {
    std::stringstream ss;
    ss << serialize<double>(2.14 + 1.0);
    double r;
    ss >> read_serialized(r);
    EXPECT_DOUBLE_EQ(r, 3.14);
}

TEST(AuxTests, CircularRange) {
    std::vector<int> v = {1, 2, 3, 5, 8};
    auto r = circularRange(v.begin(), v.end(), 0, 30);
    EXPECT_EQ(r.first, 1);
    EXPECT_EQ(r.second, 8);
    r = circularRange(v.begin(), v.end(), 0, 9);
    EXPECT_EQ(r.first, 8);
    EXPECT_EQ(r.second, 5);
}
/*
TEST(AuxTests, MakeBounds) {
    svector<double> v1 = {1, 2};
    svector<double> v2 = {4, 5, 6};
    auto b = makeBounds(v1, v2);
    ASSERT_EQ(b.size(), 2);
    EXPECT_EQ(b[0].first, 1);
    EXPECT_EQ(b[0].second, 4);
    EXPECT_EQ(b[1].first, 2);
    EXPECT_EQ(b[1].second, 5);
}
*/
TEST(AuxTests, MakeBoundsList) {
    auto b = makeBounds({1, 2}, {4, 5, 6});
    ASSERT_EQ(b.size(), 2);
    EXPECT_EQ(b[0].first, 1);
    EXPECT_EQ(b[0].second, 4);
    EXPECT_EQ(b[1].first, 2);
    EXPECT_EQ(b[1].second, 5);
}

TEST(AuxTests, BoundsIsBetween) {
    auto b = makeBounds<double>({1, 2}, {4, 5});
    EXPECT_TRUE(isBetween(b, {1, 5}));
    EXPECT_FALSE(isBetween(b, {0, 3}));
    EXPECT_FALSE(isBetween(b, {2, 6}));
}

TEST(AuxTests, BoundsTruncate) {
    auto b = makeBounds<double>({1, 2}, {4, 5});
    svectord v = {0, 6};
    truncate(b, v.begin());
    EXPECT_DOUBLE_EQ(v[0], 1.0);
    EXPECT_DOUBLE_EQ(v[1], 5.0);
}
