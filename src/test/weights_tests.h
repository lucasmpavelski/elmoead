#include "helpers.h"
#include "../surrogate/obj_dec.h"

#include <iostream>
using std::cout;
using std::endl;

TEST(WeightsTests, NoWeightsFor) {
    EXPECT_EQ(WeightSet<>::sizeFor(2, 2), 3);
    EXPECT_EQ(WeightSet<>::sizeFor(3, 2), 6);
    EXPECT_EQ(WeightSet<>::sizeFor(3, 3), 10);
}

TEST(WeightsTests, UniformWeights2D) {
    const int dim = 2;
    const int no_partitions = 2;
    const int no_weights = WeightSet<>::sizeFor(dim, no_partitions);
    double given_w[3][2] = {{0, 1}, {.5, .5}, {1, 0}};

    WeightSet<double> w(dim, no_partitions, WeightSet<>::Type::UNIFORM, 0);

    for (int i = 0; i < no_weights; i++)
        EXPECT_TRUE(allDoubleEq(w[i], given_w[i], dim));
}

#include <iostream>

TEST(WeightsTests, UniformWeights3D) {
    const int dim = 3;
    const int no_partitions = 2;
    const int no_weights = WeightSet<>::sizeFor(dim, no_partitions);
    double z = 0.001;
    double given_w[6][3] = {
        { z,  z,  1},
        { z, .5, .5},
        { z,  1,  z},
        {.5,  z, .5},
        {.5, .5,  z},
        { 1,  z,  z}
    };

    WeightSet<double> w(dim, no_partitions, WeightSet<>::Type::UNIFORM, z);
    //std::cout << w << std::endl;

    for (int i = 0; i < no_weights; i++)
        EXPECT_TRUE(allDoubleEq(w[i], given_w[i], dim));
}

TEST(WeightsTests, RandomWeights) {
    const int dim = 10;
    const int no_partitions = 3;
    auto weights = WeightSet<>(dim, no_partitions,
                                       WeightSet<>::Type::RANDOM);

    for (const auto& w : weights)
        EXPECT_NEAR(std::accumulate(w, w + 10, 0.0), 1.0, 0.001);
}

TEST(WeightsTests, SelfNeighborhood) {
    const int no_weights = 5;
    const int dim = 2;
    const int no_neighbors = 3;

    double all_w[no_weights * dim] = {
        1.0, 1.0,
        1.0, 3.0,
        3.0, 2.0,
        6.0, 3.0,
        6.0, 1.0,
    };

    WeightSet<double> w(all_w, no_weights, dim);
    Neighborhood B(w, no_neighbors);

    int given_b[no_weights][no_neighbors] = {
        {0, 1, 2},
        {1, 0, 2},
        {2, 1, 0},
        {3, 4, 2},
        {4, 3, 2},
    };

    for (int i = 0; i < no_weights; i++)
        EXPECT_TRUE(containAll(B[i], no_neighbors, given_b[i], no_neighbors));
}

TEST(WeightsTests, Neighborhood) {
    const int no_weightsa = 5;
    const int dima = 2;
    double gwa[no_weightsa * dima] = {
        0, 1,
        0, 0,
        1, 0,
    };

    WeightSet<double> wa(gwa, no_weightsa, dima);

    const int no_weightsb = 6;
    const int dimb = 2;
    double gwb[no_weightsb * dimb] = {
         0.25,  0.25, // ~> wa[1]
         0.15,  1.50, // ~> wa[0]
         0.85, -0.15, // ~> wa[2]
        -0.25, -0.25, // ~> wa[1]
         1.25,  0.25, // ~> wa[2]
        -0.20,  0.80, // ~> wa[0]
    };

    WeightSet<double> wb(gwb, no_weightsb, dimb);

    const int t = 2;
    Neighborhood B(wa, wb, t);

    int given_b[no_weightsa][t] = {
        {1, 5},
        {0, 3},
        {2, 4},
    };

    for (int i = 0; i < t; i++)
        EXPECT_TRUE(containAll(B[i], dima, given_b[i], dima));
}
