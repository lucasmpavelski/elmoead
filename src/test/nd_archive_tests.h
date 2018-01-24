#include "../surrogate/mo/nd_archive.h"
#include "helpers.h"
#include <algorithm>
using std::any_of;

AssertionResult popIncludeSol(const Solution* pop, const int pop_size,
                              const Solution* sol)
{
    const int no_vars = sol->problem->no_vars;
    const int no_objs = sol->problem->no_objs;

    bool r = any_of(pop, pop + pop_size, [=](const Solution& s) {
        return allDoubleEq(sol->vars, s.vars, no_vars) &&
               allDoubleEq(sol->objs, s.objs, no_objs);
    });

    return r ? AssertionSuccess() : AssertionFailure();
}
/*
TEST(NdArchiveTests, All)
{
    const int pop_size = 7;
    double objs[pop_size][2] = {
        {1, 5},
        {3, 3},
        {5, 1},
        {2, 2},
        {6, 2},
        {0, 4},
        {6, 0},
    };

    TestProblem prob(2, 2);
    Population pop(7, prob, RANDOM_SAMPLING);
    Population nd_pop (7, prob, RANDOM_SAMPLING);

    for (int i = 0; i < pop_size; ++i)
        pop[i].setObjs(objs[i]);

    //NdArchive* nda = newNdArchive(prob, pop_size, 0, NULL);
    MinDistanceArchive nda(prob, pop_size, 0);

    nda.add(pop[0]);
    nda.add(pop[1]);
    nda.add(pop[2]);
    int no_nd = nda.getNdPopulation(nd_pop);
    EXPECT_EQ(no_nd, 3);
    EXPECT_TRUE(popIncludeSol(nd_pop, no_nd, &pop[0]));
    EXPECT_TRUE(popIncludeSol(nd_pop, no_nd, &pop[1]));
    EXPECT_TRUE(popIncludeSol(nd_pop, no_nd, &pop[2]));

    nda.add(pop[3]);
    no_nd = nda.getNdPopulation(nd_pop);
    EXPECT_EQ(no_nd, 3);
    EXPECT_TRUE(popIncludeSol(nd_pop, no_nd, &pop[0]));
    EXPECT_TRUE(popIncludeSol(nd_pop, no_nd, &pop[2]));
    EXPECT_TRUE(popIncludeSol(nd_pop, no_nd, &pop[3]));

    nda.add(pop[4]);

    no_nd = nda.getNdPopulation(nd_pop);
    EXPECT_EQ(no_nd, 3);
    EXPECT_TRUE(popIncludeSol(nd_pop, no_nd, &pop[0]));
    EXPECT_TRUE(popIncludeSol(nd_pop, no_nd, &pop[2]));
    EXPECT_TRUE(popIncludeSol(nd_pop, no_nd, &pop[3]));

    nda.add(pop[5]);

    no_nd = nda.getNdPopulation(nd_pop);
    EXPECT_EQ(no_nd, 3);
    EXPECT_TRUE(popIncludeSol(nd_pop, no_nd, &pop[2]));
    EXPECT_TRUE(popIncludeSol(nd_pop, no_nd, &pop[3]));
    EXPECT_TRUE(popIncludeSol(nd_pop, no_nd, &pop[5]));

    nda.add(pop[6]);

    no_nd = nda.getNdPopulation(nd_pop);
    EXPECT_EQ(no_nd, 4);
    EXPECT_TRUE(popIncludeSol(nd_pop, no_nd, &pop[6]));
    EXPECT_TRUE(popIncludeSol(nd_pop, no_nd, &pop[2]));
    EXPECT_TRUE(popIncludeSol(nd_pop, no_nd, &pop[3]));
    EXPECT_TRUE(popIncludeSol(nd_pop, no_nd, &pop[5]));

    freePopulation(pop);
    freePopulation(nd_pop);
    delete prob;
}*/
