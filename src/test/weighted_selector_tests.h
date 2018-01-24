#include "../surrogate/weighted_selector.h"

class WeightedSelectorTest : public testing::Test
{
protected:
    double** in_w;
    double** out_w;
    Solution* pop;

    virtual void SetUp() {
        in_w = NULL;
        out_w = NULL;
        pop = NULL;
   //     prob = NULL;
  //      sol.vars = vars;
  //      sol.objs = objs;
    };

    virtual void TearDown() {
    };
};

TEST_F(WeightedSelectorTest, Select)
{
    /*in_w = uniformWeights(2, 5);
    out_w = uniformWeights(2, 2);*/

    TestProblem prob(2, 2);
    Population pop(6, prob);

    double pop_objs[6][2] = {
        { 0, 55},
        {12, 48},
        {24, 36},
        {33, 22},
        {44, 11},
        {60,  0},
    };

    for (int i = 0; i < 6; ++i)
        pop[i].setObjs(pop_objs[i]);
    //plotPopulation(pop, 6);

    WeightSet<double> in_w(2, 5);
    WeightSet<double> out_w(2, 2);
    WeightedSelector<double> ws(in_w, out_w);

    AggrProblem<double> ap(INVERTED_TCHEBYCHEFF, 2);
    const double ideal[] = {0, 0};
    const double nadir[] = {100, 100};
    ap.updateReferencePoints(ideal);
    ap.updateReferencePoints(nadir);

    std::vector<int> ranks(6);

    ws.rankSolutions(in_w, pop, ap, ranks);

   /* WeightedSelector* ws = newWeightedSelector(in_w, 6, out_w, 3, 2, NULL);
    //printWeightedSelector(stdout, ws);

    const double ideal[] = {0, 0};
    const double nadir[] = {100, 100};
    int ranks[6];
    rankSolutions(ws, pop, INVERTED_TCHEBYCHEFF, ideal, nadir, ranks);*/

    int granks[] = {0, 3, 4, 1, 2, 5};
    for (int i = 0; i < 6; ++i) {
        EXPECT_EQ(ranks[i], granks[i]);
    }

/*
    EXPECT_TRUE(allEq(ranks, granks, 6));

    freeWeightedSelector(ws);


    freeMatrix(in_w);
    freeMatrix(out_w);*/
}
