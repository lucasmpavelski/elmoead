#include "../surrogate/train_archive.h"
#include "../surrogate/mo.h"

#include "helpers.h"

class TrainArchiveTests : public testing::Test
{
protected:
    Solution* pop;
    MOProblem* prob;

    virtual void SetUp() {
        prob = new TestProblem(2, 2);
    };

    virtual void TearDown() {
        delete prob;
    };
};

// TEST_F(TrainArchiveTests, AddToTrainArchive)
// {
//     WeightSet<double> w(2, 4, WeightSet<>::Type::UNIFORM, 0);
//     Population pop(5, *prob);

//     double** ins = newMatrix<double>(5, 2);
//     double** outs = newMatrix<double>(5, 2);

//     for (int i = 0; i < 5; ++i)
//     {
//         auto wi = w[i];
//         for (int j = 0; j < 2; ++j)
//             pop[i].objs[j] = wi[j] * 1.6;
//     }

//     cout << w << endl;
//     TrainingArchive<double> ta(WeightSet<double>(2, 4), INVERTED_TCHEBYCHEFF,
//                                std::move(pop));

//     //plotPopulation(ta.solutions, ta.size);
//     //printPopulationObjs(stdout, ta.solutions, ta.size);
//     //printf("\n");

//     for (int r = 0; r < 5; ++r)
//     {
//         for (int j = 0; j < 2; ++j)
//         {
//             auto wr = w[r];
//             outs[r][j] = wr[j] * 1.1;
//         }

//         Solution sol(prob, ins[r], outs[r]);
//         const Solution* p[1] = {&sol};
//         ta.add(p, 1);

//         EXPECT_TRUE(allDoubleEq(ta.solutions[r].objs, sol.objs, 2));

//         //plotPopulation(ta.solutions, ta.size);
//         //printPopulationObjs(stdout, ta.solutions, ta.size);
//         //printf("\n");
//     }

//     freeMatrix(ins);
//     freeMatrix(outs);
// }

TEST_F(TrainArchiveTests, AddToTrainArchiveFull) {

    double objs[4][2] = {
        { 0.0,  1.1}, // enter
        { 1.1,  0.0}, // enter
        {0.875,  0.875}, // dont enter
        { 0.75,  0.75}, // dont enter
    };
    double in[2] = {0};

    //double* ins[4]  = {in,in,in,in};
    //double* outs[4] = {objs[0],objs[1],objs[2],objs[3]};

    Population pop(4, *prob);
    for (int i = 0; i < 4; ++i)
        pop[i].setObjs(objs[i]);

    TrainingArchive<double> ta(WeightSet<double>(2,3), INVERTED_TCHEBYCHEFF,
                               std::move(pop));

    double new_out[2] = {1, 1};
    Solution s;
    s.problem = prob;
    s.vars = in;
    s.objs = new_out;
    ta.add(s);

    // new_out should replace {0.875,  0.875}, not { 0.75,  0.75}
    double replaced[2] = {0.875,  0.875};

    for (int r = 0; r < 4; ++r)
        EXPECT_FALSE(allDoubleEq(ta.solutions[r].objs, replaced, 2));
}
