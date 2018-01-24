#include <cstdio>
#include <cmath>

#include "helpers.h"
#include "../surrogate/mo.h"
#include "../surrogate/aux.h"

TEST(PopulationTests, Writing) {
    TestProblem prob(3, 1);
    Population pop(2, prob);

    pop[0].setVars({0, 1, 2}); pop[0].setObjs({3});
    pop[1].setVars({4, 5, 6}); pop[1].setObjs({7});

    std::stringstream ss;
    ss << pop;
    std::string r = ss.str();
    std::size_t found = 0;
    found = r.find("2"); EXPECT_TRUE(found != std::string::npos);
    found = r.find("0", found + 1); EXPECT_TRUE(found != std::string::npos);
    found = r.find("1", found + 1); EXPECT_TRUE(found != std::string::npos);
    found = r.find("2", found + 1); EXPECT_TRUE(found != std::string::npos);
    found = r.find(";", found + 1); EXPECT_TRUE(found != std::string::npos);
    found = r.find("3", found + 1); EXPECT_TRUE(found != std::string::npos);
    found = r.find("\n", found + 1); EXPECT_TRUE(found != std::string::npos);
    found = r.find("4", found + 1); EXPECT_TRUE(found != std::string::npos);
    found = r.find("5", found + 1); EXPECT_TRUE(found != std::string::npos);
    found = r.find("6", found + 1); EXPECT_TRUE(found != std::string::npos);
    found = r.find(";", found + 1); EXPECT_TRUE(found != std::string::npos);
    found = r.find("7", found + 1); EXPECT_TRUE(found != std::string::npos);
}

TEST(PopulationTests, Reading) {
    TestProblem prob(3, 1);
    Population pop(1, prob);

    std::stringstream ss(
                "3\n"
                "0 1 2; 3\n"
                "4 5 6; 7\n"
                "8 9 10; 11\n");
    ss >> pop;
    EXPECT_DOUBLE_EQ(pop.size(), 3);
    EXPECT_DOUBLE_EQ(pop[0].vars[0], 0);
    EXPECT_DOUBLE_EQ(pop[0].vars[1], 1);
    EXPECT_DOUBLE_EQ(pop[0].vars[2], 2);
    EXPECT_DOUBLE_EQ(pop[0].objs[0], 3);
    EXPECT_DOUBLE_EQ(pop[1].vars[0], 4);
    EXPECT_DOUBLE_EQ(pop[1].vars[1], 5);
    EXPECT_DOUBLE_EQ(pop[1].vars[2], 6);
    EXPECT_DOUBLE_EQ(pop[1].objs[0], 7);
    EXPECT_DOUBLE_EQ(pop[2].vars[0], 8);
    EXPECT_DOUBLE_EQ(pop[2].vars[1], 9);
    EXPECT_DOUBLE_EQ(pop[2].vars[2], 10);
    EXPECT_DOUBLE_EQ(pop[2].objs[0], 11);
}

TEST(PopulationTests, Serialization) {
    TestProblem prob(3, 1);
    Population popw(2, prob);

    popw[0].setVars({0, 1, 2}); popw[0].setObjs({3});
    popw[1].setVars({4, 5, 6}); popw[1].setObjs({7});

    std::stringstream ss;
    ss << serialize(popw);

    Population pop(1, prob);
    ss >> read_serialized(pop);

    EXPECT_DOUBLE_EQ(pop.size(), 2);
    EXPECT_DOUBLE_EQ(pop[0].vars[0], 0);
    EXPECT_DOUBLE_EQ(pop[0].vars[1], 1);
    EXPECT_DOUBLE_EQ(pop[0].vars[2], 2);
    EXPECT_DOUBLE_EQ(pop[0].objs[0], 3);
    EXPECT_DOUBLE_EQ(pop[1].vars[0], 4);
    EXPECT_DOUBLE_EQ(pop[1].vars[1], 5);
    EXPECT_DOUBLE_EQ(pop[1].vars[2], 6);
    EXPECT_DOUBLE_EQ(pop[1].objs[0], 7);
}

TEST(PopulationTests, Hypervolume) {
    TestProblem prob(1, 2);
    Population pop(5, prob);
    pop.objs() <<
        1, 4,
        3, 4,
        2, 2,
        3, 1,
        6, 1;
    double ref[2] = {5, 5};
    EXPECT_DOUBLE_EQ(pop.hypervolume(ref), 12.0);
}

TEST(PopulationTests, NdSorting) {
    TestProblem prob(1, 2);
    Population pop(6, prob);
    pop.objs() <<
        1, 3,
        2, 3,
        3, 3,
        2, 2,
        2, 3,
        3, 1;
    ematrixd old_objs = pop.objs();
    std::vector<int> fronts;

    pop.nonDominated(fronts);

    ASSERT_EQ(fronts.size(), pop.size());
    EXPECT_EQ(fronts[0], 0);
    EXPECT_EQ(fronts[1], 1);
    EXPECT_EQ(fronts[2], 2);
    EXPECT_EQ(fronts[3], 0);
    EXPECT_EQ(fronts[4], 1);
    EXPECT_EQ(fronts[5], 0);

    pop.nonDominatedSort(fronts);

    ASSERT_EQ(fronts.size(), pop.size());
    EXPECT_EQ(fronts[0], 0);
    EXPECT_EQ(fronts[1], 0);
    EXPECT_EQ(fronts[2], 0);
    EXPECT_EQ(fronts[3], 1);
    EXPECT_EQ(fronts[4], 1);
    EXPECT_EQ(fronts[5], 2);

    // front 0
    EXPECT_TRUE(pop.objs().row(0) == old_objs.row(0) ||
                pop.objs().row(0) == old_objs.row(3) ||
                pop.objs().row(0) == old_objs.row(5));
    EXPECT_TRUE(pop.objs().row(1) == old_objs.row(0) ||
                pop.objs().row(1) == old_objs.row(3) ||
                pop.objs().row(1) == old_objs.row(5));
    EXPECT_TRUE(pop.objs().row(2) == old_objs.row(0) ||
                pop.objs().row(2) == old_objs.row(3) ||
                pop.objs().row(2) == old_objs.row(5));

    // front 1
    EXPECT_TRUE(pop.objs().row(3) == old_objs.row(1) ||
                pop.objs().row(3) == old_objs.row(4));
    EXPECT_TRUE(pop.objs().row(4) == old_objs.row(1) ||
                pop.objs().row(4) == old_objs.row(4));

    // front 2
    EXPECT_TRUE(pop.objs().row(5) == old_objs.row(2));
}

TEST(PopulationTests, Sort) {
    TestProblem prob(1, 1);
    Population pop(4, prob);

    pop[0].setVars({0}); pop[0].setObjs({3});
    pop[1].setVars({1}); pop[1].setObjs({7});
    pop[2].setVars({2}); pop[2].setObjs({1});
    pop[3].setVars({3}); pop[3].setObjs({2});

    int ranks[4] = {3, 7, 1, 2};
    pop.sort(ranks);

    EXPECT_DOUBLE_EQ(pop[0].objs[0], 1);
    EXPECT_DOUBLE_EQ(pop[1].objs[0], 2);
    EXPECT_DOUBLE_EQ(pop[2].objs[0], 3);
    EXPECT_DOUBLE_EQ(pop[3].objs[0], 7);

    EXPECT_DOUBLE_EQ(pop[0].vars[0], 2);
    EXPECT_DOUBLE_EQ(pop[1].vars[0], 3);
    EXPECT_DOUBLE_EQ(pop[2].vars[0], 0);
    EXPECT_DOUBLE_EQ(pop[3].vars[0], 1);
}

TEST(PopulationTests, SortRank) {
    TestProblem prob(1, 1);
    Population pop(4, prob);

    pop[0].setVars({0}); pop[0].setObjs({3});
    pop[1].setVars({1}); pop[1].setObjs({7});
    pop[2].setVars({2}); pop[2].setObjs({1});
    pop[3].setVars({3}); pop[3].setObjs({2});

    pop.sort([](const Solution& a, const Solution& b) {
        return a.objs[0] < b.objs[0];
    });

    EXPECT_DOUBLE_EQ(pop[0].objs[0], 1);
    EXPECT_DOUBLE_EQ(pop[1].objs[0], 2);
    EXPECT_DOUBLE_EQ(pop[2].objs[0], 3);
    EXPECT_DOUBLE_EQ(pop[3].objs[0], 7);

    EXPECT_DOUBLE_EQ(pop[0].vars[0], 2);
    EXPECT_DOUBLE_EQ(pop[1].vars[0], 3);
    EXPECT_DOUBLE_EQ(pop[2].vars[0], 0);
    EXPECT_DOUBLE_EQ(pop[3].vars[0], 1);
}

/*
class PopulationTest : public testing::Test
{
protected:
    MOProblem* prob;
    static const int no_vars = 2;
    static const int no_objs = 2;
    double lb[2];
    double ub[2];

    Solution* pop;
    static const int pop_size = 3;

    static void evaluate(void* problem, const double ins[],
                                   double outs[]) {
        outs[0] = ins[0] + ins[1];
        outs[1] = ins[0] - ins[1];
    };

    virtual void SetUp() {
        lb[0] = 0; lb[1] = 1;
        ub[0] = 1; ub[1] = 2;
        prob = newMOProblem("test", no_vars, no_objs, lb, ub, NULL, evaluate,
                            validateProbTruncating, NULL);
        pop = newSampledPopulation(pop_size, prob, LATIN_HYP_SAMPLING, NULL);
    }

    virtual void TearDown() {
        freeMOProblem(prob);
        freePopulation(pop);
    }
};

TEST_F(PopulationTest, Constructor) {
    for (int i = 0; i < pop_size; ++i)
    {
        EXPECT_EQ(pop[i].problem, prob);
        ASSERT_TRUE(pop[i].vars != NULL);
        EXPECT_TRUE(pop[i].objs != NULL);
        for (int j = 0; j < no_vars; ++j)
            EXPECT_TRUE(isWithin(pop[i].vars[j], prob->lower_bounds[j],
                                 prob->upper_bounds[j]));
        EXPECT_TRUE(pop[i].objs != NULL);
        for (int j = 0; j < no_vars; ++j)
            EXPECT_TRUE(isnan(pop[i].objs[j]));
    }
}

TEST_F(PopulationTest, Print) {
    double v = 0.0;
    for (int i = 0; i < pop_size; ++i)
        for (int j = 0; j < no_vars; ++j)
            pop[i].vars[j] = v++;
    for (int i = 0; i < pop_size; ++i)
        for (int j = 0; j < no_objs; ++j)
            pop[i].objs[j] = v++;
    FILE* f = tmpfile();
    printPopulation(f, pop, pop_size);
    char* f_str = fileToString(f);
    fclose(f);
    EXPECT_TRUE(isSubString(f_str,
                            "0.000000 1.000000; 6.000000 7.000000\n"
                            "2.000000 3.000000; 8.000000 9.000000\n"
                            "4.000000 5.000000; 10.000000 11.000000"));
    free(f_str);
}

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
    free(f_str);
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

TEST_F(PopulationTest, Eval) {
    evaluatePopulation(pop, pop_size);
    double o[no_objs];
    for (int i = 0; i < pop_size; ++i)
    {
        prob->evaluate(prob, pop[i].vars, o);
        EXPECT_TRUE(allDoubleEq(pop[i].objs, o, no_objs));
    }
}

TEST_F(PopulationTest, Validate) {
    double invalid_vars[pop_size][no_vars];
    for (int i = 0; i < pop_size; ++i)
    {
        for (int j = 0; j < no_vars; ++j)
            invalid_vars[i][j] = (RNG::realUniform<double>() - 0.5) * 1000.0;
        memcpy(pop[i].vars, invalid_vars[i], sizeof(double) * no_vars);
    }
    validatePopulation(pop, pop_size);
    double v[no_vars];
    for (int i = 0; i < pop_size; ++i)
    {
        prob->validate(prob, invalid_vars[i]);
        EXPECT_TRUE(allDoubleEq(pop[i].vars, invalid_vars[i], no_vars));
    }
}

TEST_F(PopulationTest, VarsMinMax) {
    //srand(time(NULL));
    for (int i = 0; i < no_vars; ++i)
    {
        int dt[pop_size];
        sampleUniqInt(pop_size * i, pop_size * (i + 1), pop_size, dt);
        for (int j = 0; j < pop_size; ++j)
            pop[j].vars[i] = dt[j];
    }

    //printPopulation(stdout, pop, pop_size);

    double min[no_vars];
    double max[no_vars];

    popMinMaxVars(pop, pop_size, min, max);

    for (int i = 0; i < no_vars; ++i)
    {
        EXPECT_DOUBLE_EQ(min[i], pop_size * i);
        EXPECT_DOUBLE_EQ(max[i], pop_size * (i + 1) - 1.0);
    }
}

TEST_F(PopulationTest, IdealNadir) {
    //srand(time(NULL));
    for (int i = 0; i < no_objs; ++i)
    {
        int dt[pop_size];
        sampleUniqInt(pop_size * i, pop_size * (i + 1), pop_size, dt);
        for (int j = 0; j < pop_size; ++j)
            pop[j].objs[i] = dt[j];
    }

    //printPopulation(stdout, pop, pop_size);

    double ideal[no_vars];
    double nadir[no_vars];

    popIdealNadir(pop, pop_size, ideal, nadir);

    for (int i = 0; i < no_vars; ++i)
    {
        EXPECT_DOUBLE_EQ(ideal[i], pop_size * i);
        EXPECT_DOUBLE_EQ(nadir[i], pop_size * (i + 1) - 1.0);
    }
}

TEST_F(PopulationTest, Normalize) {
    MOProblem* prob = newMOProblem("test", 6, 6, NULL, NULL, NULL, NULL, NULL, NULL);
    Solution* mpop = newPopulationForProblem(3, prob);

    double d[3][6] = {
        {1, 1, 2, 2, 3, 3},
        {2, 3, 1, 3, 1, 2},
        {3, 2, 3, 1, 2, 1},
    };

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 6; ++j)
            mpop[i].vars[j] = d[i][j];
        for (int j = 0; j < 6; ++j)
            mpop[i].objs[j] = d[i][j];
    }

    popNormalizeVars(mpop, 3,  0.0, 1.0);
    popNormalizeObjs(mpop, 3, -1.0, 0.0);

    double dn[3][6] = {
        {0.0, 0.0, 0.5, 0.5, 1.0, 1.0},
        {0.5, 1.0, 0.0, 1.0, 0.0, 0.5},
        {1.0, 0.5, 1.0, 0.0, 0.5, 0.0},
    };

    for (int i = 0; i < pop_size; ++i)
    {
        for (int j = 0; j < no_vars; ++j)
            EXPECT_DOUBLE_EQ(mpop[i].vars[j], dn[i][j]);
        for (int j = 0; j < no_objs; ++j)
            EXPECT_DOUBLE_EQ(mpop[i].objs[j], dn[i][j] - 1.0);
    }

    freePopulation(mpop);
    freeMOProblem(prob);
}

*/
