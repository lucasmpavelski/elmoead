#include "surrogate/configurable_test.h"

int main(int argc, char *argv[])
{
    using std::unique_ptr;

    ConfigurableTest ct(argc, argv);

    auto prob     = unique_ptr<MOProblem>{ct.getProblem()};
    auto moead    = unique_ptr<MOEAD>{ct.getMOEAD(prob->noVars(), prob->noObjs())};
    auto callback = ct.getCallback();
    auto no_evals = ct.getl("no_evals", 200);

    Population pop(moead->initPopSize(), *prob, LATIN_HYP_SAMPLING);
    pop.evaluateAll();
    moead->optimize(no_evals - moead->initPopSize(), *prob, pop, callback);


    size_t no_nd = pop.nonDominatedSort();
    pop.resize(no_nd);

    /* write result in file */
    auto front_fn = ct.front_fn;
    std::ofstream final_pop_file(front_fn.c_str());
    final_pop_file << pop.objs();
    final_pop_file.close();

    const auto hyp_ref = ct.gets("hyp_ref", "5.0 5.0");
    auto hyp_ref_vec = str2Vector<double>(hyp_ref);
    throw_assert(hyp_ref_vec.size() == prob->noObjs(),
        "hyp_ref dimension should be igual to the number of objectives");

    /* print time and hipervolume */
    fmt::print("hipervolume: {}\ntime: {}\nfront: {}",
        pop.hypervolume(hyp_ref_vec.data()), Clock::tellTime().count(), front_fn);

    return EXIT_SUCCESS;
}

//    Clock::begin();
//
//    cout << "seed: " << RNG::seed(//2150554808//8790;//3103;////18008;//3103;//
//                                  ) << endl;
//
//
//    const char* final_pop_fn = (argc > 1) ? argv[1] : "final_pop.dat";
////    const char* inner_pop_fn = (argc > 2) ? argv[2] : "inner_pop.dat";
//
//    // problem params
///*
//
//    initTkinter();
//    const char* prob_name = (argc > 3) ? argv[3] : "1PLWm7";
//    Proteinpp protein = buildProtein(prob_name);
//    PSPProblempp pprob(protein);
//    MOProblem& prob = pprob;
//
////    const char* prob_name = "adp";
////    std::unique_ptr<MOProblem> prob(new AirfoilDesignProblem());
//
//
//
///*/
//
//    const char* prob_name = (argc > 3) ? argv[3] : "WFG1";
//    const int no_prob_vars = (argc > 4) ? strtol(argv[4], NULL, 10) : 24;
//    std::unique_ptr<MOProblem> prob{getBenchmark(prob_name, no_prob_vars)};
////*/
//    const int no_vars = prob->noVars();
//    const int no_objs = prob->noObjs();
//    const double hv_ref[] = {5.0, 5.0}; // hipervolume reference
//
//    int no_evals = 5000;//no_vars * 200;
////         if (no_vars ==  8)                  no_evals =  5000;
////    else if (no_vars == 30 || no_vars == 10) no_evals = 25000;
////    else if (no_vars == 60 || no_vars == 20) no_evals = 50000;
//         if (no_vars ==  8)                  no_evals =  200;
//    else if (no_vars == 30 || no_vars == 10) no_evals = 1000;
//    else if (no_vars == 60 || no_vars == 20) no_evals = 2000;
//
//    // moea/d params
//    const int no_partitions = 99;
//    const int no_neighbors = 20;
//    const AggregationFunctionType aggr_func = PENALTY_BOUNDARY_INTERSECTION;
//    const PMOperator pmo(no_vars);
//    const SBXoverOperator sbxo(1.0, 20.0);
//    PopulationCallback cb =
////            NullCallback;
//                                 PlotCallback(5);
//                                 //LogCallback(final_pop_fn);
//
//    // moead-de params
//    const double real_b = 0.9;
//    const int update_limit = 1;
//    std::vector<DEOperator> de_ops = {
//        DEOperator(DEType::RAND_1_BIN, 1.0, 0.5),
//      //  DEOperator(DEType::RAND_2_BIN, 1.0, 0.5),
//      //  DEOperator(DEType::NON_LINEAR),
//    };
//
//    OperatorSelection<DEOperator> operators(de_ops); // const operator
//    //ProbabilityMatching<DEOperator> operators(de_ops, ProbabilityMatching<DEOperator>::AvgAbs);
//    //FRRMAB<DEOperator> operators(de_ops, noWeightsFor(no_objs, no_partitions));
//
//    WeightSet<double> weights(no_objs, no_partitions, WeightSet<>::Type::UNIFORM);
////    MOEAD moead(weights, no_neighbors, aggr_func, pmo, sbxo, cb);
//    MOEAD_DE moead(weights, no_neighbors, aggr_func, real_b, update_limit, pmo, operators, cb);
//    //MOEAD_DRA moead(weights, no_neighbors, aggr_func, real_b, update_limit, pmo, operators, cb);
//    //MOEAD_STM moead(weights, no_neighbors, aggr_func, real_b, pmo, operators, cb);
//
//    Population pop(moead.initPopSize(), *prob, LATIN_HYP_SAMPLING);
//    pop.evaluateAll();
//
//    moead.optimize(no_evals - pop.size(), *prob, pop);
//
//    const double hipervolume = pop.hypervolume(hv_ref);
//
//    cout << "hypervolume: " << std::setprecision(10) << hipervolume << "\n"
//         << "time: " << Clock::tellTime().count() << endl;
//
//    int no_nd = pop.nonDominatedSort();
//    pop.resize(no_nd);
//
////    Plotter plt;
////    pop.plotObjs(plt);
//
//    // write result in file
//    if (final_pop_fn)
//    {
//        std::ofstream final_pop_file(final_pop_fn);
//        final_pop_file << pop.objs() << endl;
//        final_pop_file.close();
//    }
//
//    //interpret(final_pop_fn, *protein, pop.solsData(), no_nd);
//
//    return 0;
//}


/*
FILE* pp;
FILE* moead_ipf;
Solution* pop_cp;

void moeadCallbackPop(const int gen, MOEADConfig* moead,
                      const MOProblem* prob, const Solution* pop)
{
    const int pop_size = moead->pop_size;

    int i;
    //for (i = 0; i < pop_size; ++i)
    //    copySolution(&pop[i], &pop_cp[i]);
    //const int no_nd = nonDominatedSort(pop_cp, pop_size, NULL);
    //plotPopulationCmd(pp, pop_cp, no_nd);

    plotPopulationCmd(pp, pop, pop_size);

    printf("gen: %d\n", gen);
    fflush(stdout);
}

int main(int argc, char** argv)
{
    int i;

    const long seed = (argc >= 2) ? strtol(argv[1], NULL, 10) : time(NULL);
    srand(seed);

    const char* final_pop_fn  = (argc >= 3) ? argv[2] : "final_pop.dat";
    const char* inner_pops_fn = (argc >= 4) ? argv[3] : "inner_pop.dat";

    const char* prob_name = "UF1";
    const int no_prob_vars = 30;
    MOProblem* prob = getBenchmark(prob_name, no_prob_vars);

    const int no_vars = prob->no_vars;
    const int no_objs = prob->no_objs;
    const double hv_ref[] = {5.0, 5.0}; // hipervolume reference

    pp = initPlotPopulationProc();

    MOEADConfig moead;
    moead.max_aval         = 1000              ;
    moead.partitions       = 299                 ;
    moead.no_neighbors     = 20                  ;
    moead.real_b           = 0.90                ;
    moead.update_limit     = 10                  ;
    moead.sampling_method  = LATIN_HYP_SAMPLING  ;
    moead.aggregation_type = INVERTED_TCHEBYCHEFF;
    moead.de_operator      = DEOperator::Type::DE_1_BIN;
    moead.var_adapt        = false               ;
    moead.de_cr            = 1.0                 ;
    moead.de_f             = 0.5                 ;
    moead.selecao_reward   = 1                   ;

    initMOEAD(prob, &moead);
    const int pop_size = moead.pop_size;

    Solution* moead_pop = initMOEADSTMPopulation(&moead, prob, NULL);
    pop_cp = newPopulationForProblem(pop_size, prob);

    moead_ipf = fopen(inner_pops_fn, "w");
    STMData* stm_data = newSTMData(pop_size * 2, pop_size);
    //MOEADSTMoptimize(&moead, prob, moead_pop, moeadCallbackPop, stm_data);
    MOEADDEoptimize(&moead, prob, moead_pop, moeadCallbackPop, stm_data);
    fclose(moead_ipf);

    const int no_nd = nonDominatedSort(moead_pop, pop_size, NULL);

    FILE* final_pop_file = fopen(final_pop_fn, "w");
    printPopulationObjs(final_pop_file, moead_pop, no_nd);
    fclose(final_pop_file);

    deinitPlotPopulationProc(&pp);

    freePopulation(moead_pop);
    freePopulation(pop_cp);
    delete prob;
    deinitMOEAD(&moead);
    freeSTMData(stm_data);

    return EXIT_SUCCESS;
}*/
