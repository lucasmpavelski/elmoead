#include "surrogate/configurable_test.h"

int main(int argc, char* argv[])
{
    using fmt::print;
    using std::cout;
    using std::unique_ptr;

    auto ct = ConfigurableTest(argc, argv);

    auto prob     = unique_ptr<MOProblem>{ct.getProblem()};
    auto moead    = unique_ptr<MOEAD>{ct.getMOEAD(prob->noVars(), prob->noObjs())};
    auto elmoead  = ct.getELMOEAD(*prob, *moead);
    auto elms     = ct.getELMs(*prob);
    auto callback = ct.getCallback();

    auto hyp_ref = ct.gets("hyp_ref", "5.0 5.0");
    auto hyp_ref_vec = str2Vector<double>(hyp_ref);
    throw_assert(hyp_ref_vec.size() == prob->noObjs(),
        "hyp_ref dimension should be igual to the number of objectives");

    auto moead_no_generations = ct.getl("moead_no_generations", 100);
    const size_t no_moead_eval = moead_no_generations * moead->popSize();

    auto adapt_method = str2AdaptMethod(ct.gets("adapt_method", "NONE"));
    EstimatedMOProblem<RegELM> estimated_prob(*prob, elms, elmoead.noTrain(),
                                              adapt_method,
                                              elmoead.evalSolutions());
//    WeightSet<double> prob_sel_weights(no_objs, estimators.size()-1, ws_type);
//    const int train_cluster_size = 1.5 * train_weights.size() / sel_weights.size();
//    WeightedEstimatedMOProblem<ELMModel> estimated_prob(*prob, estimators,
//                                                        train_cluster_size,
//                                                        moead_pop,
//                                                        train_weights,
//                                                        moead.weights,
//                                                        prob_sel_weights);



    print(" --- MOEA --- \n");
    cout << *moead << '\n';
    print(" --- ELM --- \n");
    cout << elms[0] << '\n';
    print(" --- SURROGATE MOEA --- \n");
    cout << elmoead << '\n';

    Population moead_pop(moead->popSize(), estimated_prob, NO_SAMPLING);

    auto nd_pop = elmoead.run(*moead, moead_pop, no_moead_eval, estimated_prob, callback);

    /* write result in file */
    auto front_fn = ct.front_fn;
    std::ofstream final_pop_file(front_fn.c_str());
    final_pop_file << nd_pop.objs();
    final_pop_file.close();

    /* print time and hipervolume */
    print("hipervolume: {}\ntime: {}\nfront: {}",
        nd_pop.hypervolume(hyp_ref_vec.data()), Clock::tellTime().count(), front_fn);

    return EXIT_SUCCESS;
}
