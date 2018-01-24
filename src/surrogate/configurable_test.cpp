#include "configurable_test.h"
#include "benchmark.h"
#include "adp.h"

using std::string;

ConfigurableTest::ConfigurableTest(int argc, char* argv[])
    : config_fn((argc > 1) ? argv[1] : "config.ini")
    , front_fn((argc > 2) ? argv[2] : "final_pop.dat")
    , conf(config_fn.c_str())
{
    seed = getl("seed", RNG::trueRandom());
    RNG::seed(seed);
    Clock::begin();
}

MOProblem* ConfigurableTest::getProblem()
{
    auto name = gets("prob_name", "WFG1");
    if (name == "ADP")
        return new AirfoilDesignProblem();
    return getBenchmark(name.c_str(),
        getl("prob_no_vars", 24),
        getl("prob_no_objs", 2));
}

OperatorSelection<DEOperator> ConfigurableTest::getDEOperators()
{
    const string de_operator = gets("de_operator", "DE/RAND/1/BIN");
    const double de_cr = getd("de_cr", 1.0);
    const double de_f = getd("de_f", 0.5);

    /* de operator types */
    std::vector<DEOperator> de_ops;
    std::stringstream ss(de_operator);
    std::transform(std::istream_iterator<string>(ss),
        std::istream_iterator<string>(),
        std::back_inserter(de_ops),
        [de_cr, de_f](const string tp) {
            return DEOperator(str2DEType(tp), de_cr, de_f);
        });
    /* de operators adaptation */
    OperatorSelection<DEOperator> operators(de_ops); // const operator
    /* TODO: make adaptation configurable */
    //    using PMRew = ProbabilityMatching<DEOperator>::RewardType;
    //    ProbabilityMatching<DEOperator> operators(de_ops,
    //                                              //PMRew::AvgAbs,
    //                                              //PMRew::AvgNorm,
    //                                              //PMRew::ExtAbs,
    //                                              PMRew::ExtNorm
    //                                              );
    //    FRRMAB<DEOperator> operators(de_ops, weights.size() / 2);
    return operators;
}

MOEAD* ConfigurableTest::getMOEAD(size_t no_vars, size_t no_objs)
{
    const string weight_type = gets("weight_type", "UNIFORM");
    const string moead_type = gets("moead_type", "MOEA/D-DE");
    const string moead_aggr_func = gets("moead_aggr_func", "PBI");
    const long moead_no_partitions = getl("moead_no_partitions", 99);
    const long moead_no_neighbors = getl("moead_no_neighbors", 10);
    const long moead_no_updates = getl("moead_no_updates", 1);
    const double moead_real_b = getd("moead_real_b", 0.9);
    const double moead_sbx_xover = getd("moead_sbx_xover", 1.0);
    const double moead_sbx_eta = getd("moead_sbx_eta", 20.0);
    const double moead_pm_mut = getd("moead_pm_mut", 1.0 / no_vars);
    const double moead_pm_eta = getd("moead_pm_eta", 20.0);

    const auto aggr_func = str2AggrFunc(moead_aggr_func);

    auto ws_type = (weight_type == "UNIFORM") ? WeightSet<>::Type::UNIFORM : WeightSet<>::Type::RANDOM;
    const WeightSet<double> weights(no_objs, moead_no_partitions, ws_type);

    const PMOperator pmo(moead_pm_mut, moead_pm_eta);
    const SBXoverOperator sbxo(moead_sbx_xover, moead_sbx_eta);

    auto operators = getDEOperators();

    return getMoeadAlg(moead_type, weights, moead_no_neighbors,
        aggr_func, moead_real_b, moead_no_updates,
        sbxo, operators, pmo);
}

std::vector<RegELM> ConfigurableTest::getELMs(const MOProblem& prob)
{
    const string elm_act_func = gets("elm_act_func", "MULTIQUADRIC");
    const double elm_C = getd("elm_C", 1.0);
    const long elm_no_hidden = getl("elm_no_hidden", 2 * prob.noVars() + 1);
    const bool elm_norm_input = getb("elm_norm_input", false);
    const bool elm_norm_output = getb("elm_norm_output", true);

    return { RegELM(prob.noVars(), elm_no_hidden, prob.noObjs(),
        SLFN::str2ELMActFunc(elm_act_func), elm_C,
        elm_norm_input, elm_norm_output,
        prob.lower_bounds, prob.upper_bounds) };

    // return {SLFN::RBF(prob.noVars(), elm_no_hidden, prob.noObjs(), prob.lower_bounds, prob.upper_bounds) };
}

size_t ConfigurableTest::noTrainPartitions(size_t no_vars, size_t no_objs,
    WeightSet<>::Type wt)
{
    if (wt == WeightSet<>::Type::RANDOM)
        return 10 * no_vars;
    if (no_objs == 2)
        return 10 * no_vars - 1;
    size_t no_part = 1;
    while (WeightSet<>::sizeFor(no_objs, no_part) < 10 * no_vars)
        no_part++;
    return no_part;
}

ELMOEAD ConfigurableTest::getELMOEAD(MOProblem& prob, MOEAD& moead)
{
    auto no_evals = getl("no_evals", 1000);
    auto no_sel_partitions = getl("no_sel_partitions", 9);
    auto archive_eps = getd("archive_eps", 0.001);

    auto no_vars = prob.noVars();
    auto no_objs = prob.noObjs();
    auto wtype = moead.weights.type();

    WeightSet<> sel_weights(no_objs, no_sel_partitions, wtype);

    auto no_train_part = noTrainPartitions(no_vars, no_objs, wtype);
    WeightSet<> train_weights(no_objs, no_train_part, moead.weights.type());

    Population init_sol(train_weights.size(), prob, LATIN_HYP_SAMPLING);
    init_sol.evaluateAll();

    return ELMOEAD(prob, no_evals, init_sol, moead.aggrFunc(),
        moead.weights, train_weights, sel_weights, archive_eps);
}

PopulationCallback ConfigurableTest::getCallback()
{
    const string cb = gets("callback_type", "NONE");
    const size_t freq = getl("callback_frequency", 1);
    if (cb == "NONE")
        return NullCallback;
    if (cb == "PLOT")
        return PlotCallback(freq);
    if (cb == "LOG")
        return LogCallback(front_fn, freq);
    throw_assert(false, "callback_type must be (NONE | PLOT | LOG)");
    return NullCallback;
}
