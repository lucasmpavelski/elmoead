#pragma once

#include <map>

#include "aux.h"
#include "mo.h"

#include "train_archive.h"
#include "ml.h"

enum class AdaptMethod {
    NONE,
    OFFLINE_INIT,
    ONLINE_MSE
};

AdaptMethod str2AdaptMethod(const std::string& s);

template <typename Regressor>
class EstimatedMOProblem : public MOProblem {
    MOProblem& orig_prob;
    std::vector<Regressor>& regressors;
    const AdaptMethod adapt_method;

    // caches
    int max_no_samples;
    ematrixd norm_all_vars;
    ematrixd norm_all_objs;

public:
    EstimatedMOProblem(MOProblem& orig_prob, std::vector<Regressor>& regressors,
                       const int max_no_samples, AdaptMethod am,
                       const Population& init_pop)
        : MOProblem(orig_prob)
        , orig_prob(orig_prob)
        , regressors(regressors)
        , adapt_method(am)
        , max_no_samples(max_no_samples)
        , norm_all_vars(max_no_samples, this->no_vars)
        , norm_all_objs(this->no_objs, max_no_samples)
    {

        if (adapt_method == AdaptMethod::OFFLINE_INIT)
            adapt(init_pop, init_pop.size());
    }

    virtual ~EstimatedMOProblem() {}

    void evaluate(const double x[], double f[]) final override
    {
        regressors[0].predict(x, f);
    }

    void validate(double x[]) final override { orig_prob.validate(x); }

    void train(const TrainingArchive<double>& train_archive)
    {

        //norm_all_objs = train_archive.solutions.objs().transpose();
        //regressors[0].train(train_archive.solutions.varsDataArr(),
        //                    train_archive.size(),
        //                    norm_all_objs.data());
        auto vars = train_archive.solutions.vars();
        auto objs = train_archive.solutions.objs();

        if (adapt_method == AdaptMethod::ONLINE_MSE)
        {
            using namespace SLFN;
            auto objs2 = ematrixd(objs.rows(), objs.cols());
            std::vector<ELMActFunc> functions = {
                ELMActFunc::SIGMOID,
                ELMActFunc::MULTIQUADRIC,
                ELMActFunc::GAUSSIAN
            };
            ELMActFunc best_act_func = functions[0];
            int best_C = -5;
            double min_error = numeric_limits<double>::infinity();
            for (auto& f : functions) {
                for (int c = -5; c <= 5; c++) {
                    regressors[0].setActFunc(f);
                    regressors[0].setC(pow(2.0, c));
                    regressors[0].train(vars, objs);
                    regressors[0].predictAll(vars, objs2);
                    double error = (objs - objs2).squaredNorm();
                    if (error < min_error) {
                        best_act_func = f;
                        best_C = c;
                        min_error = error;
                    }
                }
            }
            regressors[0].setActFunc(best_act_func);
            regressors[0].setC(pow(2.0, best_C));
            cout << "adapt: " << best_act_func << " " << best_C << '\n';
        }
        regressors[0].train(vars, objs);
    }

    void adapt(const Population& pop, const int pop_size)
    {
        const int no_objs = pop.noObjs();
        const int no_train = pop_size / 2;
        const int no_test = pop_size / 2;

        ematrixd train_vars = pop.vars().topRows(no_train);
        ematrixd train_objs = pop.objs().topRows(no_train);

        ematrixd test_vars = pop.vars().bottomRows(no_test);
        ematrixd test_objs = pop.objs().bottomRows(no_test);

        cout << train_vars << endl;
        cout << train_objs << endl;

        //  elm.adapt(train_vars, train_objs, test_vars, test_objs);

        //        Normalizer<double> norm(-1, 1);
        //        norm.findRangeAndNormalize(pop.varsDataArr(),
        //                                   pop.varsDataArr() + no_vars * pop_size);

        //        AgvNormalizer<double> norm2;
        //        ematrixd t_objs = pop.objs().transpose();
        //        for (int i = 0; i < no_objs; ++i)
        //        {
        //            norm2.findRangeAndNormalize(t_objs.row(i).data(),
        //                                       t_objs.row(i).data() + pop_size);
        //        }
        //        pop.objs() = t_objs.transpose();

        ematrixd pred_objs(no_test, no_objs);

        using namespace SLFN;
        const svector<ELMActFunc> act_funcs = {
            ELMActFunc::SIGMOID,
            ELMActFunc::MULTIQUADRIC,
            ELMActFunc::GAUSSIAN
        };

        const int ce_min = -10;
        const int ce_max = 10;

        const int total_cmins = ce_max - ce_min + 1;
        ematrixd all_scores(total_cmins, act_funcs.size());
        evectord errors(30);

        //        ematrixd mses(total_cmins, act_funcs.size());
        //        ematrixd stds(total_cmins, act_funcs.size());
        const int no_data = no_test * no_objs;

        for (int af = 0; af < act_funcs.size(); ++af) {
            const auto act_fun = act_funcs[af];
            regressors[0].setActFunc(act_fun);

            for (int ce = ce_min; ce <= ce_max; ++ce) {
                const double C = std::pow(2.0, ce);
                regressors[0].setC(C);

                for (int test = 0; test < 30; ++test) {
                    regressors[0].train(train_vars, train_objs);

                    regressors[0].predictAll(test_vars, pred_objs);

                    errors(test) = (test_objs - pred_objs).squaredNorm() / no_data;
                }
                double error_mean = errors.sum() / errors.size();
                double error_std = std::sqrt((errors.array() - error_mean).abs2().sum() / 29.);

                all_scores(ce - ce_min, af) = error_mean + error_std;
                //                mses(ce - ce_min, af) = error_mean;
                //                stds(ce - ce_min, af) = error_std;
            }
        }
        cout << all_scores << endl;
        //        cout << "mses" << mses << endl;
        //        cout << "stds" << stds << endl;

        std::ptrdiff_t ce_idx, af_idx;

        //    all_scores.col(0).minCoeff(&ce_idx);
        //    cout << "sig:" << ce_min + ce_idx << endl;
        //    all_scores.col(1).minCoeff(&ce_idx);
        //    cout << "mul:" << ce_min + ce_idx << endl;
        //    all_scores.col(2).minCoeff(&ce_idx);
        //    cout << "gau:" << ce_min + ce_idx << endl;

        all_scores.minCoeff(&ce_idx, &af_idx);
        cout << "elm_params: " << ce_min + ce_idx << " " << af_idx << endl;

        regressors[0].setC(std::pow(2.0, ce_min + ce_idx));
        regressors[0].setActFunc(act_funcs[af_idx]);


        regressors[0].train(train_vars, train_objs);

        regressors[0].predictAll(test_vars, pred_objs);

        //cout << test_objs << pred_objs << endl;
        //exit(0);
    }
};

template <typename Regressor>
class WeightedEstimatedMOProblem : public MOProblem {
    MOProblem& orig_prob;
    std::vector<Regressor> regressors;

    // caches
    size_t max_no_samples;
    ematrixd train_vars;
    ematrixd train_objs;

    Neighborhood train_associations;
    std::map<const double*, int> inverse_associations;

public:
    WeightedEstimatedMOProblem(MOProblem& orig_prob,
                               std::vector<Regressor>& regressors,
                               const int train_cluster_size,
                               const Population& pop,
                               const WeightSet<double>& train_weights,
                               const WeightSet<double>& eval_weights,
                               const WeightSet<double>& sel_weights)
        : MOProblem(orig_prob)
        , orig_prob(orig_prob)
        , regressors(regressors)
        , max_no_samples(train_cluster_size)
        , train_vars(max_no_samples, this->no_vars)
        , train_objs(this->no_objs, train_cluster_size)
        , train_associations(sel_weights, train_weights, train_cluster_size)
    {
        warning(eval_weights.size() % sel_weights.size() != 0,
                "not all inputs are mapped!");
        Neighborhood eval_associations(sel_weights, eval_weights);
        for (int i = 0; i < eval_associations.no_weights_in; ++i) {
            for (int j = 0; j < eval_associations.no_neighbors; ++j) {
                const int idx = eval_associations[i][j];
                inverse_associations.insert(std::make_pair(pop[idx].vars, i));
            }
        }
    }

    void evaluate(const double x[], double f[]) final override
    {
        const int r_idx = inverse_associations[x];
        regressors[r_idx].predict(x, f);
    }

    void validate(double x[]) final override { orig_prob.validate(x); }

    void train(const TrainingArchive<double>& train_archive)
    {
        for (int i = 0; i < train_associations.no_weights_in; ++i) {
            for (int j = 0; j < train_associations.no_neighbors; ++j) {
                const int idx = train_associations[i][j];
                train_vars.row(j) = train_archive.solutions.vars().row(idx);
                train_objs.col(j) = train_archive.solutions.objs().row(idx);
            }

            regressors[i].train(train_vars.data(),
                                train_associations.no_neighbors,
                                train_objs.data());
        }
    }
};
