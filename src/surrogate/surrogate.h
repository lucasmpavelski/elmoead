#pragma once

#include "aux.h"
#include "mo.h"
#include "moead_algs.h"
#include "train_archive.h"
#include "weighted_selector.h"

class ELMOEAD {
private:
    MOProblem& orig_prob;
    const size_t no_max_eval;
    const size_t no_init_sol;
    NdArchive nd_archive;
    TrainingArchive<double> train_archive;
    WeightedSelector<double> selector;
    //Plotter plt;

public:
    ELMOEAD(MOProblem& orig_prob,
            const int no_max_eval,
            const Population& init_sol,
            const AggregationFunctionType aggr_func,
            const WeightSet<>& moead_weights,
            const WeightSet<>& train_weights,
            const WeightSet<>& sel_weights,
            const double eps)
        : no_max_eval(no_max_eval)
        , orig_prob(orig_prob)
        , no_init_sol(train_weights.size())
        , nd_archive(orig_prob, no_max_eval, eps)
        , train_archive(train_weights, aggr_func, init_sol)
        , selector(moead_weights, sel_weights)
    {
        nd_archive.add(train_archive.solutions);
    }

    size_t noTrain() const { return train_archive.size(); }
    size_t noSelect() const { return selector.noSets(); }
    const Population& evalSolutions() const { return nd_archive.all(); }
    size_t noEvalSolutions() const { return nd_archive.size(); }

    template <typename EstimatedMOProblem, typename Callback>
    Population run(MOEAD& moead, Population& moead_pop,
                   const int no_moead_evals,
                   EstimatedMOProblem& estimated_prob,
                   Callback cb)
    {
        Population nd_pop(estimated_prob);

        nd_archive.setProblem(&estimated_prob);
        moead_pop.setProblem(&estimated_prob);

        // init selector
        const int no_sets = selector.noSets();
        const int no_ranked = selector.noRanked();
        std::vector<int> ranked_solutions(no_ranked);

        size_t no_eval = no_init_sol;
        cb(no_eval, nd_archive.all());
        while (no_eval < no_max_eval) {
            estimated_prob.train(train_archive);

            // FIND AN APPROXIMATION
            nd_archive.getNdPopulation(nd_pop);
            //nd_pop.plotObjs(plt);
            reinitMOEADPopulation(nd_pop, moead_pop);

            moead.optimize(no_moead_evals, estimated_prob, moead_pop, NullCallback);
            selector.rankSolutions(moead.weights, moead_pop, moead.aggr_prob,
                                   ranked_solutions);

            Solution const* cands[no_ranked];
            bool reselect = false;
            int no_selected = 0;
            do {
                no_selected = 0;
                reselect = false;
                for (int i = 0; i < no_ranked; ++i) {
                    const int idx = ranked_solutions[i];
                    Solution& candidate = moead_pop[idx];

                    if (!nd_archive.contain(candidate)) {
                        orig_prob.evaluate(candidate.vars, candidate.objs);
                        no_eval++;
                        cands[no_selected] = &candidate;

                        //if (no_eval == no_init_sol + int(no_evals_remain / 3.))
                        //    estimated_prob.adapt(train_archive.solutions, train_archive.size());
                        //if (no_eval == no_init_sol + int(2 * no_evals_remain / 3.))
                        //    estimated_prob.adapt(train_archive.solutions, train_archive.size());
                        nd_archive.add(candidate);
                        //train_archive.add(candidate);
                        cb(no_eval, nd_archive.all());

                        no_selected++;

                        if ((no_selected >= no_sets) || (no_eval >= no_max_eval))
                            break;
                    }
                }
                train_archive.add(cands, no_selected);

                if (no_selected == 0) {
                    nd_archive.setEps(nd_archive.getEps() / 2.0);
                    reselect = true;
                }
            } while (reselect);
        };

        // return non-dominated solutions (with original problem reference)
        nd_archive.getNdPopulation(nd_pop);
        nd_pop.setProblem(&orig_prob);

        return nd_pop;
    }

    void reinitMOEADPopulation(Population& given, Population& moead_pop)
    {
        const int pop_size = moead_pop.size();
        const int no_given = given.size();
        //(given.size() < 0.9 * pop_size) ? given.size() : 0.9 * pop_size;

        given.shuffle();
        for (int i = 0; (i < no_given) && (i < pop_size); ++i)
            moead_pop[i] = given[i];

        if (no_given < pop_size) {
            double inv_n = 1.0 / no_given;
            evectord mean = given.vars().colwise().sum() * inv_n;
            evectord std = ((given.vars().rowwise() - mean.transpose())
                                .colwise()
                                .squaredNorm() * inv_n).cwiseSqrt();
            evectord lb = mean - std;
            evectord ub = mean + std;
            moead_pop.sampleVars(LATIN_HYP_SAMPLING, moead_pop.begin() + no_given,
                                 moead_pop.end(), lb.data(), ub.data());
            for (int i = no_given; i < pop_size; ++i)
                moead_pop[i].validateEvaluate();
        }
    }

    friend std::ostream& operator<<(std::ostream& os, const ELMOEAD& em)
    {
        return os << "problem: " << em.orig_prob.getName() << '\n'
                  << "no_max_eval: " << em.no_max_eval << '\n'
                  << "archive_eps: " << em.nd_archive.getEps() << '\n';
    }
};
