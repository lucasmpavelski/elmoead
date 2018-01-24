#include "moead.h"

void MOEAD::evolution(const Population& pop, const int idx, Solution& result,
                      const int* candidates_idxs, const int no_candidates)
{
    const auto no_vars = result.noVars();
    const auto& lb = result.problem->lower_bounds;
    const auto& ub = result.problem->upper_bounds;

    double* out = result.vars;
    const int tournament_size = 2;

    auto t_eval = [this,&pop](const int& a, const int& b) {
        return this->subProbEval(pop, a) < this->subProbEval(pop, b);
    };

    int b1 = tournament(candidates_idxs, no_candidates, tournament_size, t_eval);
    int b2 = tournament(candidates_idxs, no_candidates, tournament_size, t_eval);

    double dvars[no_vars];

    sbx_operator(lb.data(), ub.data(), no_vars,
                 pop[b1].vars, pop[b2].vars,
                 out, dvars);

    result.validate();
    pm_operator.mutate(lb.data(), ub.data(), no_vars, out);
    result.evaluate();
}

int MOEAD::update(Population& pop, const Solution& ind,
                  const int* candidates_idxs, const int no_candidates)
{
    for (int i = 0; i < no_candidates; ++i)
    {
        const int k = candidates_idxs[i];
        const double f_k = fitness_aux[k];
        const double f_new = subProbEval(ind, k);

        if (f_new < f_k)
        {
            pop[k] = ind;
            fitness_aux[k] = f_new;
        }
    }
    return no_candidates;
}

void MOEAD::sortPopulationByWeights(Population& pop)
{
    const int pop_size = pop.size();
    for (int i = 0; i < pop_size; ++i)
    {
        double best_f = subProbEval(pop, i);
        int best_f_idx = i;

        for (int j = i + 1; j < pop_size; ++j)
        {
            const double f = subProbEval(pop[j], i);

            if (f < best_f)
            {
                best_f = f;
                best_f_idx = j;
            }
        }

        if (i != best_f_idx)
            pop[best_f_idx].swap(pop[i]);
            //std::swap(, );
    }
}

void MOEAD::initReferencePoints(const Population& pop)
{
    aggr_prob.resetReferencePoints();
    for (int i = 0; i < pop.size(); i++)
        aggr_prob.updateReferencePoints(pop[i].objs);
    for (int i = 0; i < pop_size; i++)
        fitness_aux[i] = subProbEval(pop, i);
}

bool MOEAD::updateReferencePoints(const Solution& sol)
{
    return aggr_prob.updateReferencePoints(sol.objs);
}

void MOEAD::optimize(const int no_evals, MOProblem& prob, Population& pop,
                     PopulationCallback callback)
{
    AllocSolution tmp(prob);

    //Step 1.4: Initialize the reference points
    initReferencePoints(pop);
    //sortPopulationByWeights(pop, initPopSize());
    //stmSort(pop, weights, ideal, nadir, type, stm_data);

    for (int i = 0; i < pop_size; ++i)
        fitness_aux[i] = subProbEval(pop, i);

    //Step 2: Reproduction and Update
    int eval_counter = 0;
    while(eval_counter < no_evals)
    {
        for (int i = 0; i < pop_size; ++i)
        {
            //Step 2.1: Select mating parents from the neighborhood of x(i) ;
            const int* candidates = neighbors[i];
            const int no_candidates = neighbors.no_neighbors;
            evolution(pop, i, tmp, candidates, no_candidates);
            eval_counter++;
            // Step 2.3: Update the reference point f
            if (updateReferencePoints(tmp))
            {
                for (int j = 0; j < pop_size; ++j)
                    fitness_aux[j] = subProbEval(pop, j);
            }
            // Step 2.4: Update the population.
            update(pop, tmp, candidates, no_candidates);
        }
        callback(eval_counter, pop);
    }
}
