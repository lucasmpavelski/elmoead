#include "moead_dra.h"

void MOEAD_DRA::rememberObjs(const Population &pop)
{
    for (int i = 0; i < pop_size; ++i)
    {
        for (int j = 0; j < no_objs; ++j)
            previous_objs(i, j) = pop[i].objs[j];
    }
}

void MOEAD_DRA::tourSelection(const int depth)
{
    // selection based on utility
    // select extremes
    assert(no_objs == 2);
    selected[0] = 0;
    selected[1] = pop_size - 1;
    int no_selected = no_objs;

    for (int i = 1; i < pop_size - 1; ++i)
        candidates[i - 1] = i;
    int no_candidates = pop_size - 2;

    while (no_selected < pop_size / 5)
    {
        int best_idx = rand() % no_candidates;
        int best_i = candidates[best_idx];
        for (int i = 1; i < depth; ++i)
        {
            const int cand_idx = rand() % no_candidates;
            const int cand_i = candidates[best_idx];
            if (utility[cand_i] > utility[best_i])
            {
                best_idx = cand_idx;
                best_i = cand_i;
            }
        }
        no_candidates--;
        std::swap(candidates[best_idx], candidates[no_candidates]);
        selected[no_selected] = best_i;
        no_selected++;
    }
}

void MOEAD_DRA::compUtility(const Population& pop)
{
    for (int i = 0; i < pop_size; ++i)
    {
        const double f_new = subProbEval(pop, i);
        const double f_old = subProbEval(previous_objs.row(i).data(), i);
        const double delta = (f_old - f_new) / f_old;
        if (delta > 0.001)
        {
            utility[i] = 1.0;
        }
        else
        {
            // uti = 0.95*(1.0+delta/0.001)*utility_[n];
            const double uti = (0.95 + delta * (0.05 / 0.001)) * utility[i];
            utility[i] = uti < 1.0 ? uti : 1.0;
        }
    }
}

void MOEAD_DRA::optimize(const int no_evals, MOProblem& prob, Population& pop,
                         PopulationCallback callback)
{
    AllocSolution tmp(prob);

    rememberObjs(pop);

    //Step 1.4: Initialize the reference points
    initReferencePoints(pop);

    //Step 2: Reproduction and Update
    int eval_counter = 0, gen = 0;
    while(eval_counter < no_evals)
    {
        tourSelection(10);
        for (int i = 0; i < pop_size / 5; ++i)
        {
            //Step 2.1: Select mating parents from the neighborhood of x(i)
            int s = selected[i];
            const bool use_neighborhood = (RNG::realUniform<double>() < realB());
            const int* candidates = use_neighborhood ? neighbors[s] : P.data();
            const int no_candidates = use_neighborhood ?
                        neighbors.no_neighbors : pop_size;

            evolution(pop, s, tmp, candidates, no_candidates);
            eval_counter++;

            update(pop, tmp, candidates, no_candidates);
        }

        //operators.updateProbabilities();

        callback(eval_counter, pop);
        gen++;
        if (gen % 50 == 0)
        {
            compUtility(pop);
            rememberObjs(pop);
        }
    }
}
