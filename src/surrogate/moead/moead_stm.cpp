#include "moead_stm.h"

void MOEAD_STM::optimize(const int no_evals, MOProblem& prob,
                            Population& pop, PopulationCallback callback)
{
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
            int s = selected[i];
            //Step 2.1: Select mating parents from the neighborhood of x(i)
            const bool use_neighborhood = (RNG::realUniform<double>() < realB());
            const int* candidates = use_neighborhood ? neighbors[s] : P.data();
            const int no_candidates = use_neighborhood ? neighbors.no_neighbors : popSize();

            evolution(pop, s, pop[pop_size + i], candidates, no_candidates);
            eval_counter++;
        }
        stm_data.sort(pop, weights, aggr_prob);

        callback(eval_counter, pop);
        gen++;
        if (gen % 50 == 0)
        {
            compUtility(pop);
            rememberObjs(pop);
        }
    }
}
