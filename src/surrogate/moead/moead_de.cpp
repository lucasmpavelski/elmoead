#include "moead_de.h"


void MOEAD_DE::optimize(const int no_evals, MOProblem& prob, Population& pop,
                        PopulationCallback callback)
{
    AllocSolution tmp(prob);

    //Step 1.4: Initialize the reference points
    initReferencePoints(pop);

    //sortPopulationByWeights(pop);
    //stm_data.sort(pop, weights, aggr_prob);

    for (int i = 0; i < pop_size; ++i)
        fitness_aux[i] = subProbEval(pop, i);
    //operators.reset(std::min_element(fitness_aux.begin(), fitness_aux.end());

    //Step 2: Reproduction and Update
    int eval_counter = 0, gen = 0;
    while (eval_counter < no_evals)
    {
        for (int i = 0; i < pop_size; ++i)
        {
            //Step 2.1: Select mating parents from the neighborhood of x(i)
            const bool use_neighborhood = (RNG::realUniform<double>() < real_b);
            const int* candidates = use_neighborhood ? neighbors[i] : P.data();
            const int no_candidates = use_neighborhood ? neighbors.no_neighbors : pop_size;

            evolution(pop, i, tmp, candidates, no_candidates);
            eval_counter++;

            // Step 2.4: Update the population.
            update(pop, tmp, candidates, no_candidates);
        }
        operators.update();
        callback(eval_counter, pop);
    }
}

void MOEAD_DE::evolution(const Population& pop, const int idx,
                            Solution& result, const int candidates_idxs[],
                            const int no_candidates)
{
    const int no_vars = result.noVars();
    const MOProblem::dvector& lb = result.problem->lower_bounds;
    const MOProblem::dvector& ub = result.problem->upper_bounds;

    const DEOperator& de_operator = operators.selectOperator();

    const int no_parents = DEOperator::noParentsFor(de_operator.type());
    const double* parents[no_parents];

    DEOperator::selectParents(de_operator.type(), pop.solsData(), candidates_idxs,
                              no_candidates, idx, parents);

    de_operator(parents, no_vars, result.vars);
    result.validate();
    pm_operator.mutate(lb.data(), ub.data(), no_vars, result.vars);
    result.validateEvaluate();
    // Step 2.3: Update the reference point f
    if (updateReferencePoints(result))
    {
        for (int i = 0; i < pop_size; ++i)
            fitness_aux[i] = subProbEval(pop, i);
    }
};

int MOEAD_DE::update(Population& pop, const Solution& ind,
                        const int candidates_idxs[], const int no_candidates)
{
    std::vector<int> perm(candidates_idxs, candidates_idxs + no_candidates);
    std::shuffle(perm.begin(), perm.end(), RNG::engine);
    int count = 0;
    for (const auto& k : perm)
    {
        const auto f_new = subProbEval(ind, k);
        if (f_new < fitness_aux[k])
        {
            operators.feedback(f_new, fitness_aux[k]);
            fitness_aux[k] = f_new;
            pop[k] = ind;
            count++;
            if (count >= update_limit)
                break;
        }
    }
    return count;
}
