#include "stm.h"

#include <algorithm>
using std::sort;
using std::iota;
using std::fill;
using std::find_if;


Solution *STMData::sort(Population& pop, const WeightSet<double>& w,
                        const AggrProblem<double>& prob)
{

    problemPreferences(pop, w, prob, psi_p);

    //computeSolsPreferences(sols, pop_size, w, no_probs, ideal, nadir,
    //                       delta_aux, psi_x);

    //stm(pop_size, no_probs, psi_p, psi_x, propose_counts, free_probs, f_x,
    //    prob,  w, sols, selected);

    stm(pop, w, prob, psi_p, propose_counts, free_probs, f_x, selected);

    pop.sort(selected.data());

    return pop.solsData();
}

void stm(const Population& pop, const WeightSet<double>& w,
         const AggrProblem<double>& prob, const ematrixi &psi_p,
         evectori& propose_counts, evectori& free_probs, evectori& f_x,
         evectori& selected)
{
    const int pop_size = pop.size();
    const int no_problems = w.size();

    std::iota(free_probs.data(), free_probs.data() + no_problems, 0);
    f_x.fill(0);
    propose_counts.fill(0);
    selected.fill(no_problems + 1);


    int no_free_problems = no_problems;
    while (no_free_problems > 0)
    {
        // randomly chose a free subproblem
        const int free_problem_idx = RNG::intUniform(no_free_problems - 1);
        const int free_problem = free_probs[free_problem_idx];

        // find most prefered solution of free_problem that it didnt tried
        // to match yet
        const int propose_next = propose_counts[free_problem];
        const int most_prefered = psi_p(free_problem, propose_next);
        propose_counts[free_problem]++;

        // propose to match
        // if it doesnt have a partner
        if (!f_x[most_prefered])
        {
            // match <3
            selected[most_prefered] = free_problem;
            f_x[most_prefered] = true;
            free_probs[free_problem_idx] = free_probs[no_free_problems-1];
            no_free_problems--;
        }
        else
        {
            // find most_prefered current partner
            const int current_partner = selected[most_prefered];

            // see if she prefers him to her actual partner
            double const* sol_objs = pop[most_prefered].objs;
            const double curr_pref = prob.affinity(sol_objs,
                                                            w[current_partner]);
            const double prob_pref = prob.affinity(sol_objs,
                                                            w[free_problem]);

            if (prob_pref < curr_pref)
            {
                // forget the other
                free_probs[free_problem_idx] = current_partner;
                // match <3
                selected[most_prefered] = free_problem;
            }
        }
    }
}

int* stm(const int pop_size, const int no_problems,
         const ematrixi& psi_p, const int* const* psi_x,
         int propose_counts[], int free_probs[], bool dum[],
         int selected[])
{
    svectorb f_x(pop_size);

    if (selected == NULL)
        selected = new int[no_problems];

    std::iota(free_probs, free_probs + no_problems, 0);
    using std::fill_n;
    fill_n(propose_counts, no_problems, 0              );
    fill_n(selected      , pop_size   , no_problems + 1);
    fill_n(f_x.begin()   , pop_size   , false          );

    int no_free_problems = no_problems;
    while (no_free_problems > 0)
    {
        // randomly chose a free subproblem
        const int free_problem_idx = RNG::intUniform(no_free_problems - 1);
        const int free_problem = free_probs[free_problem_idx];

        // find most prefered solution of free_problem that it didnt tried
        // to match yet
        const int propose_next = propose_counts[free_problem];
        const int most_prefered = psi_p(free_problem, propose_next);
        propose_counts[free_problem]++;

        // propose to match if it doesnt have a partner
        if (!f_x[most_prefered])
        {
            // match <3
            selected[most_prefered] = free_problem;
            f_x[most_prefered] = true;
            free_probs[free_problem_idx] = free_probs[no_free_problems-1];
            no_free_problems--;
        }
        else
        {
            // find most_prefered current partner
            const int current_partner = selected[most_prefered];

            // see if she prefers him to her actual partner
            const int* her_preferences = psi_x[most_prefered];
            const int* better_pretender =
                    find_if(her_preferences, her_preferences + no_problems,
                            [free_problem,current_partner](int p) {
                return (p == free_problem) || (p == current_partner);
            });

            if (*better_pretender == free_problem)
            {
                // forget the other
                free_probs[free_problem_idx] = current_partner;
                // match <3
                selected[most_prefered] = free_problem;
            }
        }
    }

    return selected;
}
