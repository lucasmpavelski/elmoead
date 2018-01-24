#pragma once

#include "../mo.h"
#include "../aux.h"
#include "../obj_dec.h"

class STMData {
private:
    // preference ranks
    ematrixi psi_p;
    //int** psi_x;

    // aux
    evectori propose_counts;
    evectori f_x;
    evectori free_probs;

    // return
    evectori selected;

public:
    STMData(const int pop_size, const int no_problems) :
        psi_p(no_problems, pop_size),
        propose_counts(no_problems),
        f_x(pop_size),
        free_probs(no_problems),
        selected(pop_size)
    {};

    virtual ~STMData() {}

    Solution* sort(Population &pop, const WeightSet<double>& w,
                   const AggrProblem<double>& prob);

};

void stm(const Population& pop, const WeightSet<double>& w,
         const AggrProblem<double>& prob, const ematrixi &psi_p,
         evectori &propose_counts, evectori &free_probs, evectori& f_x,
         evectori &selected);

int* stm(const int pop_size, const int no_problems,
         const ematrixi& psi_p, const int* const* psi_x,
         int propose_counts[], int free_probs[], bool f_x[],
         int selected[]);

template<typename T>
void solutionPreferences(const Population& pop, const WeightSet<T>& weights,
                         const AggrProblem<T>& aggr_prob, ematrixi& psi_x)
{
    const int no_sols = pop.size();
    const int no_prob = weights.size();
    psi_x.resize(no_sols, no_prob);
    svector<Indexed<T>> delta_x(no_prob);
    for (int i = 0; i < no_sols; ++i)
    {
        const auto& sol = pop[i];
        for (int j = 0; j < no_prob; ++j)
            delta_x[j] = Indexed<T>(j, aggr_prob.affinity(sol.objs, weights[j]));
        sort(delta_x.begin(), delta_x.end(), Indexed<T>::compareByVal);
        for (int j = 0; j < no_prob; ++j)
            psi_x(i, j) = delta_x[j].idx;
    }
}

template<typename T>
void problemPreferences(const Population& pop, const WeightSet<T>& weights,
                        const AggrProblem<T>& aggr_prob, ematrixi& psi_p)
{
    const int no_sols = pop.size();
    const int no_prob = weights.size();
    psi_p.resize(no_prob, no_sols);
    svector<Indexed<T>> delta_p(no_sols);
    for (int i = 0; i < no_prob; ++i)
    {
        const auto& w = weights[i];
        for (int j = 0; j < no_sols; ++j)
            delta_p[j] = Indexed<T>(j, aggr_prob.eval(pop[j].objs, w));
        sort(delta_p.begin(), delta_p.end(), Indexed<T>::compareByVal);
        for (int j = 0; j < no_sols; ++j)
            psi_p(i, j) = delta_p[j].idx;
    }
}
