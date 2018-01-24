#pragma once

#include <cstring>
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>

#include <Eigen/Core>

#include "aux.h"
#include "obj_dec.h"
#include "mo.h"
#include "moead/stm.h"

template <typename T>
class TrainingArchive {
public:
    using size_t = std::size_t;

private:
    //STMData stm;
public:
    typedef typename WeightSet<T>::vector vector;

    TrainingArchive(const WeightSet<T>& w,
                    const AggregationFunctionType aggr_func,
                    const Population& init_pop)
        : weights(w)
        , aggr_problem(aggr_func, w.dimension())
        , solutions(init_pop)
        , fitness_aux(init_pop.size())
        //, stm(init_pop.size() + 10, init_pop.size())
    {
        int i = 0;
        for (auto& s : solutions) {
            s.setVars(init_pop[i].vars);
            s.setObjs(init_pop[i].objs);
            aggr_problem.updateReferencePoints(s.objs);
            i++;
        }
        //STMData stm0(init_pop.size(), weights.size());
        //stm0.sort(solutions, weights, aggr_problem);
        updateFitness();
    }

    size_t size() const { return weights.size(); }
    size_t noVars() const { return solutions.noVars(); }
    size_t noObjs() const { return solutions.noObjs(); }

    void updateFitness()
    {
        for (int i = 0; i < size(); ++i) {
            T best_f = aggr_problem.eval(solutions[i].objs, weights[i]);
            int best_f_idx = i;
            for (int j = i + 1; j < size(); ++j) {
                const T new_f = aggr_problem.eval(solutions[j].objs, weights[i]);
                if (new_f < best_f) {
                    best_f = new_f;
                    best_f_idx = j;
                }
            }

            if (best_f_idx != i)
                solutions[i].swap(solutions[best_f_idx]);

            fitness_aux[i] = best_f;
        }
    }

//    template <class Itr>
//    void add(Itr b, Itr e) {
//        size_t no_new = std::distance(b, e);
//        curr = 0;
//        solutions.resize(size() + no_new);
//        for (auto i = b; i != e; i++)
//        {
//            aggr_problem.updateReferencePoints(i->objs);
//            solutions[size() + curr++] = *i;
//        }
//        STMData stm0(size() + no_new, size());
//        stm0.sort(solutions, weights, aggr_problem);
//    }

    void add(const Solution& new_sol)
    {
        Solution const* p[1] = { &new_sol };
        add(p, 1);
    }

    void add(Solution const* new_sols[], const int no_new_sols)
    {
        aggr_problem.resetReferencePoints();
        for (const auto& s : solutions)
            aggr_problem.updateReferencePoints(s.objs);

        bool ref_change = false;
        for (int i = 0; i < no_new_sols; ++i) {
            if (aggr_problem.updateReferencePoints(new_sols[i]->objs))
                ref_change = true;
        }

        if (ref_change) {
            updateFitness();
            //for (int i = 0; i < size; ++i)
            //    fitness_aux[i] = aggr_problem.eval(solutions[i].objs, weights[i]);
        }

        std::vector<int> assigned(no_new_sols, -1);

        for (int i = 0; i < no_new_sols; ++i) {
            const Solution& sol = *new_sols[i];
            /*
            double max_dec = -std::numeric_limits<double>::infinity();
            int max_dec_idx = -1;
            for (int j = 0; j < size(); ++j)
            {
                if (contain(assigned.begin(), assigned.begin() + i, j))
                    continue;
                const T new_f = aggr_problem.eval(sol.objs, weights[j]);
                const T old_f = fitness_aux[j];
                const T dec = (old_f - new_f);
                if (dec > max_dec)
                {
                    max_dec = dec;
                    max_dec_idx = j;
                }
            }

            fitness_aux[max_dec_idx] = aggr_problem.eval(sol.objs, weights[max_dec_idx]);
            solutions[max_dec_idx] = sol;

            assigned[i] = max_dec_idx;
/*/
            T best_f = std::numeric_limits<double>::infinity(); //aggr_problem.eval(sol.objs, weights[0]);
            int best_f_idx = -1;
            for (int j = 0; j < size(); ++j)
            {
              //  if (assigned[j] == -1)
                {
                    const T new_f = aggr_problem.eval(sol.objs, weights[j]);
                    if (new_f < best_f) {
                        best_f = new_f;
                        best_f_idx = j;
                    }
                }
            }

            int selected_idx = -1;
            if (best_f < fitness_aux[best_f_idx]) {
                // if the solution is better than current solution associated with the
                // weight, replace it
                selected_idx = best_f_idx;
            } else {
                // else, find the closest solution and replace it
                selected_idx = closestSolutionIdx(sol);
                best_f = aggr_problem.eval(sol.objs, weights[selected_idx]);
            }
            assigned[i] = selected_idx;
            fitness_aux[selected_idx] = best_f;
            solutions[selected_idx].setVars(sol.vars);
            solutions[selected_idx].setObjs(sol.objs);
            //*/
        }
    }

    /*   const int size;
    const int no_vars;
    const int no_objs;*/
    WeightSet<T> weights;
    AggrProblem<T> aggr_problem;
    Population solutions;
    vector fitness_aux;

private:
    int closestSolutionIdx(const Solution& s)
    {
        const int no_objs = noObjs();
        auto r = std::min_element(solutions.begin(), solutions.end(),
                                  [&s, no_objs](const Solution& a, const Solution& b) {
            return sqrEuclideanDist(a.objs, s.objs, no_objs) <
                   sqrEuclideanDist(b.objs, s.objs, no_objs);
        });
        return std::distance(solutions.begin(), r);
    }
};

/*
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double** inputs;
    double** outputs;
    int no_ins;
    int no_outs;

    int no_weights;
    double** weights;
    double* scores;
    double* ideal;
    double* nadir;
    AggregationFunctionType type;
} TrainArchive;

void initTrainArchive(TrainArchive* ta, const int no_ins, const int no_outs,
                      const int no_partitions, const AggregationFunctionType type,
                      double const* const init_ins[],
                      double const* const init_outs[]);
void deinitTrainArchive(TrainArchive* ta);

void getTrainArchiveInputs(const TrainArchive* ta, double** inputs);
void getTrainArchiveOutputs(const TrainArchive* ta, const int idx,
                            double* outputs);

void addToTrainArchive(TrainArchive* ta, const double inputs[],
                       const double outputs[]);

#ifdef __cplusplus
}
#endif

*/
