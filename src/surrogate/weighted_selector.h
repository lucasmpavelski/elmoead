#pragma once

#include <assert.h>
#include <math.h>

#include "aux.h"
#include "mo.h"
#include "obj_dec.h"

template< typename T >
class WeightedSelector
{
public:
    WeightedSelector(const WeightSet<T>& sol_weights,
                     const WeightSet<T>& sel_weights) :
        sel_weights(sel_weights),
        associations(sel_weights, sol_weights,
                     sol_weights.size() / sel_weights.size()),
        fitness_aux(sol_weights.size() / sel_weights.size()),
        ranks_aux(sol_weights.size() / sel_weights.size())
    {};

    int setSize()  const { return associations.no_neighbors; };
    int noSets()   const { return sel_weights.size();        };
    int noRanked() const { return setSize() * noSets();      };

    void rankSolutions(const WeightSet<T>& sol_weights, const Population& sols,
                       const AggrProblem<T>& aggr_prob,
                       std::vector<int>& ranked_idxs) {
        const int no_sel_weights = sel_weights.size();
        const int w_per_ws = associations.no_neighbors;

        for (int i = 0; i < no_sel_weights; ++i)
        {
            /*printf("\nsel_weight: ");
            printDblVector(stdout, ws->sel_weights[i], no_objs);
            printf("\nsols: \n");*/

            const int* assocs = associations[i];
            for (int j = 0; j < w_per_ws; ++j)
            {
                const int idx = assocs[j];

                /*printf("  %d (f:[", idx);
                printDblVector(stdout, sols[idx].objs, no_objs);
                printf("],w:[");
                printDblVector(stdout, w, no_objs);
                printf("],ideal[");
                printDblVector(stdout, ideal, no_objs);
                printf("])\n");
                fflush(stdout);*/

                fitness_aux[j] = aggr_prob.eval(sols[idx].objs, sol_weights, idx);
                //indexes[i][j] = idx;
            }

            /*printf("fitnesses: ");
            printDblVector(stdout, fitness, w_per_ws);*/

            //sortIdx(fitness, indexes[i], w_per_ws, w_per_ws);
            std::iota(ranks_aux.begin(), ranks_aux.end(), 0);
            std::sort(ranks_aux.begin(), ranks_aux.end(),
                      [this](const int& a, const int& b) {
                return this->fitness_aux[a] < this->fitness_aux[b];
            });

            for (int j = 0; j < w_per_ws; ++j)
                ranked_idxs[j * no_sel_weights + i] = assocs[ranks_aux[j]];

            /*printf("ranked: ");
            printIntVector(stdout, indexes[i], w_per_ws);*/
        }
    };

    friend std::ostream& operator<<(std::ostream& os,
                                    const WeightedSelector& ws) {
        os << "weighted_selector:\n"
           << "  sel_weights: \n" << ws.sel_weights << "\n"
           << "  associations: \n" << ws.associations << "\n";
        return os;
    };


private:
    WeightSet<T> sel_weights;
    Neighborhood associations;

    std::vector<int> ranks_aux;
    std::vector<T> fitness_aux;
};
