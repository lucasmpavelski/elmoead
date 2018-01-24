#pragma once

#include <cstring>
#include <unordered_set>

#include "../aux.h"

#include "population.h"

class NdArchive {
    Population solutions;
    size_t act_size, max_size;
    double eps;
    std::unordered_set<size_t> nd_set;

public:
    NdArchive(MOProblem& prob,
              const size_t max_samples = 100,
              const double eps = 0.001)
        : solutions(1, prob)
        , act_size(0)
        , max_size(max_samples)
        , eps(eps)
        , nd_set()
    {
    }

    NdArchive(const Population& init_sols,
              const size_t max_samples = 100,
              const double eps = 0.001)
        : solutions(init_sols)
        , act_size(init_sols.size())
        , max_size(max_samples)
        , eps(eps)
        , nd_set()
    {
        //solutions.resize(max_samples);
    }

    size_t size()       const { return act_size;      }
    size_t getMaxSize() const { return max_size;      }
    double getEps()     const { return eps;           }
    size_t getNoNd()    const { return nd_set.size(); }
    const Population& all() const { return solutions; }

    void setEps(double neps) { eps = neps; }
    void setProblem(MOProblem* p) { solutions.setProblem(p); }

    size_t getNdPopulation(Population& pop) const
    {
        pop.resize(getNoNd());
        size_t i = 0;
        for (const auto& idx : nd_set)
            pop[i++] = solutions[idx];
        return getNoNd();
    }

    bool contain(const Solution& sol) const;
    bool add(const Solution& sol);
    size_t add(const Population& pop);
};
