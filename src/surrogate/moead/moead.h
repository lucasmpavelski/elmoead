#pragma once

#include "../aux.h"
#include "../mo.h"
#include "../obj_dec.h"
#include "../var.h"

class MOEAD {
public:
    using dmatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using dvector = Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>;
    using ivector = Eigen::Matrix<int, 1, Eigen::Dynamic, Eigen::RowMajor>;

    const size_t no_objs;
    const PMOperator pm_operator;
    const SBXoverOperator sbx_operator;
    const size_t pop_size;
    AggrProblem<double> aggr_prob;
    WeightSet<double> weights;
    Neighborhood neighbors;
    ivector P;
    std::vector<double> fitness_aux;

    MOEAD(const WeightSet<double>& weights,
          const size_t no_neighbors,
          const AggregationFunctionType aggr_func,
          const PMOperator pm_operator,
          const SBXoverOperator sbx_operator)
        : no_objs(weights.dimension())
        , pm_operator(pm_operator)
        , sbx_operator(sbx_operator)
        , pop_size(weights.size())
        , weights(weights)
        , neighbors(weights, weights, no_neighbors)
        , aggr_prob(aggr_func, weights.dimension())
        , P(weights.size())
        , fitness_aux(weights.size())
    {
        std::iota(P.data(), P.data() + pop_size, 0);
    }

    virtual ~MOEAD(){}

    AggregationFunctionType aggrFunc() const { return aggr_prob.aggr_func; }

    void initReferencePoints(const Population& pop);
    bool updateReferencePoints(const Solution& sol);

    virtual int initPopSize() const { return pop_size; }
    int popSize() const { return pop_size; }

    virtual void optimize(const int no_evals, MOProblem& prob, Population& pop,
                          PopulationCallback callback);

    virtual void evolution(const Population& pop, const int idx, Solution& result,
                           const int candidates_idxs[], const int no_candidates);

    virtual int update(Population& pop, const Solution& ind,
                       const int* candidates_idxs, const int no_candidates);

    inline double subProbEval(const Population& pop, const int& idx) {
        return subProbEval(pop[idx], idx);
    }

    inline double subProbEval(const Solution& sol, const int& idx) {
        //return weights2.eval(aggr_func, idx, sol.objs, ideal.get(), nadir.get());
        return subProbEval(sol.objs, idx);
    };

    inline double subProbEval(const double objs[], const int& idx) {
        //return aggregationFunction(aggr_func, objs, no_objs, weights[idx],
        //                           ideal.get(), nadir.get());
        return aggr_prob.eval(objs, weights[idx]);
    };

    void sortPopulationByWeights(Population& pop);


    friend std::ostream& operator<<(std::ostream& os, const MOEAD& moead)
    {
        os << "MOEAD :\n"
           << "  pop_size: " << moead.pop_size << "\n"
           << "  no_objs: " << moead.no_objs << "\n"
           << "  no_neighbors: " << moead.neighbors.no_neighbors << "\n"
           << moead.aggr_prob
           << moead.sbx_operator
           << moead.pm_operator;
    }
};
