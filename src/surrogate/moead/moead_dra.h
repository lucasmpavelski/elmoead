#pragma once

#include "moead_de.h"
#include "stm.h"

class MOEAD_DRA : public MOEAD_DE
{
    dvector utility;
    dmatrix previous_objs;
    ivector candidates;

public:

    MOEAD_DRA(const WeightSet<double>& weights, const int no_neighbors,
              const AggregationFunctionType aggr_func,
              const double real_b, const int update_limit,
              const PMOperator pm_operator,
              OperatorSelection<DEOperator>& operators) :
        MOEAD_DE(weights, no_neighbors, aggr_func,
                 real_b, update_limit, pm_operator, operators),
        utility(pop_size),
        selected(pop_size / 5),
        candidates(pop_size - 2),
        previous_objs(pop_size, no_objs)
    {
        utility.fill(1.0);
    }

    virtual ~MOEAD_DRA() {}

    virtual int initPopSize() const { return pop_size + pop_size / 5; }

    virtual void optimize(const int no_evals, MOProblem& prob, Population& pop,
                          PopulationCallback callback);

    // DRA specific
    ivector selected; // tour selection "return" TODO: refactor
    void rememberObjs(const Population& pop);
    void tourSelection(const int depth);
    void compUtility(const Population& pop);

    friend std::ostream& operator<<(std::ostream& os, const MOEAD_DRA& moead) {
        operator<<(os, static_cast<const MOEAD_DE&>(moead));
        os << "MOEAD_DRA:\n";
        return os;
    }
};
