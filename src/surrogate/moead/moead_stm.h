#pragma once

#include "moead_dra.h"
#include "stm.h"

class MOEAD_STM : public MOEAD_DRA
{
    STMData stm_data;

public:
    MOEAD_STM(const WeightSet<double>& weights, const int no_neighbors,
              const AggregationFunctionType aggr_func,
              const double real_b, const PMOperator pm_operator,
              OperatorSelection<DEOperator>& operators,
              const PopulationCallback cb = NullCallback) :
        MOEAD_DRA(weights, no_neighbors, aggr_func,
                  real_b, 0, pm_operator, operators),
        stm_data(initPopSize(), pop_size)
    {}

    virtual ~MOEAD_STM() {}

    virtual void optimize(const int no_evals, MOProblem& prob, Population& pop,
                          PopulationCallback callback);

    friend std::ostream& operator<<(std::ostream& os, const MOEAD_STM& moead) {
        operator<<(os, static_cast<const MOEAD_DE&>(moead));
        os << "MOEAD_STM:\n";
        return os;
    }
};
