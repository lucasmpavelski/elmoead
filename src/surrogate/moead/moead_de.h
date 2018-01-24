#pragma once

#include "../aux.h"
#include "../var.h"
#include "../aos.h"

#include "moead.h"

class MOEAD_DE : public MOEAD {
private:
    const int update_limit;
    const double real_b;
    OperatorSelection<DEOperator> operators;

public:
    MOEAD_DE(const WeightSet<double>& weights, const int no_neighbors,
        const AggregationFunctionType aggr_func,
        const double real_b, const int update_limit,
        const PMOperator pm_operator,
        OperatorSelection<DEOperator>& operators)
        : MOEAD(weights, no_neighbors, aggr_func, pm_operator,
              SBXoverOperator())
        , real_b(real_b)
        , update_limit(update_limit)
        , operators(operators)
    {
    }

    virtual int initPopSize() const override { return pop_size; }
    double realB() const { return real_b; }

    virtual void optimize(const int no_evals, MOProblem& prob, Population& pop,
        PopulationCallback callback) override;

    virtual void evolution(const Population& pop, const int idx,
        Solution& result, const int candidates_idxs[],
        const int no_candidates);

    virtual int update(Population& pop, const Solution& ind,
        const int candidates_idxs[], const int no_candidates);

    friend std::ostream& operator<<(std::ostream& os, const MOEAD_DE& moead)
    {
        operator<<(os, static_cast<const MOEAD&>(moead));
        os << "MOEAD_DE:\n";
        os << "  real_b: " << moead.real_b << endl;
        os << "  update_limit: " << moead.update_limit << endl;
        os << moead.operators;
        return os;
    }
};
