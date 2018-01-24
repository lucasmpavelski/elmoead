#pragma once

#include "moead/moead.h"
#include "moead/moead_de.h"
#include "moead/moead_dra.h"
#include "moead/moead_stm.h"

inline MOEAD* getMoeadAlg(std::string type,
                             WeightSet<> weights,
                             size_t no_neighbors,
                             AggregationFunctionType aggr_func,
                             double real_b,
                             size_t update_limit,
                             SBXoverOperator sbx,
                             OperatorSelection<DEOperator> operators,
                             PMOperator pmo)
{
    if (type == "MOEA/D")
        return new MOEAD(weights, no_neighbors, aggr_func, pmo, sbx);
    if (type == "MOEA/D-DE")
        return new MOEAD_DE(weights, no_neighbors, aggr_func, real_b, update_limit, pmo, operators);
    if (type == "MOEA/D-DRA")
        return new MOEAD_DRA(weights, no_neighbors, aggr_func, real_b, update_limit, pmo, operators);
    if (type == "MOEA/D-STM")
        return new MOEAD_STM(weights, no_neighbors, aggr_func, real_b, pmo, operators);
    return nullptr;
}
