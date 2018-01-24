#include "../surrogate/moead_algs.h"
#include "../surrogate/benchmark.h"


TEST(MOEADTests,Optimize)
{
    MOProblem* prob = getBenchmark("ZDT1", 30);
    const int no_partitions = 299;
    const int no_neighbors = 60;
    const AggregationFunctionType aggr_func = TCHEBYCHEFF;

/*
    MOEAD moead(prob->no_objs, no_partitions, no_neighbors,
                aggr_func,
                PMOperator(prob->no_vars));
*/
    const string de_operator = "DE/RAND/1/BIN";
    const double de_cr = 1.0;
    const double de_f = 0.5;

    /* de operator types */
    std::vector<DEOperator> de_ops;
    std::stringstream ss(de_operator);
    std::transform(std::istream_iterator<string>(ss),
        std::istream_iterator<string>(),
        std::back_inserter(de_ops),
        [de_cr, de_f](const string tp) {
            return DEOperator(str2DEType(tp), de_cr, de_f);
        });
    /* de operators adaptation */
    OperatorSelection<DEOperator> operators(de_ops);

    MOEAD_DE moead(WeightSet<>(prob->no_objs, no_partitions), no_neighbors,
                   aggr_func, 0.9, 6,
                   PMOperator(static_cast<int>(prob->no_vars), 20.0),
                   operators
        //          ,MOEAD::PlotCallback()
                   );

/*
    MOEAD_STM moead(prob->no_objs, no_partitions, no_neighbors,
                   aggr_func, 0.9,
                   PMOperator(prob->no_vars),
                   DEOperator(DEOperator::Type::DE_1_BIN)
  //                ,MOEAD::PlotCallback()
                   );
*/

    Population pop(moead.initPopSize(), *prob, LATIN_HYP_SAMPLING);
    pop.evaluateAll();
    moead.optimize(6000 - pop.size(), *prob, pop, NullCallback);
}
