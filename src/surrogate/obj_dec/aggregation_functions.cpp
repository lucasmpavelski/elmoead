#include "aggregation_functions.h"

double aggregationFunction(const AggregationFunctionType type,
                           const double objs[], const int no_objs,
                           const double weights[], const double ideal[],
                           const double nadir[])
{
    switch(type)
    {
        case WEIGHTED_SUM:
            return weightedSum(objs, no_objs, weights);
        break;
        case TCHEBYCHEFF:
            return tchebycheff(objs, no_objs, weights, ideal);
        break;
        case INVERTED_TCHEBYCHEFF:
        case ACHIEVEMENT_SCALARIZING:
            return invertedTchebycheff(objs, no_objs, weights, ideal);
        break;
        case AUGMENTED_TCHEBYCHEFF:
            return augmentedTchebycheff(objs, no_objs, weights, ideal);
        break;
        case NORMALIZED_TCHEBYCHEFF:
            return normalizedTchebycheff(objs, no_objs, weights, ideal, nadir);
        break;
        case NORMALIZED_INVERTED_TCHEBYCHEFF:
            return normalizedInvertedTchebycheff(objs, no_objs, weights, ideal, nadir);
        break;
        case PENALTY_BOUNDARY_INTERSECTION:
            return penaltyBoundaryIntersection(objs, no_objs, weights, ideal);
        break;
        case AUGMENTED_ACHIEVEMENT_SCALARIZING:
            return augmentedAchievementScalarizing(objs, no_objs, weights, ideal);
        break;
        default:
            throw_assert(false, "unknown aggregation function type");
        break;
    }
    return nan("");
}

double weightedSum(const double* obj, const int no_objs,
                   const double weights[])
{
    double r = obj[0] * weights[0];
    for (int i = 1; i < no_objs; ++i)
    {
        //double w = (weights[n] <= 0) ? 0.00001 : weights[n];
        r += obj[i] * weights[i];
    }
    return r;
}

double tchebycheff(const double obj[], const int no_objs,
                   const double weights[], const double ideal[])
{
    static const double eps = 0.00001;
    double max_fun = fabs(obj[0] - ideal[0]) * (weights[0] == 0? eps : weights[0]);
    for (int i = 1; i < no_objs; i++)
    {
        const double feval = fabs(obj[i] - ideal[i]) * (weights[i] == 0? eps : weights[i]);
        if (feval > max_fun)
            max_fun = feval;
    }
    return max_fun;
}

double invertedTchebycheff(const double obj[], const int no_objs,
                           const double weights[], const double ideal[])
{
    static const double eps = 0.00001;
    double max_fun = (obj[0] - ideal[0]) / (weights[0] == 0 ? eps : weights[0]);
    for (int i = 1; i < no_objs; i++)
    {
        const double w = weights[i] == 0 ? eps : weights[i];
        const double feval = (obj[i] - ideal[i]) / w;
        if (feval > max_fun)
            max_fun = feval;
    }
    return max_fun;
}

double augmentedTchebycheff(const double obj[], const int no_objs,
                            const double weights[], const double ideal[])
{
    const double weighted_sum = weightedSum(obj, no_objs, weights);
    const double tche = invertedTchebycheff(obj, no_objs, weights, ideal);
    return tche + 0.01 * weighted_sum;
}

double normalizedTchebycheff(const double obj[], const int no_objs,
                             const double weights[], const double ideal[],
                             const double nadir[])
{
    static const double eps = 0.00001;
    const double w0 = weights[0] == 0 ? eps : weights[0];
    double max_fun = (obj[0] - ideal[0]) * w0 / (nadir[0] - ideal[0]);
    for (int i = 1; i < no_objs; ++i)
    {
        const double wi = weights[i] == 0 ? eps : weights[i];
        const double fi = (obj[i] - ideal[i]) * wi / (nadir[i] - ideal[i]);
        if (fi > max_fun)
            max_fun = fi;
    }
    return max_fun;
}

double normalizedInvertedTchebycheff(const double obj[], const int no_objs,
                                     const double weights[],
                                     const double ideal[], const double nadir[])
{
    static const double eps = 0.00001;
    const double w0 = weights[0] == 0 ? eps : weights[0];
    double max_fun = (obj[0] - ideal[0]) / ((nadir[0] - ideal[0]) * w0);
    for (int i = 1; i < no_objs; ++i)
    {
        const double w = weights[i] == 0 ? eps : weights[i];
        const double feval = (obj[i] - ideal[i]) / ((nadir[i] - ideal[i]) * w);
        if (feval > max_fun)
            max_fun = feval;
    }
    return max_fun;
}

double penaltyBoundaryIntersection(const double obj[], const int no_objs,
                                   const double weights[], const double ideal[])
{
    // difference beween current point and reference point
    double d1 = 0.0, w_norm = 0.0;
    for (int i = 0; i < no_objs; ++i)
    {
        d1 += (obj[i] - ideal[i]) * weights[i];
        w_norm += weights[i] * weights[i];
    }
    d1 = fabs(d1) / sqrt(w_norm);

    // distance to the search direction norm
    double d2 = 0.0;
    for (int i = 0; i < no_objs; ++i)
    {
        const double d = obj[i] - (ideal[i] + d1 * (weights[i]));
        d2 += d * d;
    }
    d2 = sqrt(d2);

    return d1 + PBI_THETA * d2;
}

double augmentedAchievementScalarizing(const double obj[], const int no_objs,
                                       const double weights[], const double ideal[])
{
    static const double eps = 0.00001;
    double max_fun = (obj[0] - ideal[0]) / (weights[0] == 0 ? eps : weights[0]);
    double sum = max_fun;
    for (int i = 1; i < no_objs; i++)
    {
        const double w = weights[i] == 0 ? eps : weights[i];
        const double feval = (obj[i] - ideal[i]) / w;
        sum += feval;
        if (feval > max_fun)
            max_fun = feval;
    }
    return max_fun + PBI_THETA * sum;
}
