#pragma once

#include <assert.h>
#include <stdbool.h>
#include <float.h>

#include "immintrin.h"

#include "../aux.h"
#include "weights.h"

#define PBI_THETA 5.0

typedef enum
{
    WEIGHTED_SUM,
    TCHEBYCHEFF,
    INVERTED_TCHEBYCHEFF,
    AUGMENTED_TCHEBYCHEFF,
    NORMALIZED_TCHEBYCHEFF,
    NORMALIZED_INVERTED_TCHEBYCHEFF,
    PENALTY_BOUNDARY_INTERSECTION,
    ACHIEVEMENT_SCALARIZING,
    AUGMENTED_ACHIEVEMENT_SCALARIZING
} AggregationFunctionType;

inline AggregationFunctionType str2AggrFunc(const std::string& str) {
    if (str == "WS"            || str == "WEIGHTED_SUM"                     ) return WEIGHTED_SUM;
    if (str == "TCHE"          || str == "TCHEBYCHEFF"                      ) return TCHEBYCHEFF;
    if (str == "INV_TCHE"      || str == "INVERTED_TCHEBYCHEFF"             ) return INVERTED_TCHEBYCHEFF;
    if (str == "AUG_TCHE"      || str == "AUGMENTED_TCHEBYCHEFF"            ) return AUGMENTED_TCHEBYCHEFF;
    if (str == "NORM_TCHE"     || str == "NORMALIZED_TCHEBYCHEFF"           ) return NORMALIZED_TCHEBYCHEFF;
    if (str == "NORM_INV_TCHE" || str == "NORMALIZED_INVERTED_TCHEBYCHEFF"  ) return NORMALIZED_INVERTED_TCHEBYCHEFF;
    if (str == "PBI"           || str == "PENALTY_BOUNDARY_INTERSECTION"    ) return PENALTY_BOUNDARY_INTERSECTION;
    if (str == "AUG_AS"        || str == "AUGMENTED_ACHIEVEMENT_SCALARIZING") return AUGMENTED_ACHIEVEMENT_SCALARIZING;

    throw_assert(false,
        "Aggregation function must be one of: \n"
        "(WS           | WEIGHED_SUM                      ) | \n"
        "(TCHE         | TCHEBYCHEFF                      ) | \n"
        "(INV_TCHE     | INVERTED_TCHEBYCHEFF             ) | \n"
        "(AUG_TCHE     | AUGMENTED_TCHEBYCHEFF            ) | \n"
        "(NORM_TCHE    | NORMALIZED_TCHEBYCHEFF           ) | \n"
        "(NORM_INV_TCHE| NORMALIZED_INVERTED_TCHEBYCHEFF  ) | \n"
        "(PBI          | PENALTY_BOUNDARY_INTERSECTION    ) | \n"
        "(AUG_AS       | AUGMENTED_ACHIEVEMENT_SCALARIZING)\n"
    );

    return TCHEBYCHEFF;
}

template<AggregationFunctionType type=TCHEBYCHEFF>
class AggrFunc
{    
    static constexpr double eps = 0.000001;

    static double epsReplace(const double w) {
        return (w < eps) ? eps : w;
    };

public:

    double apply(const double objs[], const int no_objs, const double w[],
                 const double ideal[]=nullptr, const double nadir[]=nullptr) {
        return tchebycheff(objs, no_objs, w, ideal);
    };

    AggrFunc() {};

    static double weightedSum(const double* obj, const int no_objs,
                              const double w[]) {
        double r = obj[0] * w[0];
        for (int i = 1; i < no_objs; ++i)
        {
            //double w = (weights[n] <= 0) ? 0.00001 : weights[n];
            r += obj[i] * w[i];
        }
        return r;
    };

    static double tchebycheff(const double obj[], const int no_objs,
                              const double weights[], const double ideal[]) {
        double max_fun = fabs(obj[0] - ideal[0]) *
                (weights[0] == 0 ? eps : weights[0]);
        for (int i = 1; i < no_objs; i++)
        {
            const double feval = fabs(obj[i] - ideal[i]) *
                    (weights[i] == 0? eps : weights[i]);
            if (feval > max_fun)
                max_fun = feval;
        }
        return max_fun;
    };

    static double invertedTchebycheff(const double obj[], const int no_objs,
                                      const double w[], const double ideal[]) {
        double max_fun = (obj[0] - ideal[0]) / (w[0] == 0 ? eps : w[0]);
        for (int i = 1; i < no_objs; i++)
        {
            const double wi = w[i] == 0 ? eps : w[i];
            const double feval = (obj[i] - ideal[i]) / wi;
            if (feval > max_fun)
                max_fun = feval;
        }
        return max_fun;
    };

    static double augmentedTchebycheff(const double obj[], const int no_objs,
                                       const double weights[],
                                       const double ideal[]) {
        const double weighted_sum = weightedSum(obj, no_objs, weights);
        const double tche = invertedTchebycheff(obj, no_objs, weights, ideal);
        return tche + 0.01 * weighted_sum;
    };

    static double normalizedTchebycheff(const double obj[], const int no_objs,
                                        const double w[], const double ideal[],
                                        const double nadir[]) {
        const double w0 = w[0] == 0 ? eps : w[0];
        double max_fun = (obj[0] - ideal[0]) * w0 / (nadir[0] - ideal[0]);
        for (int i = 1; i < no_objs; ++i)
        {
            const double wi = w[i] == 0 ? eps : w[i];
            const double fi = (obj[i] - ideal[i]) * wi / (nadir[i] - ideal[i]);
            if (fi > max_fun)
                max_fun = fi;
        }
        return max_fun;
    };

    static double normalizedInvertedTchebycheff(const double obj[],
                                                const int no_objs,
                                                const double w[],
                                                const double ideal[],
                                                const double nadir[]) {
        const double w0 = w[0] == 0 ? eps : w[0];
        double max_fun = (obj[0] - ideal[0]) / ((nadir[0] - ideal[0]) * w0);
        for (int i = 1; i < no_objs; ++i)
        {
            const double wi = w[i] == 0 ? eps : w[i];
            const double feval = (obj[i] - ideal[i]) /
                    ((nadir[i] - ideal[i]) * wi);
            if (feval > max_fun)
                max_fun = feval;
        }
        return max_fun;
    };

    static double penaltyBoundaryIntersection(const double obj[],
                                              const int no_objs,
                                              const double w[],
                                              const double ideal[],
                                              const double penalty=5.0) {
        double d1 = 0.0, w_norm = 0.0;
        for (int i = 0; i < no_objs; ++i)
        {
           d1 += (obj[i] - ideal[i]) * w[i];
           w_norm += w[i] * w[i];
        }
        d1 = fabs(d1) / sqrt(w_norm);
        double d2 = 0.0; // distance to the search direction norm
        for (int i = 0; i < no_objs; ++i)
        {
           const double d = obj[i] - (ideal[i] + d1 * (w[i]));
           d2 += d * d;
        }
        d2 = sqrt(d2);
        return d1 + penalty * d2;
    };

    static double augmentedAchievementScalarizing(const double obj[],
                                                  const int no_objs,
                                                  const double w[],
                                                  const double ideal[],
                                                  const double penalty) {
        static const double eps = 0.00001;
        double max_fun = (obj[0] - ideal[0]) /
                (w[0] == 0 ? eps : w[0]);
        double sum = max_fun;
        for (int i = 1; i < no_objs; i++)
        {
            const double wi = (w[i] == 0) ? eps : w[i];
            const double feval = (obj[i] - ideal[i]) / wi;
            sum += feval;
            if (feval > max_fun)
                max_fun = feval;
        }
        return max_fun + penalty * sum;
    };
};

template<>
struct AggrFunc<WEIGHTED_SUM>
{
    double apply(const double obj[], const int no_objs, const double w[],
                 const double ideal[]=nullptr, const double nadir[]=nullptr) {
        return AggrFunc<>::weightedSum(obj, no_objs, w);
    };
};

template<>
struct AggrFunc<INVERTED_TCHEBYCHEFF>
{
    double apply(const double obj[], const int no_objs, const double w[],
                 const double ideal[], const double nadir[]=nullptr) {
        return AggrFunc<>::invertedTchebycheff(obj, no_objs, w, ideal);
    };
};

template<>
struct AggrFunc<NORMALIZED_TCHEBYCHEFF>
{
    double apply(const double obj[], const int no_objs, const double w[],
                 const double ideal[]=nullptr, const double nadir[]=nullptr) {
        return AggrFunc<>::normalizedTchebycheff(obj, no_objs, w, ideal,
                                                         nadir);
    };
};

template<>
struct AggrFunc<NORMALIZED_INVERTED_TCHEBYCHEFF>
{
    double apply(const double obj[], const int no_objs, const double w[],
                 const double ideal[]=nullptr, const double nadir[]=nullptr) {
        return AggrFunc<>::normalizedInvertedTchebycheff(obj, no_objs, w, ideal,
                                                         nadir);
    };
};

template<>
struct AggrFunc<AUGMENTED_TCHEBYCHEFF>
{
    double apply(const double obj[], const int no_objs, const double w[],
                 const double ideal[]=nullptr, const double nadir[]=nullptr) {
        return AggrFunc<>::augmentedTchebycheff(obj, no_objs, w, ideal);
    };
};

template<>
struct AggrFunc<PENALTY_BOUNDARY_INTERSECTION>
{
    const double penalty;
    AggrFunc(const double penalty=5.0) : penalty(penalty) {};

    double apply(const double obj[], const int no_objs, const double w[],
                 const double ideal[]=nullptr, const double nadir[]=nullptr) {
        return AggrFunc<>::penaltyBoundaryIntersection(obj, no_objs, w, ideal,
                                                       penalty);
    };
};

template<>
struct AggrFunc<AUGMENTED_ACHIEVEMENT_SCALARIZING>
{
    const double penalty;
    AggrFunc(const double penalty=5.0) : penalty(penalty) {};

    double apply(const double obj[], const int no_objs, const double w[],
                 const double ideal[]=nullptr, const double nadir[]=nullptr) {
        return AggrFunc<>::augmentedAchievementScalarizing(obj, no_objs, w,
                                                           ideal, penalty);
    };
};

double aggregationFunction(const AggregationFunctionType type, 
        const double objs[], const int no_objs, const double weights[], 
        const double ideal[], const double nadir[]);

double weightedSum(const double* obj, const int no_objs,
                   const double weights[]);

double tchebycheff(const double obj[], const int no_objs,
                   const double weights[], const double ideal[]);

double invertedTchebycheff(const double obj[], const int no_objs,
                           const double weights[], const double ideal[]);

double augmentedTchebycheff(const double obj[], const int no_objs,
                            const double weights[], const double ideal[]);

double normalizedTchebycheff(const double obj[], const int no_objs,
                             const double weights[], const double ideal[],
                             const double nadir[]);

double normalizedInvertedTchebycheff(const double obj[], const int no_objs,
                                     const double weights[],
                                     const double ideal[], const double nadir[]);

double penaltyBoundaryIntersection(const double obj[], const int no_objs,
                                   const double weights[], const double ideal[]);

double augmentedAchievementScalarizing(const double obj[], const int no_objs,
                                       const double weights[],
                                       const double ideal[]);

//template <AggregationFunctionType aft, typename T=double>
template <typename T=double>
struct AggrProblem
{
    using vector = Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>;

    AggrProblem(const AggregationFunctionType aggr_func, const size_t no_objs) :
        aggr_func (aggr_func),
        no_objs   (no_objs),
        ideal     (no_objs),
        nadir     (no_objs),
        norm_scale(no_objs),
        norm_trans(no_objs)
    {
        resetReferencePoints();
    };

    inline T const* idealPtr() const { return ideal.data(); };
    inline T const* nadirPtr() const { return nadir.data(); };

    inline T* idealPtr() { return ideal.data(); };
    inline T* nadirPtr() { return nadir.data(); };

    inline void resetReferencePoints() {
        const T inf = numeric_limits<T>::infinity();
        ideal.fill(inf);
        nadir.fill(-inf);
        norm_scale.fill(1.0);
        norm_trans.fill(0.0);
    };

    inline bool updateReferencePoints(T const objs[]) {
        bool changed = false;
        for (int i = 0; i < no_objs; ++i)
        {
            const T v = objs[i];
            if (v < ideal[i])
            {
                ideal[i] = v;
                changed = true;
            }
            if (v > nadir[i])
            {
                nadir[i] = v;
                changed = true;
            }
        }
        if (changed)
            updateNormCaches();
        return changed;
    };

    inline T eval(const T objs[], const WeightSet<T>& w_set,
                  const int idx) const {
        return eval(objs, w_set[idx]);
    }

    inline T eval(const T objs[], const T w[]) const {
        return //aggr_func.apply(objs, no_objs, w, ideal.data(), nadir.data());

               aggregationFunction(aggr_func, objs, no_objs, w, ideal.data(),
                                   nadir.data());
    }

    inline T affinity(const T objs[], T const w[]) const {
        T norm_objs[no_objs], proj[no_objs];
        for (int i = 0; i < no_objs; ++i)
            norm_objs[i] = objs[i] * norm_scale[i] + norm_trans[i];
        projection(norm_objs, w, no_objs, proj);
        return sqrEuclideanDist(norm_objs, proj, no_objs);
    }

    static bool aggregationFunctionUseNadir(const AggregationFunctionType af) {
        return af == NORMALIZED_TCHEBYCHEFF ||
               af == NORMALIZED_INVERTED_TCHEBYCHEFF;
    };

    const AggregationFunctionType aggr_func;
    const size_t no_objs;
    vector ideal;
    vector nadir;
    vector norm_scale;
    vector norm_trans;

    friend std::ostream& operator<<(std::ostream& os, const AggrProblem& ap) {
        return os << "AggrProblem:\n"
                  //<< "  aggr_func: " << ap.aggr_func << "\n";
                  << "  no_objs: " << ap.no_objs << "\n"
                  << "  current_ideal: " << ap.ideal.transpose() << "\n"
                  << "  current_nadir: " << ap.nadir.transpose();
    };

private:
    void updateNormCaches() {
        norm_scale = (nadir - ideal).cwiseInverse();
        norm_trans = - ideal.cwiseProduct(norm_scale);
    };
};
