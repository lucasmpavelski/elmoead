#pragma once

#include <vector>
#include <algorithm>

#include "rand.h"
#include "error.h"

enum SamplingMethod {
    NO_SAMPLING,
    RANDOM_SAMPLING, 
    LATIN_HYP_SAMPLING
};

void randomSampling(const int no_rows, const int no_cols,
        const double* min_col, const double* max_col, double* output);

void latinHypercubeSampling(const int no_rows, const int no_cols,
        const double* min_col, const double* max_col, double* output);

void sample(const SamplingMethod method, const int no_rows, 
        const int no_cols, const double* min_col, const double* max_col, 
        double* output);

template<SamplingMethod sm>
void sample(const int no_rows,  const int no_cols, const double* min_col,
            const double* max_col, double* output);
