#include "sampling_methods.h"

void colsRange(const int no_cols, const double* min_col, 
        const double* max_col, double* range)
{
    for (int d = 0; d < no_cols; d++)
        range[d] = max_col[d] - min_col[d];   
}

void sample(const SamplingMethod method, const int no_rows, 
        const int no_cols, const double* min_col, const double* max_col, 
        double* output)
{
    switch(method) 
    {
        case NO_SAMPLING :
        break;
        case RANDOM_SAMPLING :
            randomSampling(no_rows, no_cols, min_col, max_col, output);
        break;
        case LATIN_HYP_SAMPLING :
            latinHypercubeSampling(no_rows, no_cols, min_col, max_col, output);
        break;
        default :
            throw_assert(false, "Unknown sampling method");
        break;
    }
}

void randomSampling(const int no_rows, const int no_cols, 
        const double* min_col, const double* max_col, double* output)
{
    for (int i = 0; i < no_rows; ++i)
    {
        for (int d = 0; d < no_cols; ++d)
            output[i * no_cols + d] = RNG::realUniform(min_col[d], max_col[d]);
    }
}

void latinHypercubeSampling(const int no_samples, const int no_dims,
        const double* min_col, const double* max_col, double* output)
{
    std::vector<double> grid(no_dims);
    for (int d = 0; d < no_dims; d++)
        grid[d] = (max_col[d] - min_col[d]) / no_samples;

    std::vector<int> a(no_samples);
    std::iota(a.begin(), a.end(), 0); // a = 0, 1, ..., no_rows

    for (int d = 0; d < no_dims; d++)
    {
        std::shuffle(a.begin(), a.end(), RNG::engine);

        for (int i = 0; i < no_samples; i++)
        {
            const double lower_bound = min_col[d] + grid[d] * i;
            const double upper_bound = min_col[d] + grid[d] * (i + 1);
            const int idx = a[i] * no_dims + d;
            output[idx] = RNG::realUniform(lower_bound, upper_bound);
        }
    }
}

template<>
void sample<NO_SAMPLING>(const int no_rows,  const int no_cols,
                         const double* min_col, const double* max_col,
                         double* output)
{}

template<>
void sample<RANDOM_SAMPLING>(const int no_rows,  const int no_cols,
                         const double* min_col, const double* max_col,
                         double* output)
{
    randomSampling(no_rows, no_cols, min_col, max_col, output);
}

template<>
void sample<LATIN_HYP_SAMPLING>(const int no_rows,  const int no_cols,
                         const double* min_col, const double* max_col,
                         double* output)
{
    latinHypercubeSampling(no_rows, no_cols, min_col, max_col, output);
}

