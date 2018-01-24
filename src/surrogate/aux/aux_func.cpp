#include "aux_func.h"

void sortIdx(double x[], int idx[], const int n, const int m)
{
    int i;
    for (i = 0; i < m; i++)
    {
        double min_x = x[i];
        int min_x_idx = i;
        int j;
        for (j = i + 1; j < n; j++)
        {
            if (x[j] < min_x)
            {
                min_x = x[j];
                min_x_idx = j;
            }
        }

        if (min_x_idx != i)
        {
            std::swap(x[i], x[min_x_idx]);
            std::swap(idx[i], idx[min_x_idx]);
        }
    }
}

double sqr(const double& base)
{
    return base * base;
}

double intpow(double base, int exp)
{
    double result = 1.0;
    while (exp != 0)
    {
        if ((exp & 1) == 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }

    return result;
}

int binomial(const int n, const int k)
{
    int new_k = (k > n - k) ? n - k : k;
    int acc = 1;
    int i;
    for (i = 1; i < new_k + 1; ++i)
    {
        acc *= n - (new_k - i);
        acc = acc / i;
    }
    return acc;
}

double sqrEuclideanDist(const double a[], const double b[], const int n)
{
    double v = (a[0] - b[0]) * (a[0] - b[0]);
    int i;
    for (i = 1; i < n; i++)
        v += (a[i] - b[i]) * (a[i] - b[i]);
    return v;
}

double euclideanDist(const double a[], const double b[], const int n)
{
    return sqrt(sqrEuclideanDist(a, b, n));
}

double* orders_aux;
static int orders_compar(const void *a, const void *b)
{
    const double aa = orders_aux[*((const int*)a)];
    const double bb = orders_aux[*((const int*)b)];
    if (aa < bb)
        return -1;
    if (aa > bb)
        return 1;
    // if (aa == bb)
    return 0;
}

int* ranks(const double x[], const int n, int r[])
{
    if (r == NULL)
        r = Malloc(int, n);
    orders_aux = (double*)x;
    int i;
    for (i = 0; i < n; ++i)
        r[i] = i;
    qsort(r, n, sizeof(int), orders_compar);
    orders_aux = NULL;
    return r;
}


double norm(const double x[], const int n)
{
    double r = x[0] * x[0];
    int i;
    for (i = 1; i < n; ++i)
        r += x[i] * x[i];
    return sqrt(r);
}

double* projection(const double v1[], const double v2[], const int n,
                   double ep[])
{
    if (ep == NULL)
        ep = Malloc(double, n);

    double v2_norm = v2[0] * v2[0];
    double inner = v1[0] * v2[0];

    int i;
    for (i = 1; i < n; ++i)
    {
        inner += v1[i] * v2[i];
        v2_norm += v2[i] * v2[i];
    }

    const double s = inner / v2_norm;
    for (i = 0; i < n; ++i)
        ep[i] = s * v2[i];

    return ep;
}

double stddev(const double x[], const int n)
{
    const double m = mean(x, n);
    return stdWithMean(x, n, m);
}

double stdWithMean(const double x[], const int n, const double m)
{
    double sum = (x[0] - m) * (x[0] - m);
    int i;
    for (i = 1; i < n; ++i)
        sum += (x[i] - m) * (x[i] - m);
    return sqrt(sum / n);
}

void colwiseMean(const double* const x[], const int n, const int no_cols,
                 double cm[])
{
    for (int i = 0; i < no_cols; ++i)
        cm[i] = x[0][i];
    for (int i = 1; i < n; ++i)
    {
        for (int j = 0; j < no_cols; ++j)
            cm[j] += x[i][j];
    }
    for (int i = 0; i < no_cols; ++i)
        cm[i] = cm[i] / n;
}

void colwiseMeanStd(const double* const x[], const int n, const int no_cols,
                    double cm[], double cstd[])
{
    bool free_mean = cm == NULL;
    if (cm == NULL)
        cm = Malloc(double, n);

    colwiseMean(x, n, no_cols, cm);

    for (int i = 0; i < no_cols; ++i)
        cstd[i] = (x[0][i] - cm[i]) * (x[0][i] - cm[i]);
    for (int i = 1; i < n; ++i)
    {
        for (int j = 0; j < no_cols; ++j)
            cstd[j] += (x[i][j] - cm[j]) * (x[i][j] - cm[j]);
    }
    for (int i = 0; i < no_cols; ++i)
        cstd[i] = sqrt(cstd[i] / n);

    if (free_mean)
        free(cm);
}
