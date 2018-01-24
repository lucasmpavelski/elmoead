#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>


template<typename Itr1, typename Itr2, typename T>
bool contain(Itr1 beg, Itr2 end, const T& val)
{
    return std::any_of(beg, end, [&val](const T& v) { return v == val; });
}

template<typename T, typename Itr1, typename Itr2>
std::pair<T,T> circularRange(Itr1 b, Itr2 e, const T lb, const T ub)
{
    const int n = std::distance(b, e);
    std::sort(b, e);
    int li = 0, ui = n - 1;
    const T first = *(b + li), last = *(b + ui);
    T gap = (ub - last) + (lb - first);
    for (int i = 0; i < n - 1; i++)
    {
        T other_gap = *(b + i + 1) - *(b + i);
        if (other_gap > gap)
        {
            gap = other_gap;
            li = i + 1;
            ui = i;
        }
    }
    return {*(b + li), *(b + ui)};
}

template <typename Itr>
inline typename std::iterator_traits<Itr>::value_type sum(Itr b, Itr e)
{
    return std::accumulate(b, e, std::iterator_traits<Itr>::value_type(0));
}

template <typename T>
inline T sum(const T x[], const size_t n)
{
    return std::accumulate(x, x + n, T(0));
}

template <typename T>
inline T mean(const T x[], const size_t n)
{
    T s = sum(x, n);
    return s / T(n);
}

template <typename T>
inline void normalizeToBounds(T x[], const size_t n, const T new_min,
                              const T new_max)
{
    T old_min, old_max;
    old_min = old_max = x[0];
    for (int i = 1; i < n; ++i)
    {
        if (x[i] < old_min)
            old_min = x[i];
        if (x[i] > old_max)
            old_max = x[i];
    }
    const T scale = (new_max - new_min) / (old_max - old_min);
    for (int i = 0; i < n; ++i)
        x[i] = (x[i] - old_min) * scale + new_min;
}

template <typename RandomAccessIterator, typename Predicate>
void insertionSort(RandomAccessIterator begin, RandomAccessIterator end,
                   const Predicate& p)
{
    for (auto i = begin; i != end; ++i)
        std::rotate(std::upper_bound(begin, i, *i, p), i, i + 1);
}

template <typename RandomAccessIterator>
void insertionSort(RandomAccessIterator begin, RandomAccessIterator end)
{
    typedef typename std::iterator_traits<RandomAccessIterator>::value_type val;
    insertionSort(begin, end, std::less<val>());
}

template <typename Itr1, typename Itr2>
inline std::pair<typename std::iterator_traits<Itr1>::value_type,
                 typename std::iterator_traits<Itr1>::value_type> meanStd(Itr1 b, Itr2 e)
{
    using T = typename std::iterator_traits<Itr1>::value_type;
    const std::size_t n = std::distance(b, e);
    T sum = 0, sqr_sum = 0;
    for (; b != e; b++) {
        sum += *b;
        sqr_sum += *b * *b;
    }
    return {sum / n, sqrt((sqr_sum - sum * sum / n) / (n - 1))};
}

template<class It>
class iterable
{
    It m_first, m_last;

public:
    iterable() = default;
    iterable(It first, It last) : m_first(first), m_last(last) {}
    It begin() const { return m_first; }
    It end() const { return m_last; }

};

template<class It>
static inline iterable<It> make_iterable(It a, It b)
{
    return iterable<It>(a, b);
}

// macros
#define Malloc(type,n) (type *)malloc((n) * sizeof(type))
#define StrEq(s1,s2) (!strcmp ((s1),(s2)))

void sortIdx(double x[], int idx[], const int n, const int m);

double sqr(const double& base);
double intpow(double base, int exp);

int binomial(const int n, const int k);
double sqrEuclideanDist(const double a[], const double b[], const int n);
double euclideanDist(const double a[], const double b[], const int n);
int* ranks(const double x[], const int n, int r[]);
double norm(const double x[], const int n);
double* projection(const double v1[], const double v2[], const int n,
                   double ep[]);
double stddev(const double x[], const int n);
double stdWithMean(const double x[], const int n, const double mean);
void colwiseMean(double const* const x[], const int n, const int no_cols,
                 double cm[]);
void colwiseMeanStd(const double* const x[], const int n, const int no_cols,
                    double cm[], double cstd[]);

