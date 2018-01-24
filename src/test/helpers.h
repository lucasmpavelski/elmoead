#ifndef HELPERS_H
#define HELPERS_H

#include <cstring>
#include <cctype>

#include <string>
using std::stringstream;

#include "gtest/gtest.h"
using ::testing::AssertionResult;
using ::testing::AssertionSuccess;
using ::testing::AssertionFailure;
using ::testing::internal::Double;

#include <chrono>
#include <utility>
/*
template <typename Itr1, typename Itr2>
AssertionResult allEq(Itr1 fst, Itr1 lst, Itr2 ofst)
{
    auto r = std::mismatch(fst, lst, ofst);
    if (r.first != lst)
    {
        const size_t i = std::distance(fst, r.first);
        return AssertionFailure() << i + 1
                                  << ((i == 0) ? "st" :
                                      (i == 1) ? "nd" :
                                      (i == 2) ? "rd" : "th")
                                  << " elements differ ("
                                  << r.first << " and " << r.second;
    }
    return AssertionSuccess();
}
*/
template <typename T>
AssertionResult allEq(const T* a, const T* b, const int n)
{
    //return allEq(a, a + n, b);

    if (a == b) return AssertionSuccess();
    for (int i = 0; i < n; i++)
    {
        if (a[i] != b[i])
        {
            return AssertionFailure() << i + 1
                                      << ((i == 0) ? "st" :
                                          (i == 1) ? "nd" :
                                          (i == 2) ? "rd" : "th")
                                      << " elements differ ("
                                      << a[i] << " and " << b[i];
        }
    }
    return AssertionSuccess();
}


template <typename T>
AssertionResult allDoubleEq(const T* a, const T* b, const int n)
{
    if (a == b) return AssertionSuccess();
    for (int i = 0; i < n; i++)
    {
        Double a_i(a[i]);
        Double b_i(b[i]);
        if (!a_i.AlmostEquals(b_i))
        {
            return AssertionFailure() << i + 1
                                      << ((i == 0) ? "st" :
                                          (i == 1) ? "nd" :
                                          (i == 2) ? "rd" : "th")
                                      << " elements differ ("
                                      << a[i] << " and " << b[i];
        }
    }
    return AssertionSuccess();
}

template <typename T>
AssertionResult contain(const T* set, const T& el, const int n)
{
    if (n <= 0) return AssertionFailure() << "the set is empty";
    bool contain = false;
    for (int j = 0; (j < n) && (!contain); ++j) {
        if (set[j] == el)
            contain = true;
    }
    if (!contain)
    {
        stringstream ss;
        ss << set[0];
        for (int j = 1; j < n; ++j)
            ss << set[j] << ",";
        return AssertionFailure() << el << " not in {" << ss.str() << "}";
    }
    return AssertionSuccess();
}

template <typename T>
AssertionResult containAll(const T* set1, const int n1, const T* set2, const int n2)
{
    for (int i = 0; i < n2; ++i)
    {
        const T& el = set2[i];
        AssertionResult r = contain(set1, el, n1);
        if (!r)
        {
            stringstream ss;
            ss << set1[0];
            for (int j = 1; j < n1; ++j)
                ss << "," << set1[j];
            return AssertionFailure() << el << " not in {" << ss.str() << "}";
        }

    }
    return AssertionSuccess();
}

AssertionResult isSubString(const char* haystack, const char* needle)
{
    if (strstr(haystack, needle) != NULL)
        return AssertionSuccess();
    return AssertionFailure() << "\"" << needle << "\" not found in \""
                              << haystack << "\"";
}

template <typename T>
AssertionResult isWithin(const T& x, const T& a, const T& b)
{
    if ((x >= a) && (x <= b))
        return AssertionSuccess();
    return AssertionFailure() << x << " isnt within [" << a << "," << b << "]";
}

#endif // HELPERS_H
