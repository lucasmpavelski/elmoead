#pragma once

#include <vector>
#include <utility>
#include <algorithm>

template <class T>
using bound = std::pair<T, T>;

template <class T>
using bounds = std::vector<bound<T>>;

template <class Iter1, class Iter2>
bounds<typename std::iterator_traits<Iter1>::value_type>
makeBounds(Iter1 lb, Iter1 le, Iter2 ub, Iter2 ue)
{
    using T = typename std::iterator_traits<Iter1>::value_type;
    bounds<T> r;
    r.reserve(std::min(std::distance(lb, le), std::distance(ub, ue)));
    while (lb != le && ub != ue)
    {
        r.emplace_back(*lb, *ub);
        ++lb;
        ++ub;
    }
    return r;
}
/*
template <class Seq>
bounds<typename std::remove_reference<Seq>::type>
makeBounds(Seq&& low, Seq&& up)
{
    using std::begin;
    using std::end;
    return makeBounds(begin(low), end(low), begin(up), end(up));
}
*/
template <class T>
bounds<T> makeBounds(std::initializer_list<T> low, std::initializer_list<T> up)
{
    using std::begin;
    using std::end;
    return makeBounds(begin(low), end(low), begin(up), end(up));
}

template <class T>
bool isBetween(const bound<T>& bd, const T& el)
{
    return (el >= bd.first) && (el <= bd.second);
}

template <class T, class Iter>
bool isBetween(const bounds<T>& bds, Iter b, Iter e)
{
    using std::begin;
    using std::end;
    return std::equal(begin(bds), end(bds), b, isBetween<T>);
}

template <class T, class Seq>
bool isBetween(const bounds<T>& bds, Seq&& seq)
{
    using std::begin;
    using std::end;
    return isBetween(bds, begin(seq), end(seq));
}

template <class T1, class T2>
bool isBetween(const bounds<T1>& bds, std::initializer_list<T2>&& seq)
{
    using std::begin;
    using std::end;
    return isBetween(bds, begin(seq), end(seq));
}

template <class T>
T truncate(const bound<T>& bd, const T& el)
{
    if (el < bd.first)
        return bd.first;
    if (el > bd.second)
        return bd.second;
    return el;
}

template <class T, class Iter>
void truncate(const bounds<T>& bds, Iter b)
{
    using std::begin;
    using std::end;
    std::transform(begin(bds), end(bds), b, b, truncate<T>);
}

template <class T>
T reflect(const bound<T>& bd, const T& x)
{
    while (x < bd.first)
      x = bd.second - std::fabs(x - bd.first);
    while (x > bd.second)
      x = bd.first + std::fabs(x - bd.second);
    return x;
}

template <class T, class Iter>
void reflect(const bounds<T>& bds, Iter b)
{
    using std::begin;
    using std::end;
    std::transform(begin(bds), end(bds), b, b, reflect<T>);
}

