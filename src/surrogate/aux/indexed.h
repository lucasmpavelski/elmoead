#pragma once

#include <stddef.h>
#include <iostream>

template <typename T>
struct Indexed {
    using size_t = std::size_t;

    size_t idx;
    T val;

    Indexed(const int& idx = 0, const T& val = 0)
        : idx(idx)
        , val(val)
    {
    }

    static inline bool compareByIdx(const Indexed<T>& a, const Indexed<T>& b)
    {
        return a.idx < b.idx;
    }

    static inline bool compareByVal(const Indexed<T>& a, const Indexed<T>& b)
    {
        return a.val < b.val;
    }

    friend std::ostream& operator<<(std::ostream& o, const Indexed<T>& iv)
    {
        return o << iv.idx << ": " << iv.val;
    }
};
