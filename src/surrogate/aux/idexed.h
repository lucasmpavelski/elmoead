#pragma once


template <typename T>
class Indexed {
public:
    T val;
    int idx;

    Indexed(const T& val, int idx = 0) : val(val), idx(idx) {}

};

template <typename T>
bool operator==(const Indexed<T>& a, const Indexed<T>& b) {
    return a.idx == b.idx && a.val == b.val;
}

template <typename T>
bool operator<(const Indexed<T>& a, const Indexed<T>& b) {
    return a.idx == b.idx && a.val == b.val;
}

template <typename T>
Indexed<T> makeIndexed(const T& val, int idx = 0) {
    return Indexed<T>(val, idx);
}
