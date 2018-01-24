#pragma once

#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <vector>

#include "../aux.h"
#include "weights.h"

struct Neighborhood
{
    typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> imatrix;
    typedef Eigen::Matrix<int,Eigen::Dynamic,1> ivector;
    typedef Eigen::Matrix<double,Eigen::Dynamic,1> dvector;

    size_t no_weights_in;
    size_t no_neighbors;
    imatrix neighbors;
//    ivector inverse_map;

    template <typename T>
    Neighborhood(const WeightSet<T>& in, const WeightSet<T>& out,
                 const size_t t) :
        no_weights_in(in.size()),
        no_neighbors(t),
        neighbors(neighborhoodFor(in, out, t))//,
//        inverse_map(buildInverse(neighbors))
    {};

    template <typename T>
    Neighborhood(const WeightSet<T>& in, const WeightSet<T>& out) :
        Neighborhood{in, out, out.size() / in.size()}
    {};

    template <typename T>
    Neighborhood(const WeightSet<T>& inout, const size_t t) :
        Neighborhood{inout, inout, t}
    {};

    Neighborhood(int* data, const size_t no_weights_in, const size_t t) :
        no_weights_in(no_weights_in),
        no_neighbors(t),
        neighbors(Eigen::Map<imatrix>(data, no_weights_in, t))
    {};

    inline const int* operator[](const size_t i) const {
        return neighbors.row(i).data();
    };

    template <typename T>
    static imatrix neighborhoodFor(const WeightSet<T>& win,
                                   const WeightSet<T>& wout, const size_t t) {
        const size_t no_ins = win.size();
        const size_t no_outs = wout.size();

        assert(t <= no_outs);

        imatrix B(win.size(), t);
        dvector x(no_outs);
        std::vector<int> idx(no_outs);

        for (int i = 0; i < no_ins; ++i)
        {
            for (int j = 0; j < no_outs; ++j)
                x[j] = (win.weights.row(i) - wout.weights.row(j)).squaredNorm();

            std::iota(idx.begin(), idx.end(), 0);
            std::partial_sort(idx.begin(), idx.begin() + t, idx.end(),
                              [&x](const int& a, const int& b) {
                return x[a] < x[b];
            });

            std::copy(idx.begin(), idx.begin() + t, B.row(i).data());
        }

        return B;
    };

    /*ivector buildInverse(const imatrix& B) {
        ivector r(B.size());
        for (int i = 0; i < B.rows(); ++i)
            for (int j = 0; j < B.cols(); ++j)
                r[B(i, j)] = i;
        return r;
    };*/

    friend std::ostream& operator<<(std::ostream& os, const Neighborhood& n) {
        return os << "neighborhood(" << n.no_weights_in << "x" << n.no_neighbors
                  << ")\n" << n.neighbors << "\n";
    };
};

