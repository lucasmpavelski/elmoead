#pragma once

#include <iostream>
#include <algorithm>
#include <vector>

#include "../aux.h"

template <typename T = double>
class WeightSet {
public:
    using matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using vector = Eigen::Matrix<T, 1, Eigen::Dynamic>;
    enum class Type { UNIFORM,
                      RANDOM };

private:
    size_t dim;
    size_t no_partitions;
    size_t no_weights;
    Type   mtype;
    double eps;

public:
    matrix weights;

    WeightSet(const size_t dim, const size_t no_partitions,
              const Type type = Type::UNIFORM,
              const double eps = 1e-6)
        : dim(dim)
        , no_partitions(no_partitions)
        , no_weights((type == Type::UNIFORM) ? WeightSet::sizeFor(dim, no_partitions) : no_partitions)
        , weights((type == Type::UNIFORM) ? uniformWeights(dim, no_partitions, eps) : randomWeights(dim, no_partitions))
        , mtype(type)
        , eps(eps)
    {}

    WeightSet(const T* values, const int no_weights, const int no_objs)
        : dim(no_objs)
        , no_partitions(-1)
        , no_weights(no_weights)
        , weights(Eigen::Map<const matrix>(values, no_weights, no_objs))
        , mtype(Type::RANDOM)
        , eps(0)
    {}

    void swap(WeightSet& o) noexcept
    {
        using std::swap;
        swap(o.dim, this->dim);
        swap(o.no_partitions, this->no_partitions);
        swap(o.no_weights, this->no_weights);
        swap(o.weights, this->weights);
        swap(o.eps, this->eps);
    }

    size_t size() const { return no_weights; }
    size_t dimension() const { return dim; }
    Type   type() const { return mtype; }

    inline const T* operator[](const int i) const { return weights.row(i).data(); }

    void invert() { weights.cwiseInverse(); }

    template <class ValT>
    class iterator {
    public:
        typedef iterator self_type;
        typedef ValT& reference;
        typedef ValT value_type;
        typedef ValT* ptr_type;
        typedef std::forward_iterator_tag iterator_category;
        typedef int difference_type;
        iterator(ptr_type begin, difference_type curr_pos, difference_type inc)
            : curr(begin + curr_pos * inc)
            , inc(inc){}
        self_type operator++()
        {
            curr += inc;
            return *this;
        }
        self_type operator++(int dummy) { return operator++(); }
        bool operator==(const self_type& rhs) { return curr == rhs.curr; }
        bool operator!=(const self_type& rhs) { return curr != rhs.curr; }
        ptr_type operator->() { return curr; }
        ptr_type operator*() { return curr; }

    private:
        ptr_type curr;
        difference_type inc;
    };

    iterator<T> begin() { return iterator<T>(weights.data(), 0, dim); }
    iterator<T> end() { return iterator<T>(weights.data(), size(), dim); }

    iterator<const T> begin() const { return iterator<const T>(weights.data(), 0, dim); }
    iterator<const T> end() const { return iterator<const T>(weights.data(), size(), dim); }

    static matrix uniformWeights(const int m, const int p, const double eps)
    {
        matrix w(WeightSet::sizeFor(m, p), m);
        vector acc_dims(m);
        acc_dims.fill(0.0);
        act_point = 0;
        uniformRecursive(m, p, p + 1, w, acc_dims, 0, eps);
        return w;
    }

    static matrix randomWeights(const int dim, const int no_weights)
    {
        matrix w(no_weights, dim);
        for (int i = 0; i < no_weights; ++i) {
            double s = 0.0;
            for (int j = 0; j < dim - 1; ++j) {
                const double r = RNG::realUniform(0.0, 1.0);
                w(i, j) = (1.0 - s) * (1.0 - std::pow(r, 1.0 / (dim - j)));
                s += w(i, j);
            }
            w(i, dim - 1) = 1.0 - s;
        }
        return w;
    }


    /*    inline T eval(AggregationFunctionType af, const int idx, const T* objsp,
                  const T* idealp, const T* nadirp) const {
        const Eigen::Map<const vector> objs(objsp, no_objs);
        return eval(af, idx,
                    Eigen::Map<const vector>(objs, no_objs),
                    Eigen::Map<const vector>(ideal, no_objs),
                    Eigen::Map<const vector>(nadir, no_objs));
    };

    inline T eval(AggregationFunctionType af, const int idx,
                  const Eigen::Ref<const vector>& objs,
                  const Eigen::Ref<const vector>& ideal,
                  const Eigen::Ref<const vector>& nadir) const {
        const auto& w = weights.row(idx);
        switch(af)
        {
        case WEIGHTED_SUM: {
            return w.transpose().dot(objs);
        }
        break;
        case TCHEBYCHEFF: {
            const Eigen::Map<const vector> ideal(idealp, no_objs);
            return (w.cwiseProduct(objs - ideal)).maxCoeff();
        }
        break;
        case INVERTED_TCHEBYCHEFF:{
            const Eigen::Map<const vector> ideal(idealp, no_objs);
            return (w.cwiseInverse().cwiseProduct(objs - ideal)).maxCoeff();
        }
        break;
        case AUGMENTED_TCHEBYCHEFF:{
            const Eigen::Map<const vector> ideal(idealp, no_objs);
            return (w.cwiseProduct(objs - ideal)).maxCoeff() +
                    0.01 * weights.row(idx).transpose().dot(objs);
        }
        break;
        case NORMALIZED_TCHEBYCHEFF:{
            const Eigen::Map<const vector> ideal(idealp, no_objs);
            const Eigen::Map<const vector> nadir(nadirp, no_objs);
            return (w.cwiseProduct(objs - ideal)
                    .cwiseQuotient(nadir - ideal)).maxCoeff();
        }
        break;
        case PENALTY_BOUNDARY_INTERSECTION: {
            const Eigen::Map<const vector> ideal(idealp, no_objs);
            const double w_norm = weights.row(idx).norm();
            const T d1 = weights.row(idx).cwiseProduct(objs - ideal).sum() / w_norm;
            const T d2 = (objs - ideal - d1 * weights.row(idx)).norm();
            return d1 + PBI_THETA * d2;
        }
        break;
        default:
            error(true, "unknown aggregation function type");
        break;
        }
    };
*/
    static size_t sizeFor(const int no_objs, const int no_part)
    {
        return binomial(no_objs + no_part - 1, no_part);
    }

    template <typename Tp>
    friend std::ostream& operator<<(std::ostream& os, const WeightSet<Tp>& w)
    {
        os << "weights(" << w.size() << "x" << w.dimension() << ")\n"
           << w.weights << "\n";
        return os;
    }

private:
    static int act_point;

    static void uniformRecursive(const int m, const int p, const int end,
                                 matrix& points, vector acc_dims, const int act_dim, const double meps)
    {
        if (act_dim == m - 1) {
            acc_dims[m - 1] = static_cast<T>(end - 1) / p;
            for (int i = 0; i < m; ++i)
                points(act_point, i) = (acc_dims[i] == 0) ? meps : acc_dims[i];
            act_point++;
        } else {
            for (int i = 0; i < end; ++i) {
                acc_dims[act_dim] = static_cast<T>(i) / p;
                uniformRecursive(m, p, end - i, points, acc_dims, act_dim + 1,
                                 meps);
            }
        }
    }
};

template <typename T>
int WeightSet<T>::act_point = 0;
