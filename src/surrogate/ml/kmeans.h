#pragma once

#include <cmath>
#include <algorithm>

#include "../aux.h"

using std::size_t;

class KMeans
{
public:
    using matrix  = Eigen::MatrixXd;
    using vector  = Eigen::VectorXd;
    using ivector = Eigen::VectorXi;

    matrix centers;
    ivector membership;

    KMeans(size_t dim, size_t no_centers)
        : centers{dim, no_centers}
    {}

    size_t dim()       const { return centers.rows();    }
    size_t noData()    const { return membership.size(); }
    size_t noCenters() const { return centers.cols();    }

    void apply(const matrix& data) {
        throw_assert(data.rows() == dim(), "dimension of data must be the"
                                              " same as defined! (data has " <<
                     data.rows() << " instead of " << dim() << ')');

        membership = ivector::Constant(data.cols(), -1);

        for (int i = 0; i < noCenters(); ++i)
        {
            int rnd;
            do {
                rnd = RNG::intUniform(noData() - 1);
            } while (membership(rnd) != -1);
            centers.col(i) = data.col(rnd);
            membership(rnd) = i;
        }

        bool changed{false};
        vector dists{noCenters()};
        ivector counters{noCenters()};
        matrix new_centers{dim(), noCenters()};

        do {
            changed = false;
            new_centers = matrix::Zero(dim(), noCenters());
            counters = ivector::Zero(noCenters());
            for (int i = 0; i < noData(); ++i)
            {
                int min_idx{-1};
                (centers.colwise() - data.col(i)).colwise().squaredNorm().minCoeff(&min_idx);
                if (membership(i) != min_idx)
                {
                    membership(i) = min_idx;
                    changed = true;
                }
                counters(min_idx) += 1;
                new_centers.col(min_idx) += data.col(i);
            }

            for (int i = 0; i < noCenters(); ++i)
                centers.col(i) = new_centers.col(i) / counters(i);
        } while (changed == true);
    }

    /*void resetCenters(const matrix& data);
    void resetCenters(const vector& lb, const vector& ub);

    double updateMembership(const matrix& data);

    void updateCenters(const matrix& data);

    void cluster(const matrix& data);

    const int clusterData(const int center, const int rank) const {
        return clusters(center, rank);
    };

    void classify(const vector& x, std::vector<int>& ranks);

    int classify(const vector& x);
*/
};
