#pragma once

#include <cmath>
#include <algorithm>
#include <vector>

#include "Eigen/Core"

using std::size_t;

class FuzzyCMeans
{
public:
    using matrix = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,
                                 Eigen::RowMajor>;
    using imatrix = Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic,
                                  Eigen::RowMajor>;
    using vector = Eigen::Matrix<double,1,Eigen::Dynamic>;
    using ivector = Eigen::Matrix<int,1,Eigen::Dynamic>;

    FuzzyCMeans(const size_t no_data, const size_t data_dim,
                const size_t no_clusters, const size_t no_data_per_cluster=-1,
                const double m=2.0, const double eps=0.01) :
        no_data(no_data),
        data_dim(data_dim),
        no_clusters(no_clusters),
        no_data_per_cluster(no_data_per_cluster == -1 ? no_data / no_clusters :
                                                        no_data_per_cluster),
        m(m),
        eps(eps),
        membership(no_data, no_clusters),
        centers(no_clusters, data_dim),
        clusters(no_clusters, no_data_per_cluster),
        idxs_tmp(no_data),
        pertinences_tmp(no_clusters)
    {};

    void resetCenters(const vector& lb, const vector& ub);

    int noClusters() const { return no_clusters; };
    int noDataPerCluster() const { return no_data_per_cluster; };

    double updateMembership(const matrix& data);

    void updateCenters(const matrix& data);

    void cluster(const matrix& data);

    const int clusterData(const int center, const int rank) const {
        return clusters(center, rank);
    };

    void classify(const vector& x, std::vector<int>& ranks);

    int classify(const vector& x);

private:

    double pertinence(const vector& data, const int j) const;;

    static double dist(const vector& a, const vector& b) {
        return (a - b).norm();
    };

    const size_t no_data;
    const size_t data_dim;
    const size_t no_clusters;
    const size_t no_data_per_cluster;
    const double m;
    const double eps;
    matrix membership;
    matrix centers;
    imatrix clusters;

    ivector idxs_tmp;
    vector pertinences_tmp;
};
