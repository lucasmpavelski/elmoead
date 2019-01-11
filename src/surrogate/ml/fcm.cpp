#include "fcm.h"
#include <vector>
#include <numeric>

void FuzzyCMeans::resetCenters(const FuzzyCMeans::vector &lb,
                               const FuzzyCMeans::vector &ub)
{
    centers = ((matrix::Random(no_clusters, data_dim)
                + matrix::Ones(no_clusters, data_dim)) * 0.5)
            .cwiseProduct((ub - lb).replicate(no_clusters, 1))
            + lb.replicate(no_clusters, 1);
}

double FuzzyCMeans::updateMembership(const FuzzyCMeans::matrix& data)
{
    double max_diff = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < no_data; ++i)
    {
        for (int j = 0; j < no_clusters; ++j)
        {
            double u_ij = pertinence(data.row(i), j);
            const double diff = std::abs(u_ij - membership(i, j));
            if (diff > max_diff)
                max_diff = diff;
            membership(i, j) = u_ij;
        }
    }
    return max_diff;
}

void FuzzyCMeans::updateCenters(const FuzzyCMeans::matrix &data)
{
    vector tmp(data_dim);
    for (int j = 0; j < no_clusters; ++j)
    {
        double quotient = 0.0;
        tmp.fill(0.0);
        for (int i = 0; i < no_data; ++i)
        {
            const double u_ij_m = pow(membership(i, j), m);
            quotient += u_ij_m;
            tmp += u_ij_m * data.row(i);
        }
        centers.row(j) = tmp / quotient;
    }
}

void FuzzyCMeans::cluster(const FuzzyCMeans::matrix& data)
{
    resetCenters(data.colwise().minCoeff(), data.colwise().maxCoeff());
    double max_diff = updateMembership(data);
    while (max_diff > eps)
    {
        updateCenters(data);
        max_diff = updateMembership(data);
    }

    for (int i = 0; i < no_clusters; ++i)
    {
        int* idxs = idxs_tmp.data();
        std::iota(idxs, idxs + no_data, 0);

        std::partial_sort(idxs, idxs + no_data_per_cluster, idxs + no_data,
                          [this,&i](const int a, const int b) {
            return this->membership(a, i) > this->membership(b, i);
        });

        for (int j = 0; j < no_data_per_cluster; ++j)
            clusters(i, j) = idxs_tmp(j);
    }
}

void FuzzyCMeans::classify(const FuzzyCMeans::vector& x, std::vector<int>& ranks)
{
    const int no_ranks = ranks.size();
    for (int i = 0; i < no_clusters; ++i)
        pertinences_tmp(i) = pertinence(x, i);
    int* idxs = idxs_tmp.data();
    std::iota(idxs, idxs + no_clusters, 0);
    std::partial_sort(idxs, idxs + no_ranks, idxs + no_clusters,
                      [this](const int a, const int b) {
        return this->pertinences_tmp(a) > this->pertinences_tmp(b);
    });
    for (int i = 0; i < no_ranks; ++i)
        ranks[i] = idxs_tmp(i);
}

int FuzzyCMeans::classify(const FuzzyCMeans::vector &x)
{
    double max_pert = pertinence(x, 0);
    int most_pert = 0;
    for (int i = 1; i < no_clusters; ++i)
    {
        const double p = pertinence(x, i);
        if (p > max_pert)
        {
            max_pert = p;
            most_pert = i;
        }
    }
    return most_pert;
}

double FuzzyCMeans::pertinence(const FuzzyCMeans::vector& data, const int j) const
{
    const double d_xi_cj = dist(data, centers.row(j));
    if (d_xi_cj == 0)
        return 1.0;
    double u = 0.0;
    for (int k = 0; k < no_clusters; ++k)
    {
        const double d_xi_ck = dist(data, centers.row(k));
        if (d_xi_ck == 0.0)
            return 0.0;
        u += pow(d_xi_ck / d_xi_cj, 1.0 / (1.0 - m));
    }
    return 1.0 / u;
}
