#ifndef PROBABILITY_MATCHING_H
#define PROBABILITY_MATCHING_H

#include <limits>

#include "../aux.h"
#include "adaptive_operator_selection.h"


template < typename OpT >
class ProbabilityMatching : public OperatorSelection<OpT>
{
public:
    enum class RewardType {
        AvgAbs, AvgNorm, ExtAbs, ExtNorm
    };
    typedef std::vector<double> real_vec;
    typedef std::vector<int> int_vec;

private:

    const double alpha;
    const double p_min;
    RewardType rew_type;
    double best_fitness;
    int chosen_strat;

    real_vec quality;
    real_vec reward;
    real_vec prob;
    real_vec S;
    int_vec nr_S;
    real_vec maior_S;

    void updateRewards();
    double updateQualities();

public:

    using OperatorSelection<OpT>::doAdapt;
    using OperatorSelection<OpT>::noOperators;

    ProbabilityMatching(const std::vector<OpT>& operators,
                        const RewardType& rew_type=RewardType::AvgAbs,
                        const double alpha=0.8, const double p_min=0.1) :
        OperatorSelection<OpT>(operators),
        rew_type(rew_type),
        alpha(alpha),
        p_min(p_min),
        quality(noOperators()),
        reward (noOperators()),
        prob   (noOperators()),
        S      (noOperators()),
        nr_S   (noOperators()),
        maior_S(noOperators())
    {
        warning(noOperators() <= 1,
                "not enough operators for pm adaptation");
        reset(std::numeric_limits<double>::infinity());
    };

    void reset(const double best_f) final override;

    OpT& selectOperator() final override;
    void feedback(const double cf, const double pf) final override;
    void update() final override;

    template< typename OpT2 >
    friend std::ostream& operator<<(std::ostream& os,
                                    ProbabilityMatching<OpT2> const& pm) {
        os << static_cast<const OperatorSelection<OpT>&>(pm)
           << "\n"
           << "  strategy: ProbabilityMatching\n"
           << "  alpha: " << pm.alpha << "\n"
           << "  p_min: " << pm.p_min << "\n"
           << "  reward: " << pm.rew_type;
    }

    friend std::ostream& operator<<(std::ostream& os, RewardType const& rt) {
        switch (rt)
        {
          case RewardType::AvgAbs : os << "AvgAbs" ; break;
          case RewardType::AvgNorm: os << "AvgNorm"; break;
          case RewardType::ExtAbs : os << "ExtAbs" ; break;
          case RewardType::ExtNorm: os << "ExtNorm"; break;
        }
    };
};

template < typename OpT >
void ProbabilityMatching<OpT>::reset(double best_f)
{
    best_fitness = best_f;
    chosen_strat = -1;
    using std::fill;
    fill(quality.begin(), quality.end(), 0.0);
    fill(prob.begin()   , prob.end()   , 1.0 / noOperators());
    fill(S.begin()      , S.end()      , 0.0);
    fill(nr_S.begin()   , nr_S.end()   , 0);
    fill(maior_S.begin(), maior_S.end(), 0.0);
}

template <typename OpT>
OpT& ProbabilityMatching<OpT>::selectOperator()
{
    const double rnd = RNG::realUniform<double>();
    double aux = 0.0;
    for (int k = 0; k < noOperators(); k++)
    {
        aux += prob[k];
        if (rnd <= aux)
        {
            chosen_strat = k;
            return this->getOperator(k);
        }
    }
    throw_assert(false, "probabilities are invalid : ["
                 << printSeq(prob.begin(), prob.end()) << "]"
                 << ", sum = " << sum(prob.begin(), prob.end()));
    return this->getOperator(0);
}

template < typename OpT >
void ProbabilityMatching<OpT>::feedback(const double cf, const double pf)
{
    if (cf < pf)
    {
        const double n = (pf - cf) / pf;
        //((best_fitness + 10e-12) / (cf + 10e-12)) * (pf - cf);
        S[chosen_strat] += n;
        nr_S[chosen_strat]++;
        if (n > maior_S[chosen_strat])
            maior_S[chosen_strat] = n;
    }
};

template < typename OpT >
void ProbabilityMatching<OpT>::update()
{
    this->updateRewards();
    const double q_total = this->updateQualities();

    if (q_total == 0.0) // no improvement, maintain the probabilities
        return;

    const double s = (1.0 - noOperators() * p_min) / q_total;
    for (int k = 0; k < noOperators(); ++k)
        prob[k] = p_min + s * quality[k];

    using std::fill;
    fill(      S.begin(),       S.end(), 0.0);
    fill(   nr_S.begin(),    nr_S.end(), 0  );
    fill(maior_S.begin(), maior_S.end(), 0.0);
}

template< typename OpT >
void ProbabilityMatching<OpT>::updateRewards()
{
    switch (rew_type)
    {
        case RewardType::AvgAbs:
        {
            for (int k = 0; k < noOperators(); ++k)
                reward[k] = nr_S[k] > 0 ? S[k] / nr_S[k] : 0.0;
            break;
        }
        case RewardType::AvgNorm:
        {
            double rew_linha[noOperators()];
            double max_rew_linha = 0.0;
            //calcula o r_linha
            for (int k = 0; k < noOperators(); ++k)
            {
                rew_linha[k] = nr_S[k] > 0 ? S[k] / nr_S[k] : 0.0;
                //acha o maior r_linha
                if (rew_linha[k] > max_rew_linha)
                    max_rew_linha = rew_linha[k];
            }
            //calcula o reward
            for (int k = 0; k < noOperators(); ++k)
                reward[k] = nr_S[k] > 0 ? rew_linha[k] / max_rew_linha : 0.0;
            break;
        }
        case RewardType::ExtAbs:
        {
            for (int k = 0; k < noOperators(); ++k)
                reward[k] = maior_S[k];
            break;
        }
        case RewardType::ExtNorm:
        {
            double rew_linha[noOperators()];
            double max_rew_linha = 0.000001;
            //calcula r_linha
            for (int k = 0; k < noOperators(); k++)
            {
                rew_linha[k] = maior_S[k];
                //acha o maior r_linha
                if (rew_linha[k] > max_rew_linha)
                    max_rew_linha = rew_linha[k];
            }
            //calcula o reward
            for (int k = 0; k < noOperators(); k++)
                reward[k] = rew_linha[k] / max_rew_linha; // WO: division by 0!
            break;
        }
    }
}

template < typename OpT >
double ProbabilityMatching<OpT>::updateQualities()
{
    double q_total = 0.0;
    for (int k = 0; k < noOperators(); ++k)
    {
        quality[k] = quality[k] + alpha * (reward[k] - quality[k]);
        q_total = q_total + quality[k];
    }
    return q_total;
}

#endif // PROBABILITY_MATCHING_H
