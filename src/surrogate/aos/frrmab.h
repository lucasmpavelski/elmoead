#ifndef FRRMAB_H
#define FRRMAB_H

#include <vector>
#include <algorithm>

#include "../aux.h"
#include "adaptive_operator_selection.h"

template<typename OpT>
class FRRMAB : public OperatorSelection<OpT>
{
public:
    using OperatorSelection<OpT>::doAdapt;
    using OperatorSelection<OpT>::noOperators;

    FRRMAB(const std::vector<OpT>& strategies, const int window_size,
           //const double scale = 5.0, const double decay = 1.0 //1st conf
           //const double scale = 1.0, const double decay = 0.5 //2nd conf
            const double scale = 5.0, const double decay = 0.5 //3rd conf
          //  const double scale = 1.0, const double decay = 1.0 //4th conf
            ) :
        OperatorSelection<OpT>(strategies),
        scale(scale),
        decay(decay),
        last_op(-1),
        fir_records(window_size, Indexed<double>(-1, 0.0)),
        fir         (noOperators()),
        frr         (noOperators()),
        reward      (noOperators()),
        num         (noOperators()),
        ranks       (noOperators()),
        not_selected(noOperators()),
        unused_operators_exist(true)
    {
        reset(0.0);
    };

    void reset(double d) final override;

    void update() final override {
        std::fill(fir.begin(), fir.end(), 0.0);
    };
    OpT& selectOperator() final override;
    void feedback(const double cf, const double pf) final override {
        if (cf < pf)
            fir[last_op] += (pf - cf) / pf;
    };

    template< typename OpT2 >
    friend std::ostream& operator<<(std::ostream& os, FRRMAB<OpT2> const& fm);

private:

    void assignCredits();

    const double scale;
    const double decay;
    int last_op;

    SlidingWindow<Indexed<double>> fir_records;

    typedef std::vector<double> real_vec;
    typedef std::vector<int> int_vec;
    typedef std::vector<bool> bool_vec;

    real_vec fir;
    real_vec frr;
    real_vec reward;
    int_vec num;
    int_vec ranks;

    bool_vec not_selected;
    bool unused_operators_exist;
};

template< typename OpT >
void FRRMAB<OpT>::reset(double d)
{
    last_op = -1;
    fir_records.clear(Indexed<double>(-1, 0.0));
    using std::fill;
    fill(fir.begin()     , fir.end()     , 0.0);
    fill(frr.begin()     , frr.end()     , 0.0);
    fill(reward.begin()  , reward.end()  , 0.0);
    fill(num.begin()     , num.end()     , 0);
    fill(not_selected.begin(), not_selected.end(), true);
    unused_operators_exist = true;
};

template<typename OpT>
OpT& FRRMAB<OpT>::selectOperator()
{
    if (last_op != -1)
        assignCredits();

    int idx = -1;

    if (unused_operators_exist)
    {
        idx = RNG::intUniform(0, noOperators() - 1);
        not_selected[idx] = false;
        unused_operators_exist = std::any_of(not_selected.begin(),
                                             not_selected.end(),
                                             [](const bool b){ return b; });
    }
    else
    {
        double max_opt = -numeric_limits<double>::infinity();
        const double log_sum = 2.0 * log(sum(num.data(), num.size()));
        for (int i = 0; i < noOperators(); ++i)
        {
            if (num[i] > 0)
            {
                const double opt = frr[i] + scale * sqrt(log_sum / num[i]);
                if (opt > max_opt)
                {
                    max_opt = opt;
                    idx = i;
                }
            }
        }
        throw_assert(idx != -1, "invalid operators usage counters in MAB : ["
                     << printSeq(num.begin(), num.end()) << "]");
    }
    last_op = idx;
    return this->getOperator(idx);
}

template< typename OpT >
void FRRMAB<OpT>::assignCredits()
{
    std::fill(reward.begin(), reward.end(), 0.0);
    std::fill(num.begin(), num.end(), 0);

    fir_records.append(Indexed<double>(last_op, fir[last_op]));
    for (const auto& el : fir_records)
    {
        if (el.idx >= 0)
        {
            reward[el.idx] += el.val;
            num[el.idx]++;
        }
    }

    std::iota(ranks.begin(), ranks.end(), 0);
    std::sort(ranks.begin(), ranks.end(), [this](const int& a, const int& b) {
        return this->reward[a] > this->reward[b];
    });

    last_op = -1;
    double sum = 0.0;
    for (int i = 0; i < noOperators(); ++i)
    {
        frr[i] = intpow(decay, ranks[i]) * reward[i];
        sum += frr[i];
    }
    if (sum > 0.0)
    {
        for (int i = 0; i < noOperators(); ++i)
            frr[i] = frr[i] / sum;
    }
}

template< typename OpT >
std::ostream& operator<<(std::ostream& os, FRRMAB<OpT> const& fm)
{
    os << static_cast<OperatorSelection<OpT> const&>(fm)
       << "\n"
       << "  scale: " << fm.scale << "\n"
       << "  decay: " << fm.decay << "\n"
       << "  sliding_window_size: " << fm.fir_records.size();
    return os;
};


#endif // FRRMAB
