#pragma once

#include <cstring>
#include <cstdlib>
#include <algorithm>
using std::fill_n;

#include "../aux.h"

template <typename T>
class AgvNormalizer {
    T smean, sstddev;

public:
    AgvNormalizer()
        : smean(0)
        , sstddev(1)
    {
    }

    void setOldRange(const T l, const T u)
    {
        smean = l;
        sstddev = u;
    }

    template <typename Itr>
    void findOldRange(const Itr& begin, const Itr& end)
    {
        auto r = meanStd(begin, end);
        setOldRange(r.first, r.second);
    }

    T normalize(T v) const noexcept { return (v - smean) / sstddev; }

    T denormalize(T v) const noexcept { return sstddev * v + smean; }

    template <typename Itr1, typename Itr2, typename OutItr>
    void findRangeAndNormalize(Itr1 beg, Itr2 end, OutItr out)
    {
        findOldRange(beg, end);
        normalizeAll(beg, end, out);
    }

    template <typename Itr1, typename Itr2>
    void findRangeAndNormalize(Itr1 beg, Itr2 end)
    {
        findRangeAndNormalize(beg, end, beg);
    }

    template <typename Itr1, typename Itr2>
    void normalizeAll(Itr1 beg, Itr2 end) const { normalizeAll(beg, end, beg); }

    template <typename Itr1, typename Itr2, typename OutItr>
    void normalizeAll(Itr1 beg, Itr2 end, OutItr out) const
    {
        std::transform(beg, end, out, [this](const T& v) -> T {
            return this->normalize(v);
        });
    }

    template <typename Itr1, typename Itr2, typename OutItr>
    void denormalizeAll(Itr1 beg, Itr2 end, OutItr out) const
    {
        std::transform(beg, end, out, [this](const T& v) -> T {
            return this->denormalize(v);
        });
    }

    template <typename Itr1, typename Itr2>
    void denormalizeAll(Itr1 beg, Itr2 end) const
    {
        denormalizeAll(beg, end, beg);
    }

    template <typename T2>
    friend std::ostream& operator<<(std::ostream& os, const AgvNormalizer<T2>& n)
    {
        return os << "AvgNormalizer:\n"
                  << "  mean: " << n.smean << '\n'
                  << "  std: " << n.sstddev << '\n';
    }
};

template <typename T>
class Normalizer {
public:
    Normalizer(const T new_low = 0, const T new_up = 1,
               const T old_low = 0, const T old_up = 1)
        : new_low(new_low)
        , new_up(new_up)
        , old_low(old_low)
        , old_up(old_up)
    {
        update();
        throw_assert(new_low <= new_up, "wrong new range for normalization");
        throw_assert(old_low <= old_up, "wrong old range for normalization");
    }

    void setNewLow(const T v)
    {
        new_low = v;
        throw_assert(new_low <= new_up, "wrong new range for normalization");
        update();
    }

    void setNewUp(const T v)
    {
        new_up = v;
        throw_assert(new_low <= new_up, "wrong new range for normalization");
        update();
    }

    void setNewRange(const T l, const T u)
    {
        new_low = l;
        new_up = u;
        throw_assert(new_low <= new_up, "wrong new range for normalization");
        update();
    }

    void setOldLow(const T v)
    {
        old_low = v;
        throw_assert(old_low <= old_up, "wrong old range for normalization");
        update();
    }

    void setOldUp(const T v)
    {
        old_up = v;
        throw_assert(old_low <= old_up, "wrong old range for normalization");
        update();
    }

    void setOldRange(const T l, const T u)
    {
        old_low = l;
        old_up = u;
        throw_assert(old_low <= old_up, "wrong old range for normalization");
        update();
    }

    T updateOldLow(const T v)
    {
        if (old_low > v)
            setOldLow(v);
        return old_low;
    }

    T updateOldUp(const T v)
    {
        if (old_up < v)
            setOldUp(v);
        return old_up;
    }

    template <typename Itr>
    void findOldRange(const Itr& begin, const Itr& end)
    {
        auto r = std::minmax_element(begin, end);
        setOldRange(*r.first, *r.second);
    }

    T truncateToNewRange(const T v) const noexcept
    {
        if (v < new_low)
            return new_low;
        if (v > new_up)
            return new_up;
        return v;
    }

    T truncateToOldRange(const T v) const noexcept
    {
        if (v < old_low)
            return old_low;
        if (v > old_up)
            return old_up;
        return v;
    }

    T normalize(T v) const noexcept { return range_ratio * v + b; }

    T denormalize(T v) const noexcept { return i_range_ratio * v + i_b; }

    template <typename Itr1, typename Itr2>
    void findRangeAndNormalize(Itr1 beg, Itr2 end)
    {
        findRangeAndNormalize(beg, end, beg);
    }

    template <typename Itr1, typename Itr2, typename OutItr>
    void findRangeAndNormalize(Itr1 beg, Itr2 end, OutItr out)
    {
        findOldRange(beg, end);
        normalizeAll(beg, end, out);
    }

    template <typename Itr1, typename Itr2>
    void normalizeAll(Itr1 beg, Itr2 end) const { normalizeAll(beg, end, beg); }

    template <typename Itr1, typename Itr2, typename OutItr>
    void normalizeAll(Itr1 beg, Itr2 end, OutItr out) const
    {
        std::transform(beg, end, out, [this](const T& v) -> T {
            return this->normalize(v);
        });
    }

    template <typename Itr1, typename Itr2>
    void denormalizeAll(Itr1 beg, Itr2 end) const
    {
        denormalizeAll(beg, end, beg);
    }

    template <typename Itr1, typename Itr2, typename OutItr>
    void denormalizeAll(Itr1 beg, Itr2 end, OutItr out) const
    {
        std::transform(beg, end, out, [this](const T& v) -> T {
            return this->denormalize(v);
        });
    }

    template <typename T2>
    friend std::ostream& operator<<(std::ostream& os, const Normalizer<T2>& n)
    {
        os << "Normalizer:\n"
           << "  old_range: [" << n.old_low << ":" << n.old_up << "]\n"
           << "  new_range: [" << n.new_low << ":" << n.new_up << "]\n";
        return os;
    }

private:
    void update()
    {
        const T old_range = old_up - old_low;
        const T new_range = new_up - new_low;
        range_ratio = new_range / old_range;
        i_range_ratio = old_range / new_range;
        b = new_low - range_ratio * old_low;
        i_b = old_low - i_range_ratio * new_low;
    }

    T old_up, old_low; // original bounds
    T new_up, new_low; // normalized bounds
    T range_ratio, b; // cache for denormalization
    T i_range_ratio, i_b; // cache for denormalization
};

template <typename T>
class NonLinearNormalizer {
public:
    NonLinearNormalizer(const T new_low = 0, const T new_up = 1,
                        const T old_low = 0, const T old_up = 1,
                        const T k = 1.0)
        : new_low(new_low)
        , new_up(new_up)
        , old_low(old_low)
        , old_up(old_up)
    {
        update();
    };

    void setNewLow(const T v)
    {
        new_low = v;
        update();
    };
    void setNewUp(const T v)
    {
        new_up = v;
        update();
    };
    void setNewRange(const T l, const T u)
    {
        new_low = l;
        new_up = u;
        update();
    };

    void setOldLow(const T v)
    {
        new_low = v;
        update();
    };
    void setOldUp(const T v)
    {
        new_up = v;
        update();
    };
    void setOldRange(const T l, const T u)
    {
        old_low = l;
        old_up = u;
        update();
    };

    T updateOldLow(const T v)
    {
        if (old_low > v)
            setOldLow(v);
        return old_low;
    };
    T updateOldUp(const T v)
    {
        if (old_up < v)
            setOldUp(v);
        return old_up;
    };

    template <typename Itr>
    void findOldRange(const Itr& begin, const Itr& end)
    {
        auto r = std::minmax_element(begin, end);
        setOldRange(*r.first, *r.second);
    };

    T truncateToNewRange(const T v) const noexcept
    {
        if (v < new_low)
            return new_low;
        if (v > new_up)
            return new_up;
        return v;
    };

    T truncateToOldRange(const T v) const noexcept
    {
        if (v < old_low)
            return old_low;
        if (v > old_up)
            return old_up;
        return v;
    };

    T truncate(const T v) const noexcept
    {
        if (v < 0.0)
            return 0.0;
        if (v > 1.0)
            return 1.0;
        return v;
    };

    T normalize(T x) const noexcept
    {
        const T norm_x = x * a + b; //std::pow(truncate(x * a + b), 3.0);
        x = std::log2(pow(truncate(norm_x), 1.0 / k) + 1) * new_range + new_low;
        return x;
    };

    T denormalize(T fx) const noexcept
    {
        const T norm_fx = (fx * ia + ib);
        fx = pow(pow(2.0, truncate(norm_fx)) - 1, k) * old_range + old_low;
        return fx;
    };

    template <typename Itr1, typename Itr2>
    void findRangeAndNormalize(Itr1 beg, Itr2 end)
    {
        findRangeAndNormalize(beg, end, beg);
    };

    template <typename Itr1, typename Itr2, typename OutItr>
    void findRangeAndNormalize(Itr1 beg, Itr2 end, OutItr out)
    {
        findOldRange(beg, end);
        normalizeAll(beg, end, out);
    };

    template <typename Itr1, typename Itr2>
    void normalizeAll(Itr1 beg, Itr2 end) { normalizeAll(beg, end, beg); };

    template <typename Itr1, typename Itr2, typename OutItr>
    void normalizeAll(Itr1 beg, Itr2 end, OutItr out)
    {
        std::transform(beg, end, out, [this](const T& v) -> T {
            return this->normalize(v);
        });
    };

    template <typename Itr1, typename Itr2>
    void denormalizeAll(Itr1 beg, Itr2 end) { denormalizeAll(beg, end, beg); };

    template <typename Itr1, typename Itr2, typename OutItr>
    void denormalizeAll(Itr1 beg, Itr2 end, OutItr out)
    {
        std::transform(beg, end, out, [this](const T& v) -> T {
            return this->denormalize(v); });
    };

private:
    void update()
    {
        old_range = old_up - old_low;
        new_range = new_up - new_low;
        a = 1.0 / old_range;
        ia = 1.0 / new_range;
        b = -old_low * a;
        ib = -new_low * ia;
    };

    T k;
    T old_up, old_low; // original bounds
    T new_up, new_low; // normalized bounds
    T old_range, new_range;
    T a, b; // cache for normalization
    T ia, ib; // cache for denormalization
};
