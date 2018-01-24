#pragma once

#include <random>

struct RNG {
    static std::mt19937_64 engine;
    static std::random_device true_rand_engine;

    static long seed(long s = true_rand_engine());

    static long trueRandom() { return true_rand_engine(); }

    template <typename RealT = double>
    inline static RealT realUniform(const RealT& from = 0.0,
                                    const RealT& to = 1.0)
    {
        std::uniform_real_distribution<RealT> uniform_dist(from, to);
        return uniform_dist(engine);
    }

    template <typename IntT = int>
    inline static IntT intUniform(const IntT& from, const IntT& to)
    {
        std::uniform_int_distribution<IntT> uniform_dist(from, to);
        return uniform_dist(engine);
    }

    template <typename IntT = int>
    inline static IntT intUniform(const IntT& until = 100)
    {
        return intUniform(static_cast<IntT>(0), until);
    }

    inline static bool flipCoin()
    {
        static std::bernoulli_distribution dist(0.5);
        return dist(engine);
    }
};
