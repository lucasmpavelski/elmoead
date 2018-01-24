#pragma once

#include <chrono>
#include <utility>

struct Clock {
    using system_clock = std::chrono::system_clock;
    using time_point = std::chrono::time_point<system_clock>;
    using duration = std::chrono::duration<double>;

private:
    static time_point start;

public:
    static void begin();
    static duration tellTime();
};

template<typename TimeT = std::chrono::milliseconds>
struct Measure
{
    template<typename F>
    static typename TimeT::rep execution(F const &func)
    {
        auto start = std::chrono::system_clock::now();
        func();
        auto duration = std::chrono::duration_cast< TimeT>(
            std::chrono::system_clock::now() - start);
        return duration.count();
    }
};

template<typename TimeT = std::chrono::milliseconds>
struct Benchmark
{
    using rep = typename TimeT::rep;

    template<typename F>
    static std::pair<rep,rep> execution(F const &func, const int n)
    {
        const rep sum_m;
        const rep sum_m_sqr;
        for (int i = 0; i < n; ++i)
        {
            const rep m = Measure<TimeT>::execution(func);
            sum_m += m;
            sum_m_sqr += m * m;
        }
        const rep var = (n * sum_m_sqr - (sum_m * sum_m)) / (n * (n - 1));
        return std::make_pair(sum_m / n, sqrt(var));
    }
};
