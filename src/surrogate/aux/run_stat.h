#pragma once

#include <cmath>
#include <iostream>

class RunningStat {
    size_t m_n;
    double m_oldM, m_newM, m_oldS, m_newS;

public:
    RunningStat() : m_n(0) {}

    void clear() { m_n = 0; }

    void push(const double x) {
        m_n++;

        // See Knuth TAOCP vol 2, 3rd edition, page 232
        if (m_n == 1)
        {
            m_oldM = m_newM = x;
            m_oldS = 0.0;
        }
        else
        {
            m_newM = m_oldM + (x - m_oldM) / m_n;
            m_newS = m_oldS + (x - m_oldM) * (x - m_newM);

            // set up for next iteration
            m_oldM = m_newM;
            m_oldS = m_newS;
        }
    }

    RunningStat& operator<<(const double x) { push(x); return *this; }

    size_t noDataValues() const { return m_n; }
    double mean()         const { return (m_n > 0) ? m_newM : 0.0; }
    double variance() const {
        return ((m_n > 1) ? m_newS / (m_n - 1) : 0.0);
    }
    double standardDeviation() const { return sqrt(variance()); }


    friend std::ostream& operator<<(std::ostream& os, const RunningStat& rs) {
        return os << "no_data:" << rs.noDataValues() << '\n'
                  << "mean: " << rs.mean() << '\n'
                  << "variance: " << rs.variance() << '\n'
                  << "std_dev: " << rs.standardDeviation() << '\n';
    }
};
