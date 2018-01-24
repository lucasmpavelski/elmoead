#pragma once

#include <iostream>
#include <vector>

#include "../aux.h"

template< typename OpT >
class OperatorSelection
{
    std::vector<OpT> operators;

public:
    // lifecycle
    OperatorSelection(const std::vector<OpT>& operators) :
        operators(operators)
    {};

    virtual ~OperatorSelection() {};

    // accessors
    int noOperators() const { return operators.size(); };
    bool doAdapt() const { return noOperators() > 1; };
    OpT& getOperator(const int idx) { return operators[idx]; };

    // main interface
    virtual void reset(double d) {}; // on init algorithm
    virtual OpT& selectOperator() { return operators[0]; };
    virtual void feedback(const double cf, const double pf) {};
    virtual void update() {}; // on finish generation

    // misc
    template<typename OpT2>
    friend std::ostream& operator<<(std::ostream& os,
                                    OperatorSelection<OpT2> const& aos) {
        os << "AdaptiveOperatorSelection:\n";
        if (aos.doAdapt())
        {
            os << "  operators:";
            for (const OpT& op : aos.operators)
                os << "\n" << op;
        }
        else
        {
            os << "  constant_operator:\n" << aos.operators[0];
        }
        return os;
    }
};
