#pragma once

#include <string>
#include <cstdio>
#include <iomanip>

#include "mo.h"
#include "aux.h"

class AirfoilDesignProblem : public MOProblem
{
public:
    AirfoilDesignProblem(const std::string& log_fn = "adp.log")
        : MOProblem("adp", 12, 2,
//r_leup|r_lelo|alpha_te|beta_te| z_te|dz_te|x_up|z_up|z_xxup|x_lo|z_lo|z_xxlo
 { .0085,  .002,       7,     10,-.006,.0025, .41, .11,   -.9, .20,-.023,  .05},
 { .0126,  .004,      10,     14,-.003,.0050, .46, .13,   -.7, .26,-.015,  .20})
        , log_file(log_fn.c_str())
    {}

    // MOProblem interface
public:
    void evaluate(const double x[], double y[]) override final;

    void validate(double x[]) override final { truncateToBounds(x); }

private:

    std::ofstream log_file;

    std::string buildCmd(const double x[], const double re, const double mach,
                         const double cl);
};
