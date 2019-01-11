#include "adp.h"
#include "format.h"

#include <vector>

void AirfoilDesignProblem::evaluate(const double x[], double y[])
{
    //cerr << printSeq(x, x + noVars()) << endl;
    //throw_assert(isValid(x), "invalid x: " << printSeq(x, x + noVars()));
    auto beg = Clock::tellTime();

    auto vpar = buildCmd(x, 2040000, 0.12, 0.63);
    log_file << vpar << "\n";
    auto out = exec2Str(vpar.c_str());
    log_file << out;

    double cd_1 = std::stod(out);
    y[0] = cd_1 / 0.63;

    vpar = buildCmd(x, 1290000, 0.08, 1.05);
    log_file << vpar << "\n";
    out = exec2Str(vpar.c_str());
    log_file << out;

    double cd_2 = std::stod(out);
    y[1] = cd_2 / pow(1.05, 3./2.);

    auto td = Clock::tellTime() - beg;

    log_file << "time: " << td.count() << "\n" << endl;
}

std::string AirfoilDesignProblem::buildCmd(const double x[], const double re,
                                           const double mach, const double cl)
{
    return fmt::sprintf("adp/eval.py "
            "--r_leup"   "=%.10f "
            "--r_lelo"   "=%.10f "
            "--alpha_te" "=%.10f "
            "--beta_te"  "=%.10f "
            "--z_te"     "=%.10f "
            "--dz_te"    "=%.10f "
            "--x_up"     "=%.10f "
            "--z_up"     "=%.10f "
            "--z_xxup"   "=%.10f "
            "--x_lo"     "=%.10f "
            "--z_lo"     "=%.10f "
            "--z_xxlo"   "=%.10f "
            "--re"       "=%.10f "
            "--mach"     "=%.10f "
            "--cl"       "=%.10f",
            x[0], x[1], x[2], x[3], x[4], x[5], x[6],
            x[7], x[8], x[9], x[10], x[11], re, mach, cl);
 }
