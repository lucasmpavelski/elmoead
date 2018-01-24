#include "polynomial_mutation.h"

static double sampleY(const double y, const double yl0, const double yu0,
                      const double eta_m)
{
    const double mut_pow = 1.0 / (eta_m + 1.0);
    
    const double yl = (y < yl0) ? y : yl0;
    const double yu = (y > yu0) ? y : yu0;
    const double y_delta = yu - yl;

    double deltaq = 0;
    const double rnd = RNG::realUniform<double>();
    if (rnd <= 0.5)
    {
        const double delta1 = (y - yl) / y_delta;
        const double xy = 1.0 - delta1;
        const double val = 2.0 * rnd + (1.0 - 2.0 * rnd) * pow(xy, eta_m + 1.0);
        deltaq = pow(val, mut_pow) - 1.0;
    }
    else
    {
        const double delta2 = (yu - y) / y_delta;
        const double xy = 1.0 - delta2;
        const double val = 2.0 * (1.0 - rnd) +
                          (2.0 * rnd - 1.0) * pow(xy, eta_m + 1.0);
        deltaq = 1.0 - pow(val, mut_pow);
    }
    
    double y_new = y + deltaq * (yu - yl);
    
    if (y_new < yl)
        y_new = yl;
    if (y_new > yu)
        y_new = yu;
    
    return y_new;
}

void polynomialMutation(const double in[], const double lb[], const double ub[],
                        const int no_vars, const double mut_rate,
                        const double eta_m, double out[])
{
    for (int i = 0; i < no_vars; ++i)
    {
        if (RNG::realUniform<double>() <= mut_rate)
            out[i] = sampleY(in[i], lb[i], ub[i], eta_m);
    }
}

void SBXover(const double lb[], const double ub[], const int no_vars,
             const double xover_rate, const double eta_c,
             const double x1[], const double x2[],
             double offs1[], double offs2[])
{
    static const double eps = 10e-6;
    const double exp_ieta = (1.0 / (eta_c + 1.0));

    if (RNG::realUniform<double>() < xover_rate)
    {
        double rnd;
        double y1, y2, yL, yu;
        double c1, c2;
        double alpha, beta, betaq;
        double valueX1, valueX2;

        for (int i = 0; i < no_vars; i++)
        {
            valueX1 = x1[i];
            valueX2 = x2[i];

            if (RNG::flipCoin())
            {
                if (fabs(valueX1 - valueX2) > eps)
                {
                    if (valueX1 < valueX2)
                    {
                        y1 = valueX1;
                        y2 = valueX2;
                    }
                    else
                    {
                        y1 = valueX2;
                        y2 = valueX1;
                    }

                    yL = lb[i];
                    yu = ub[i];
                    rnd = RNG::realUniform<double>();
                    beta = 1.0 + (2.0 * (y1 - yL) / (y2 - y1));
                    alpha = 2.0 - pow(beta, -(eta_c + 1.0));

                    if (rnd <= (1.0 / alpha))
                        betaq = pow(rnd * alpha, exp_ieta);
                    else
                        betaq = pow(1.0 / (2.0 - rnd * alpha), exp_ieta);

                    c1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
                    beta = 1.0 + (2.0 * (yu - y2) / (y2 - y1));
                    alpha = 2.0 - pow(beta, -(eta_c + 1.0));

                    if (rnd <= (1.0 / alpha))
                        betaq = pow(rnd * alpha, 1.0 / (eta_c + 1.0));
                    else
                        betaq = pow(1.0 / (2.0 - rnd * alpha),
                                    1.0 / (eta_c + 1.0));

                    c2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));

                    if (c1 < yL) c1 = yL;
                    if (c2 < yL) c2 = yL;
                    if (c1 > yu) c1 = yu;
                    if (c2 > yu) c2 = yu;

                    if (RNG::flipCoin())
                    {
                        offs1[i] = c2;
                        offs2[i] = c1;
                    }
                    else
                    {
                        offs1[i] = c1;
                        offs2[i] = c2;
                    }
                }
                else
                {
                    offs1[i] = valueX1;
                    offs2[i] = valueX2;
                }
            }
            else
            {
                offs1[i] = valueX2;
                offs2[i] = valueX1;
            }
        }
    }
    else
    {
        for (int i = 0; i < no_vars; i++)
        {
            offs1[i] = x1[i];
            offs2[i] = x2[i];
        }
    }
}

int tournament(const int competitors[], const int no_comp, const int t_size,
               const std::function<bool(const int&,const int&)>& compare)
{
    int r = RNG::intUniform(no_comp - 1);
    int winner = competitors[r];
    for (int i = 0; i < t_size - 1; ++i)
    {
        r = RNG::intUniform(no_comp - 1);
        if (compare(competitors[r], winner))
            winner = competitors[r];
    }
    return winner;
}
