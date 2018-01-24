#include "de_operators.h"

/*template<>
void DEOperator2<DEType::RAND_1_BIN>::apply(const double* in[], const int no_vars,
                                          double* out) const
{
    DE_1Bin(in, no_vars, cr, f, out);
};

template<>
void DEOperator2<DEType::RAND_2_BIN>::apply(const double* in[], const int no_vars,
                                          double* out) const
{
    DE_2Bin(in, no_vars, cr, f, out);
};

template<>
void DEOperator2<DEType::NON_LINEAR>::apply(const double* in[], const int no_vars,
                                          double* out) const
{
    DE_NonLinear(in, no_vars, out);
};
template<>
const int DEOperator2<DEType::RAND_1_BIN>::no_parents = 3;
template<>
const int DEOperator2<DEType::RAND_2_BIN>::no_parents = 5;
template<>
const int DEOperator2<DEType::NON_LINEAR>::no_parents = 3;
*/

std::ostream& operator<<(std::ostream& os, const DEType& dt)
{
    switch (dt)
    {
        case DEType::RAND_1_BIN: os << "RAND_1_BIN"; break;
        case DEType::RAND_2_BIN: os << "RAND_2_BIN"; break;
        case DEType::NON_LINEAR: os << "NON_LINEAR"; break;
        default: throw_assert(false, "unknown DE operator"); break;
    }
}

void DE_1Bin(const double* in[], const int no_vars, const double cr, 
        const double f, double* out)
{
    const double* r0 = in[0];
    const double* r1 = in[1];
    const double* r2 = in[2];

    //Aplicar diferenciação, mutação, crossover
    const int j_rand = RNG::intUniform(no_vars - 1); // para garantir q pelo menos um parâmetro será modificado
    for (int j = 0; j < no_vars; ++j)
    {
        //testa se faz a mutação
        if ((RNG::realUniform<double>() <= cr) || (j == j_rand))
            out[j] = r0[j] + f * (r1[j] - r2[j]);
        else
            out[j] = r0[j];
    }
}

void DE_2Bin(const double* in[], const int no_vars, const double cr, 
        const double f, double* out)
{
    const double* r0 = in[0];
    const double* r1 = in[1];
    const double* r2 = in[2];
    const double* r3 = in[3];
    const double* r4 = in[4];

    //Aplicar diferenciação, mutação, crossover
    const int j_rand = RNG::intUniform(no_vars - 1); // para garantir q pelo menos um parâmetro será modificado
    for (int j = 0; j < no_vars; ++j)
    {
        //testa se faz a mutação
        if ((RNG::realUniform<double>() <= cr) || (j == j_rand))
            out[j] = r0[j] + f * ((r1[j] - r2[j]) + (r3[j] - r4[j]));
        else
            out[j] = r0[j];
    }
}

void DE_NonLinear(const double* in[], const int no_vars, double* out)
{
    // DE_nao_linear não usa CR nem F
    static const double P_inter = 0.75;
    
    const double* r0 = in[0];
    const double* r1 = in[1];
    const double* r2 = in[2];

    for (int j = 0; j < no_vars; ++j)
    {
        const double c0 = r0[j];
        const double c1 = (4.0 * r1[j] - 3.0 * r0[j] - r2[j]) * 0.5;
        const double c2 = (r0[j] - 2.0 * r1[j] + r2[j]) * 0.5;
        
        const double t = (RNG::realUniform<double>() <= P_inter) ?
                    RNG::realUniform(0.0, 2.0) : RNG::realUniform(2.0, 3.0);
        
        out[j] = c2 * t * t + c1 * t + c0;
    }
}
