#pragma once

#include "../aux.h"
#include "../mo.h"

#ifdef __cplusplus
extern "C" {
#endif

void DE_1Bin(const double* in[], const int no_vars, const double cr,
             const double f, double* out);

void DE_2Bin(const double* in[], const int no_vars, const double cr,
             const double f, double* out);

void DE_NonLinear(const double* in[], const int no_vars, double* out);

#ifdef __cplusplus
}
#endif

enum class DEType {
    RAND_1_BIN,
    RAND_2_BIN,
    NON_LINEAR
};

inline DEType str2DEType(const std::string& str)
{
    if (str == "DE/RAND/1/BIN") return DEType::RAND_1_BIN;
    if (str == "DE/RAND/2/BIN") return DEType::RAND_2_BIN;
    if (str == "DE/NONLINEAR" ) return DEType::NON_LINEAR;

    throw_assert(false,
        "de_operator must be (DE/RAND/1/BIN | DE/RAND/2/BIN | DE/NONLINEAR)");

    return DEType::RAND_1_BIN;
}

std::ostream& operator<<(std::ostream& os, const DEType& dt);

class DEOperator {
public:
    DEOperator(const DEType de_type = DEType::RAND_1_BIN,
               const double cr = 1.0,
               const double f = 0.5) : de_type(de_type) , cr(cr) , f(f) {}

    inline DEType type() const { return de_type; }

    inline void apply(const double* in[], const int no_vars,
                      double* out) const
    {
        switch (de_type) {
        case DEType::RAND_1_BIN:
            DE_1Bin(in, no_vars, cr, f, out);
            break;
        case DEType::RAND_2_BIN:
            DE_2Bin(in, no_vars, cr, f, out);
            break;
        case DEType::NON_LINEAR:
            DE_NonLinear(in, no_vars, out);
            break;
        default:
            throw_assert(false, "unknown DE operator");
            break;
        }
    };

    inline void operator()(const double* in[], const int no_vars,
                           double* out) const
    {
        apply(in, no_vars, out);
    };

    friend std::ostream& operator<<(std::ostream& os, const DEOperator& de_op)
    {
        os << "de_operator :\n";
        os << "  de_type: " << de_op.de_type;
        if (de_op.de_type != DEType::NON_LINEAR) {
            os << "\n  cr: " << de_op.cr << "\n";
            os << "  f: " << de_op.f;
        }
        return os;
    };

    friend bool operator==(const DEOperator& de_opa, const DEOperator& de_opb)
    {
        return (de_opa.de_type == de_opb.de_type) && (de_opa.cr == de_opb.cr) && (de_opa.f == de_opb.f);
    };

    static int noParentsFor(const DEType& de_type)
    {
        switch (de_type) {
        case DEType::RAND_1_BIN:
            return 3;
            break;
        case DEType::RAND_2_BIN:
            return 5;
            break;
        case DEType::NON_LINEAR:
            return 3;
            break;
        default:
            break;
        }
    };

    static void selectParents(const DEType& de_type, const Solution pop[],
                              const int candidates_idxs[], const int no_candidates,
                              const int p_idx, double const* parents[])
    {
        const int no_parents = noParentsFor(de_type);
        int parents_idxs[no_parents];
        sampleUniqInt(0, no_candidates, no_parents, parents_idxs);
        parents[0] = pop[p_idx].vars;
        for (int i = 1; i < no_parents; ++i)
            parents[i] = pop[candidates_idxs[parents_idxs[i - 1]]].vars;
    }

private:
    static void sampleUniqInt(const int from, const int to, const int no_samples,
                       int samples[])
    {
        throw_assert(to - from >= no_samples,
              "Invalid range for sampling that number of points");
        for (int i = 0; i < no_samples; i++)
            samples[i] = from - 1;
        int sampled = 0;
        while (sampled < no_samples) {
            const int s = RNG::intUniform(from, to - 1);
            bool found = false;
            for (int i = 0; i < sampled && !found; i++)
                found = (s == samples[i]);
            if (!found) {
                samples[sampled] = s;
                sampled++;
            }
        }
    }

    const DEType de_type;
    const double cr;
    const double f;
};

#include "polynomial_mutation.h"
struct DifferentialEvolution {
    DEOperator de_op;
    PMOperator pmo;

    DifferentialEvolution(DEOperator de_op, PMOperator pmo)
        : de_op(de_op)
        , pmo(pmo)
    {}

    void optimize(size_t max_no_evals, MOProblem& prob, Population& pop,
                  PopulationCallback callback) {

        AllocSolution child(prob);
        const size_t no_parents = de_op.noParentsFor(de_op.type());
        std::vector<const double*> parents(no_parents);
        size_t curr_no_evals = 0;

        while (curr_no_evals < max_no_evals)
        {
            for (int i = 0; i < pop.size(); ++i)
            {
                parents[0] = pop[i].vars;
                size_t im = 1;
                while (im < no_parents)
                {
                    double* selected = pop[RNG::intUniform(size_t(0), pop.size()-1)].vars;
                    if (std::find(parents.begin(), parents.end(), selected) != parents.end())
                        parents[im++] = selected;
                }

                de_op.apply(parents.data(), pop.noVars(), child.vars);
                child.validate();
                pmo.apply(prob.lower_bounds.data(), prob.upper_bounds.data(),
                          prob.noVars(), child.vars, child.vars);
                child.evaluate();
                curr_no_evals++;
 //               cout << curr_no_evals << endl;

                if (child < pop[i])
                    pop[i] = child;

                callback(curr_no_evals, pop);
            }
        }
    }
};
