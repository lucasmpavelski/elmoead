#include "aux.h"
#include "mo.h"
#include "ml.h"
#include "moead_algs.h"
#include "benchmark.h"
#include "adp.h"
#include "surrogate.h"
#include "estimated_problem.h"

#include "format.h"
#include "INIReader.h"

struct ConfigurableTest {
    using string = std::string;

    long seed;
    string config_fn, front_fn;
    INIReader conf;

    ConfigurableTest(int argc, char* argv[]);

    template <class T>
    void log(const T& v, const string& nm) { fmt::print("{} = {}\n", nm, v); }

    long getl(const string& nm, const long& def)
    {
        auto v = conf.GetInteger("", nm, def);
        log(v, nm);
        return v;
    }
    double getd(const string& nm, const double& def)
    {
        auto v = conf.GetReal("", nm, def);
        log(v, nm);
        return v;
    }
    string gets(const string& nm, const string& def)
    {
        auto v = conf.Get("", nm, def);
        log(v, nm);
        return v;
    }
    bool getb(const string& nm, const bool& def)
    {
        auto v = conf.GetBoolean("", nm, def);
        log(v, nm);
        return v;
    }

    size_t noTrainPartitions(size_t no_vars, size_t no_objs, WeightSet<>::Type wt);

    MOProblem* getProblem();
    OperatorSelection<DEOperator> getDEOperators();
    MOEAD* getMOEAD(size_t no_vars, size_t no_objs);
    std::vector<RegELM> getELMs(const MOProblem& prob);
    ELMOEAD getELMOEAD(MOProblem& prob, MOEAD& moead);
    PopulationCallback getCallback();
};
