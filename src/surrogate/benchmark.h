#pragma once

#include "benchmark/zdt_problems.h"
#include "benchmark/cec09_problems.h"
#include "benchmark/dtlz_problems.h"
#include "benchmark/wfg_problems.h"

inline MOProblem* getBenchmark(const std::string& name, size_t no_vars,
                               size_t no_objs = 2)
{
    if (name == "test" ) return new TestProblem(no_vars, no_objs);
    if (name == "ZDT1" ) return new ZDT1Problem(no_vars);
    if (name == "ZDT2" ) return new ZDT2Problem(no_vars);
    if (name == "ZDT3" ) return new ZDT3Problem(no_vars);
    if (name == "ZDT4" ) return new ZDT4Problem(no_vars);
    if (name == "ZDT6" ) return new ZDT6Problem(no_vars);
    if (name == "UF1"  ) return new UF1Problem(no_vars);
    if (name == "UF2"  ) return new UF2Problem(no_vars);
    if (name == "UF3"  ) return new UF3Problem(no_vars);
    if (name == "UF4"  ) return new UF4Problem(no_vars);
    if (name == "UF5"  ) return new UF5Problem(no_vars);
    if (name == "UF6"  ) return new UF6Problem(no_vars);
    if (name == "UF7"  ) return new UF7Problem(no_vars);
    if (name == "UF8"  ) return new UF8Problem(no_vars);
    if (name == "UF9"  ) return new UF9Problem(no_vars);
    if (name == "UF10" ) return new UF10Problem(no_vars);
    if (name == "DTLZ1") return new DTLZ1Problem(no_vars, no_objs);
    if (name == "DTLZ2") return new DTLZ2Problem(no_vars, no_objs);
    if (name == "DTLZ3") return new DTLZ3Problem(no_vars, no_objs);
    if (name == "DTLZ4") return new DTLZ4Problem(no_vars, no_objs);
    if (name == "DTLZ5") return new DTLZ5Problem(no_vars, no_objs);
    if (name == "DTLZ6") return new DTLZ6Problem(no_vars, no_objs);
    if (name == "DTLZ7") return new DTLZ7Problem(no_vars, no_objs);
    if (name == "WFG1" ) return new WFG1Problem(no_vars, no_objs);
    if (name == "WFG2" ) return new WFG2Problem(no_vars, no_objs);
    if (name == "WFG3" ) return new WFG3Problem(no_vars, no_objs);
    if (name == "WFG4" ) return new WFG4Problem(no_vars, no_objs);
    if (name == "WFG5" ) return new WFG5Problem(no_vars, no_objs);
    if (name == "WFG6" ) return new WFG6Problem(no_vars, no_objs);
    if (name == "WFG7" ) return new WFG7Problem(no_vars, no_objs);
    if (name == "WFG8" ) return new WFG8Problem(no_vars, no_objs);
    if (name == "WFG9" ) return new WFG9Problem(no_vars, no_objs);

    throw_assert(false, "benchmark named " << name << " not found.");

    return nullptr;
}
