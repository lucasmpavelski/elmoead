#include "zdt_problems.h"
/*
void EvalKNO1(void* prob, const double x[], double y[])
{
    double c = x[0] + x[1];
    double f = 20 - (11 + 3 * sin((5 * c) * (0.5 * c)) +
                          3 * sin(4 * c) +
                          5 * sin(2 * c + 2));
    double g = (M_PI / 2.0) * (x[0] - x[1] + 3.0) / 6.0;

    y[0] = 20 - (f * cos(g));
    y[1] = 20 - (f * sin(g));
}

void initKNO1Problem(MOProblem* prob)
{
    prob->name = (char*) malloc(sizeof(char) * 5);
    strcpy(prob->name, "KNO1");
    prob->no_vars = 2;
    prob->no_objs = 2;
    double* lb = new double[2];
    lb[0] = lb[1] = 0.0;
    double* ub = new double[2];
    printf("%p\n", ub);
    ub[0] = 3.0;
    ub[1] = 3.0;

    prob->lower_bounds = lb;
    prob->upper_bounds = ub;
    prob->info = NULL;
    prob->evaluate = EvalKNO1;
    prob->validate = validateProbTruncating;
}

void deinitKNO1Problem(MOProblem* prob)
{
    free(prob->name);
    free(prob->lower_bounds);
    free(prob->upper_bounds);
}

static void initAnyZDT(MOProblem* prob, const char* name, const int no_vars,
                       evaluationFunc fe, const double lb, const double ub)
{
    prob->name = new char[strlen(name) + 1];
    strcpy(prob->name, name);
    prob->no_vars = no_vars;
    prob->no_objs = 2;
    prob->lower_bounds = new double[no_vars];
    fill_n(prob->lower_bounds, no_vars, lb);
    prob->upper_bounds = new double[no_vars];
    fill_n(prob->upper_bounds, no_vars, ub);
    prob->info = NULL;
    prob->evaluate = fe;
    prob->validate =  validateProbTruncating;
}

void evalZDT1(void* p, const double x[], double y[]) { ZDT1(x, ((const MOProblem*)p)->no_vars, y); }
void evalZDT2(void* p, const double x[], double y[]) { ZDT2(x, ((const MOProblem*)p)->no_vars, y); }
void evalZDT3(void* p, const double x[], double y[]) { ZDT3(x, ((const MOProblem*)p)->no_vars, y); }
void evalZDT4(void* p, const double x[], double y[]) { ZDT4(x, ((const MOProblem*)p)->no_vars, y); }
void evalZDT6(void* p, const double x[], double y[]) { ZDT6(x, ((const MOProblem*)p)->no_vars, y); }


MOProblem* newZDTProblem(const char* name, const int no_vars)
{
    MOProblem* prob = new MOProblem();
    if (StrEq(name, "ZDT1"))
    {
        initAnyZDT(prob, name, no_vars, evalZDT1, 0.0, 1.0);
    }
    else if (StrEq(name, "ZDT2"))
    {
        initAnyZDT(prob, name, no_vars, evalZDT2, 0.0, 1.0);
    }
    else if (StrEq(name, "ZDT3"))
    {
        initAnyZDT(prob, name, no_vars, evalZDT3, 0.0, 1.0);
    }
    else if (StrEq(name, "ZDT4"))
    {
        initAnyZDT(prob, name, no_vars, evalZDT4, -5.0, 5.0);
        prob->lower_bounds[0] = 0.0;
        prob->upper_bounds[0] = 1.0;
    }
    else if (StrEq(name, "ZDT6"))
    {
        initAnyZDT(prob, name, no_vars, evalZDT6, 0.0, 1.0);
    }
    else
    {
        delete prob;
        return NULL;
    }
    return prob;
}

void freeZDTProblem(MOProblem* prob)
{
    freeMOProblem(prob);
}
*/
