#include "population.h"

/*
Solution* newPopulationForProblem(const int pop_size, MOProblem* prob)
{
    error(prob == NULL, "unknown problem");
    const int no_vars = prob->no_vars;
    const int no_objs = prob->no_objs;

    Solution* pop = new Solution[pop_size];
    Solution::VarsT** vars = newMatrix<Solution::VarsT>(pop_size, no_vars);
    Solution::ObjsT** objs = newMatrix<Solution::ObjsT>(pop_size, no_objs);

    std::fill_n(objs[0], pop_size * no_objs, nan(""));
    std::fill_n(vars[0], pop_size * no_vars, nan(""));

    for (int i = 0; i < pop_size; i++)
    {
        pop[i].idx = i;
        pop[i].vars = vars[i];
        pop[i].objs = objs[i];
        pop[i].problem = prob;
    }

    // deallocate indexes
    delete[] vars;
    delete[] objs;

    return pop;
}

Solution* newSampledPopulation(const int pop_size, MOProblem* prob,
                               const SamplingMethod method, Solution* pop)
{
    if (pop == NULL)
        pop = newPopulationForProblem(pop_size, prob);

    const int no_vars = prob->no_vars;
    const MOProblem::dvector& lower_bounds = prob->lower_bounds;
    const MOProblem::dvector& upper_bounds = prob->upper_bounds;

    sample(method, pop_size, no_vars, lower_bounds.data(), upper_bounds.data(),
           pop[0].vars);

    return pop;
}

void freePopulationContent(Solution* pop)
{
    // deallocate variables
    delete[] pop[0].vars;
    delete[] pop[0].objs;
}

void freePopulation(Solution* pop)
{
    freePopulationContent(pop);
    delete[] pop;
}

void printPopulation(FILE* out, const Solution* pop, const int pop_size)
{
    printPopulationPrefix(out, pop, pop_size, "");
}

void printPopulationPrefix(FILE* out, const Solution* pop, const int pop_size,
                           const char* prefix)
{
    for (int i = 0; i < pop_size; i++)
    {
        printSolutionPrefix(out, &pop[i], prefix);
        fprintf(out, "\n");
    }
}

void printPopulationObjs(FILE* out, const Solution* pop, const int pop_size)
{
    for (int i = 0; i < pop_size; i++)
    {
        const int no_objs = pop[i].noObjs();
        printDblVector(out, pop[i].objs, no_objs);
        fprintf(out, "\n");
    }
}

int readPopulation(FILE* in, const int pop_size, Solution* pop)
{
    error(pop == NULL, "population must be allocated");
    int i = 0;
    while ((i < pop_size) && (!feof(in)))
    {
        if (readSolution(in, &pop[i]) != NULL)
            i++; // skip wrong format
        while ((fgetc(in) != '\n') && (!feof(in)));
    }
    return i;
}

Solution* readPopulationForProblem(FILE* in, MOProblem* prob,
                                    const int pop_size, Solution* pop)
{
    if (pop == NULL)
        pop = newPopulationForProblem(pop_size, prob);
    readPopulation(in, pop_size, pop);
    return pop;
}

FILE* initPlotPopulationProc()
{
    FILE* pplot = popen("gnuplot -noraise -persist", "w");
//    fprintf(pplot, "set logscale xy\n");
    return pplot;
}

void deinitPlotPopulationProc(FILE** pplot)
{
    fprintf(*pplot, "exit\n");
    pclose(*pplot);
    *pplot = NULL;
}

//#define RESCALE

static void plotPopulationCmd2D(FILE* pout, Solution const pop[],
                                const int pop_size)
{
    FILE* tmp_file = fopen(g_plot_filename, "w");
    for (int i = 0; i < pop_size; i++)
    {
#ifdef RESCALE
        // scale = 1 - 1 / (obj + 1) => obj = -(scale - 1)
        fprintf(tmp_file, "%.6lf %.6lf\n",
                1.0 / (1.0 - pop[i].objs[0]) - 1.0,
                1.0 / (1.0 - pop[i].objs[1]) - 1.0);
#else
        fprintf(tmp_file, "%.6lf %.6lf\n", pop[i].objs[0], pop[i].objs[1]);
#endif
    }
    fclose(tmp_file);
    fprintf(pout, "plot '%s'\n", g_plot_filename);
    fflush(pout);
}

static void plotPopulationCmd3D(FILE* pout, Solution const pop[],
                                const int pop_size)
{
    FILE* tmp_file = fopen(g_plot_filename, "w");
    for (int i = 0; i < pop_size; i++)
        fprintf(tmp_file, "%.6lf %.6lf %.6lf\n", pop[i].objs[0], pop[i].objs[1],
                pop[i].objs[2]);
    fclose(tmp_file);
    fprintf(pout, "splot '%s'\n", g_plot_filename);
    fflush(pout);
}

void plotPopulationCmd(FILE* pout, const Solution pop[], const int pop_size)
{
    if (pop[0].noObjs() == 2)
        plotPopulationCmd2D(pout, pop, pop_size);
    else if (pop[0].noObjs() >= 3)
        plotPopulationCmd3D(pout, pop, pop_size);
    else
        error(true, "problem must have at least 2 objectives");
}

void plotPopulation(Solution const pop[], const int pop_size)
{
    FILE* p = initPlotPopulationProc();
    plotPopulationCmd(p, pop, pop_size);
    printf("Population plotted. Press any key to continue...\n");
    fflush(stdout);
    getchar();
    deinitPlotPopulationProc(&p);
}

void evaluatePopulation(Solution pop[], const int pop_size)
{
    for (int i = 0; i < pop_size; i++)
        pop[i].evaluate();
}

void validatePopulation(Solution pop[], const int pop_size)
{
    for (int i = 0; i < pop_size; i++)
        pop[i].validate();
}

void popMinMaxVars(const Solution pop[], const int pop_size,
                   Solution::VarsT min_vars[],
                   Solution::VarsT max_vars[])
{
    const int no_vars = pop[0].noVars();

    std::copy_n(pop[0].vars, no_vars, min_vars);
    std::copy_n(pop[0].vars, no_vars, max_vars);

    for (int i = 1; i < pop_size; ++i)
    {
        for (int j = 0; j < no_vars; ++j)
        {
            const double v = pop[i].vars[j];
            if (v < min_vars[j])
                min_vars[j] = v;
            else if (v > max_vars[j])
                max_vars[j] = v;
        }
    }
}

void popIdealNadir(Solution const pop[], const int pop_size,
                   Solution::ObjsT ideal[],
                   Solution::ObjsT nadir[])
{
    const int no_objs = pop[0].noObjs();

    std::copy_n(pop[0].objs, no_objs, ideal);
    std::copy_n(pop[0].objs, no_objs, nadir);

    for (int i = 1; i < pop_size; ++i)
    {
        for (int j = 0; j < no_objs; ++j)
        {
            const double v = pop[i].objs[j];
            if (v < ideal[j])
                ideal[j] = v;
            else if (v > nadir[j])
                nadir[j] = v;
        }
    }
}

Solution::VarsT* getVarsArr(Solution const pop[], const int pop_size,
                            Solution::VarsT *vars)
{
    const int no_vars = pop[0].noVars();
    if (vars == NULL)
        vars = new double[pop_size * no_vars];
    for (int i = 0; i < pop_size; i++)
        std::copy_n(pop[i].vars, no_vars, vars + i * no_vars);
    return vars;
}

Solution::ObjsT* getObjsArr(Solution const pop[], const int pop_size,
                            Solution::ObjsT objs[])
{
    const int no_objs = pop[0].noObjs();
    if (objs == NULL)
        objs = new double[pop_size * no_objs];
    for (int i = 0; i < pop_size; i++)
        std::copy_n(pop[i].objs, no_objs, objs + i * no_objs);
    return objs;
}

Solution::ObjsT* getObjArr(Solution const pop[], const int pop_size,
                           const int obj, Solution::ObjsT objs[])
{
    if (objs == NULL)
        objs = new double[pop_size];
    for (int i = 0; i < pop_size; i++)
        objs[i] = pop[i].objs[obj];
    return objs;
}


Solution::VarsT** getVars(Solution const pop[], const int pop_size,
                          Solution::VarsT** vars)
{
    const int no_vars = pop[0].noVars();
    if (vars == NULL)
        vars = newMatrix<Solution::VarsT>(pop_size, no_vars);
    for (int i = 0; i < pop_size; i++)
        std::copy_n(pop[i].vars, no_vars, vars[i]);
    return vars;
}

Solution::ObjsT** getObjs(Solution const pop[], const int pop_size,
                          Solution::ObjsT** objs)
{
    const int no_objs = pop[0].noObjs();
    if (objs == NULL)
        objs = newMatrix<Solution::ObjsT>(pop_size, no_objs);
    for (int i = 0; i < pop_size; i++)
        std::copy_n(pop[i].objs, no_objs, objs[i]);
    return objs;
}

Solution::VarsT* getVarsMean(const Solution* pop, const int pop_size,
                             Solution::VarsT mean[])
{
    const int no_vars = pop[0].noVars();
    if (mean == NULL)
        mean = new Solution::VarsT[no_vars];
    for (int i = 0; i < pop_size; i++)
    {
        for (int j = 0; j < no_vars; ++j)
            mean[j] += pop[i].vars[j];
    }
    for (int j = 0; j < no_vars; ++j)
        mean[j] /= pop_size;
    return mean;
}

Solution::VarsT* getVarsStd(const Solution* pop, const int pop_size,
                            const Solution::VarsT mean[],
                            Solution::VarsT std[])
{
    const int no_vars = pop[0].noVars();
    if (std == NULL)
        std = new Solution::VarsT[no_vars];
    for (int i = 0; i < pop_size; i++)
    {
        for (int j = 0; j < no_vars; ++j)
            std[j] += sqr(pop[i].vars[j] - mean[j]);
    }
    for (int j = 0; j < no_vars; ++j)
        std[j] = sqrt(std[j] / pop_size);
    return std;
}

Solution::ObjsT* getObjsMean(const Solution* pop, const int pop_size,
                             Solution::ObjsT mean[])
{
    const int no_objs = pop[0].noObjs();
    if (mean == NULL)
        mean = new Solution::ObjsT[no_objs];
    for (int i = 0; i < pop_size; i++)
    {
        for (int j = 0; j < no_objs; ++j)
            mean[j] += pop[i].objs[j];
    }
    for (int j = 0; j < no_objs; ++j)
        mean[j] /= pop_size;
    return mean;
}

Solution::ObjsT* getObjsStd(const Solution* pop, const int pop_size,
                            const Solution::ObjsT mean[], Solution::ObjsT std[])
{
    const int no_objs = pop[0].noObjs();
    if (std == NULL)
        std = new Solution::ObjsT[no_objs];
    for (int i = 0; i < pop_size; i++)
    {
        for (int j = 0; j < no_objs; ++j)
            std[j] += sqr(pop[i].objs[j] - mean[j]);
    }
    for (int j = 0; j < no_objs; ++j)
        std[j] = sqrt(std[j] / pop_size);
    return std;
}

void popNormalizeVars(Solution* pop, const int pop_size,
                      const Solution::VarsT lower,
                      const Solution::VarsT upper)
{
    const int no_vars = pop[0].noVars();
    Solution::VarsT vmin[no_vars], vmax[no_vars];
    popMinMaxVars(pop, pop_size, vmin, vmax);
    popNormalizeVarsIn(pop, pop_size, vmin, vmax, lower, upper);
}

void popNormalizeVarsIn(Solution* pop, const int pop_size,
                        const Solution::VarsT* vmin, const Solution::VarsT* vmax,
                        const Solution::VarsT lower, const Solution::VarsT upper)
{
    static const Solution::VarsT eps = 1e-12;
    const int no_vars = pop[0].noVars();
    for (int i = 0; i < no_vars; ++i)
    {
        const Solution::VarsT delta = vmax[i] - vmin[i],
                s = (delta >= eps) ? (upper - lower) / delta : 0;
        for (int j = 0; j < pop_size; ++j)
            pop[j].vars[i] = (pop[j].vars[i] - vmin[i]) * s + lower;
    }
}

void popNormalizeObjs(Solution* pop, const int pop_size,
                      const Solution::ObjsT lower,
                      const Solution::ObjsT upper)
{
    const int no_objs = pop[0].noObjs();
    Solution::ObjsT ideal[no_objs], nadir[no_objs];
    popIdealNadir(pop, pop_size, ideal, nadir);
    popNormalizeObjsIn(pop, pop_size, ideal, nadir, lower, upper);

}

void popNormalizeObjsIn(Solution* pop, const int pop_size,
                        const Solution::ObjsT* omin, const Solution::ObjsT* omax,
                        const Solution::ObjsT lower, const Solution::ObjsT upper)
{
    static const double eps = 1e-12;
    const int no_objs = pop[0].noObjs();
    for (int i = 0; i < no_objs; ++i)
    {
        const Solution::ObjsT delta = omax[i] - omin[i],
                a = (delta >= eps) ? (upper - lower) / delta : 0,
                b = - a * omin[i] + lower;
        for (int j = 0; j < pop_size; ++j)
            pop[j].objs[i] = pop[j].objs[i] * a + b;
    }
}

void shufflePopulation(Solution* pop, const int pop_size)
{
    for (int i = pop_size - 1; i >= 0; --i)
    {
        const int rnd = RNG::intUniform(i);
        swapSolutions(&pop[i], &pop[rnd]);
    }
}
*/
void sortPopulationBy(Solution pop[], const int pop_size, int ranks[])
{
    for (int i = 0; i < pop_size; i++)
    {
        int min_rank = ranks[i], min_rank_idx = i;

        for (int j = i + 1; j < pop_size; j++)
        {
            if (ranks[j] < min_rank)
            {
                min_rank = ranks[j];
                min_rank_idx = j;
            }
        }

        if (min_rank_idx != i)
        {
            swap(pop[i], pop[min_rank_idx]);
            std::swap(ranks[i], ranks[min_rank_idx]);
        }
    }
}


void Population::shuffle(Population::sols_iterator first,
                         Population::sols_iterator last)
{
    for (auto i = last - first - 1; i > 0; --i)
    {
        std::uniform_int_distribution<decltype(i)> d(0,i);
        Solution& a = first[i];
        Solution& b = first[d(RNG::engine)];
        std::swap_ranges(a.vars, a.vars + noVars(), b.vars);
        std::swap_ranges(a.objs, a.objs + noObjs(), b.objs);
    }

    //std::shuffle(from, to, RNG::engine);
}

double Population::hypervolume(const Solution::ObjsT *ref)
{
    std::vector<Solution*> ref_nd;
    AllocSolution ref_sol(*(sols[0].problem));
    ref_sol.setObjs(ref);
    for (auto& sol : sols)
        if (sol < ref_sol)
            ref_nd.push_back(&sol);
    const size_t n = ref_nd.size();
    if (n == 0) return 0;
    ematrixd ref_nd_objs(n, noObjs());
    for (int i = 0; i < ref_nd.size(); i++)
        ref_nd[i]->getObjs(ref_nd_objs.row(i).data());
    return hv(ref_nd_objs.data(), noObjs(), n, ref);
}

void Population::nonDominated(std::vector<int> &fronts) {
    int assigneds = 0;
    std::vector<bool> dominated(pop_size, false);
    std::fill(fronts.begin(), fronts.end(), -1);
    for (int l = 0; (l < size()) && (assigneds < size()); l++)
    {
        for (int i = 0; i < size(); i++)
        {
            dominated[i] = false;
            for (int j = 0; (j < i) && !dominated[i]; j++)
            {
                if (fronts[j] == -1)
                    dominated[i] = sols[j] < sols[i];
            }
            for (int j = i + 1; (j < size()) && !dominated[i]; j++)
            {
                if (fronts[j] == -1)
                    dominated[i] = sols[j] < sols[i];
            }
        }
        for (int i = 0; i < size(); i++)
        {
            if (!dominated[i] && (fronts[i] == -1))
            {
                fronts[i] = l;
                assigneds++;
            }
        }
    }
}
