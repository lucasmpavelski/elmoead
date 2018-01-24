#pragma once

#include "population.h"
#include "../aux.h"

using PopulationCallback = std::function<void(int, const Population&)>;

inline void NullCallback(int, const Population&) {}

struct PlotCallback {
    Plotter plotter;
    size_t plot_each;

    PlotCallback(const size_t plot_each = 1)
        : plotter()
        , plot_each(plot_each) {}

    void operator()(int no_evals, const Population& cpop)
    {
        Population pop = cpop;
        size_t ns = pop.nonDominatedSort();
        pop.resize(ns);
        if (no_evals % plot_each == 0)
            pop.plotObjs(plotter);

    }
};

struct LogCallback {
    enum Mode { txt, binary } mode;
    std::string base_folder;
    const size_t freq;
    bool first_time;

    LogCallback(const std::string& base_nm, size_t freq=1, Mode mode=txt)
        : base_folder(base_nm + "_logs")
        , freq(freq)
        , mode(mode)
        , first_time{true}
    {
        std::string mkdir_cmd = "mkdir -p " + base_folder;
        throw_assert(system(mkdir_cmd.c_str()) == 0,
                     "could not create log dir");
    }

    void operator()(int no_evals, const Population& pop)
    {
        if ((no_evals % freq == 0) || first_time)
        {
            first_time = false; // always log first time
            std::ofstream out(base_folder + "/" + std::to_string(no_evals) +
                              ".evals.dat");
            double ts = Clock::tellTime().count();
            switch (mode)
            {
            case binary:
                out << serialize(ts) << serialize(pop);
                break;
            case txt:
                out << ts << '\n' << pop.objs();
                break;
            }
            out.close();
        }
    }

};
