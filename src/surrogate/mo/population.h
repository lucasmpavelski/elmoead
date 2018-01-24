#ifndef POPULATION_H
#define POPULATION_H

#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>

#include "hv.h"

#include "../aux.h"
#include "../plotter.h"

#include "mo_problem.h"
#include "solution.h"

class Population {
public:
    using VarsT = Solution::VarsT;
    using ObjsT = Solution::ObjsT;
    using sols_container = std::vector<Solution>;
    using sols_iterator = sols_container::iterator;
    using const_sols_iterator = sols_container::const_iterator;
    using vars_container = Eigen::Matrix<VarsT,
                                         Eigen::Dynamic, Eigen::Dynamic,
                                         Eigen::RowMajor>;
    using objs_container = Eigen::Matrix<ObjsT,
                                         Eigen::Dynamic, Eigen::Dynamic,
                                         Eigen::RowMajor>;

private:
    size_t pop_size, no_vars, no_objs;
    sols_container sols;
    vars_container vars_data;
    objs_container objs_data;

public:
    Population(const size_t& n, MOProblem& prob,
               const SamplingMethod& method = NO_SAMPLING)
        : pop_size(n)
        , no_vars(prob.noVars())
        , no_objs(prob.noObjs())
        , sols(n)
        , vars_data(n, prob.noVars())
        , objs_data(n, prob.noObjs())
    {
        throw_assert(n > 0, "invalid population size of " << n);
        initSols(&prob);
        sampleVars(method, prob.lower_bounds.data(), prob.upper_bounds.data());
    }

    Population(MOProblem& prob, const SamplingMethod& method = NO_SAMPLING)
        : Population{ 1, prob, method }
    {
    }

    Population(const Population& other_pop)
        : Population{ *other_pop[0].problem }
    {
        resize(other_pop.size());
        vars_data = other_pop.vars_data;
        objs_data = other_pop.objs_data;
    }

    virtual ~Population(){}

    size_t size()   const { return pop_size; }
    size_t noVars() const { return no_vars; }
    size_t noObjs() const { return no_objs; }

    Solution* solsData() { return sols.data(); }
    const Solution* solsData() const { return sols.data(); }

    Solution::VarsT* varsDataArr() { return vars_data.data(); }
    Solution::ObjsT* objsDataArr() { return objs_data.data(); }
    const Solution::VarsT* varsDataArr() const { return vars_data.data(); }
    const Solution::ObjsT* objsDataArr() const { return objs_data.data(); }

    vars_container& vars() { return vars_data; }
    objs_container& objs() { return objs_data; }
    const vars_container& vars() const { return vars_data; }
    const objs_container& objs() const { return objs_data; }

    bool dataIsContiguous()
    {
        for (int i = 0; i < size(); ++i)
            if (sols[i].vars != varsDataArr() + i * noVars() || sols[i].objs != objsDataArr() + i * noObjs())
                return false;
        return true;
    }

    const Solution& operator[](const size_t& i) const { return sols[i]; }
    Solution& operator[](const size_t& i) { return sols[i]; }

    const const_sols_iterator begin() const { return sols.cbegin(); }
    sols_iterator begin() { return sols.begin(); }

    const const_sols_iterator end() const { return sols.cend(); }
    sols_iterator end() { return sols.end(); }

    MOProblem* getProblem() const { sols[0].problem; }

    void setProblem(MOProblem* prob)
    {
        for (auto& s : sols)
            s.problem = prob;
    }

    void resize(const size_t& n)
    {
        throw_assert(n > 0, "invalid population size of: " << n);
        pop_size = n;
        MOProblem* p = sols[0].problem;
        sols.resize(n);
        vars_data.conservativeResize(n, noVars());
        objs_data.conservativeResize(n, noObjs());
        initSols(p);
    }

    void sampleVars(const SamplingMethod method, const VarsT* lb,
                    const VarsT* ub)
    {
        sampleVars(method, begin(), end(), lb, ub);
    }

    void sampleVars(const SamplingMethod method, sols_iterator from,
                    sols_iterator to, const VarsT* lb, const VarsT* ub)
    {
        const size_t b = std::distance(begin(), from);
        const size_t n = std::distance(from, to);
        VarsT* dt = varsDataArr() + b * noVars();
        sample(method, n, noVars(), lb, ub, dt);
    }

    void shuffle() { shuffle(begin(), end()); }
    void shuffle(sols_iterator first, sols_iterator last);

    void validateAll()
    {
        for (auto& s : sols)
            s.validate();
    }

    void evaluateAll()
    {
        for (auto& s : sols)
            s.evaluate();
    }

    void validateAndEvaluateAll()
    {
        validateAll();
        evaluateAll();
    }

    void plotObjs(Plotter& p) const { p.plot(objsDataArr(), size(), noObjs()); }
    void plotVars(Plotter& p) const { p.plot(varsDataArr(), size(), noVars()); }

    double hypervolume(const Solution::ObjsT* ref);

    void nonDominated(std::vector<int>& fronts);

    template <typename Compare>
    void sort(Compare func)
    {
        for (int i = 0; i < size(); i++) {
            int min_idx = i;
            for (int j = i + 1; j < size(); j++)
                if (func(sols[j], sols[min_idx]))
                    min_idx = j;
            if (min_idx != i)
                swap(sols[i], sols[min_idx]);
        }
    }

    template <typename SeqT>
    void sort(SeqT ranks[])
    {
        for (int i = 0; i < size(); i++) {
            SeqT min_rank = ranks[i];
            int min_rank_idx = i;
            for (int j = i + 1; j < size(); j++) {
                if (ranks[j] < min_rank) {
                    min_rank = ranks[j];
                    min_rank_idx = j;
                }
            }
            if (min_rank_idx != i) {
                swap(sols[i], sols[min_rank_idx]);
                std::swap(ranks[i], ranks[min_rank_idx]);
            }
        }
    }

    int nonDominatedSort()
    {
        std::vector<int> fronts(pop_size);
        return nonDominatedSort(fronts);
    }

    int nonDominatedSort(std::vector<int>& fronts)
    {
        fronts.resize(pop_size);
        nonDominated(fronts);
        sort(fronts.data());
        size_t no_nd = 0;
        while ((no_nd < pop_size) && (fronts[no_nd] == 0)) no_nd++;
        return no_nd;
    }

    friend std::ostream& operator<<(std::ostream& os, const Population& pop)
    {
        os << pop.size() << "\n";
        for (const auto& s : pop.sols)
            os << s << "\n";
        return os;
    }

    friend std::istream& operator>>(std::istream& is, Population& pop)
    {
        std::string size_s;
        std::getline(is, size_s);
        pop.resize(std::stoi(size_s));
        for (auto& s : pop.sols)
            is >> s;
        return is;
    }

    friend bool operator==(const Population& popa, const Population& popb)
    {
        return popa.size() == popb.size() && popa.vars_data == popb.vars_data && popa.objs_data == popb.objs_data;
    }

private:
    void initSols(MOProblem* prob)
    {
        for (int i = 0; i < size(); ++i) {
            sols[i].idx = i;
            sols[i].problem = prob;
            sols[i].vars = vars_data.row(i).data();
            sols[i].objs = objs_data.row(i).data();
        }
    }
};

template <>
struct SerializedWriter<const Population> {
    SerializedWriter(const Population& dt)
        : dt(dt){};

    friend std::ostream& operator<<(std::ostream& os,
                                    const SerializedWriter<const Population>& ser_pop)
    {
        const Population& pop = ser_pop.dt;
        size_t size = pop.size();
        os << serialize(size);
        for (int i = 0; i < pop.size() * pop.noVars(); i++)
            os << serialize(pop.varsDataArr()[i]);
        for (int i = 0; i < pop.size() * pop.noObjs(); i++)
            os << serialize(pop.objsDataArr()[i]);
        return os;
    };

    const Population& dt;
};

template <>
struct SerializedReader<Population> {
    SerializedReader(Population& dt)
        : dt(dt){};

    friend std::istream& operator>>(std::istream& is,
                                    SerializedReader<Population>&& ser_pop)
    {
        Population& pop = ser_pop.dt;
        size_t size;
        is >> read_serialized(size);
        pop.resize(size);
        for (int i = 0; i < pop.size() * pop.noVars(); i++)
            is >> read_serialized(pop.varsDataArr()[i]);
        for (int i = 0; i < pop.size() * pop.noObjs(); i++)
            is >> read_serialized(pop.objsDataArr()[i]);
        return is;
    };

    Population& dt;
};

#endif /* POPULATION_H */
