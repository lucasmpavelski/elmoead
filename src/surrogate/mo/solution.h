#ifndef SOLUTION_H
#define SOLUTION_H

#include <iostream>
#include <algorithm>

#include "../aux.h"

#include "mo_problem.h"

struct Solution
{
    typedef double VarsT;
    typedef double ObjsT;

    MOProblem* problem;
    VarsT* vars;
    ObjsT* objs;
    int idx;

    Solution(MOProblem* problem=nullptr, VarsT* vars=nullptr,
             ObjsT* objs=nullptr, int idx = -1) noexcept :
        problem(problem),
        vars(vars),
        objs(objs),
        idx(idx)
    {}

    Solution(const Solution& s) = default;
    Solution(Solution&& s) = delete;

    Solution& operator=(const Solution& sol) noexcept
    {
        warning(problem != sol.problem,
                "assigning solutions with different problems");

        if (problem != sol.problem)
        {
            cout << problem << endl;
            cout << sol.problem << endl;
        }
        setVars(sol.vars);
        setObjs(sol.objs);
        return *this;
    }

    Solution& operator=(Solution&& sol) = delete;
    /*Solution& operator=(Solution&& sol) noexcept
    {
        warning(problem != sol.problem,
                "assigning solutions with different problems");
        vars = sol.vars;
        objs = sol.objs;
        idx = sol.idx;
        sol.vars = nullptr;
        sol.objs = nullptr;
        sol.idx = -1;
        return *this;
    }*/

    virtual ~Solution() {}

    inline void swap(Solution& sol) noexcept {
        if (this == &sol) return;
        warning(problem != sol.problem, "swapping solutions for different problems");
        std::swap_ranges(vars, vars + noVars(), sol.vars);
        std::swap_ranges(objs, objs + noObjs(), sol.objs);
    }

    inline int noVars() const { return problem? problem->noVars() : 0; }
    inline int noObjs() const { return problem? problem->noObjs() : 0; }

    template<typename Itr>
    void getVars(Itr p) const { std::copy_n(vars, noVars(), p); }
    template<typename Itr>
    void getObjs(Itr p) const { std::copy_n(objs, noObjs(), p); }

    void setVars(const VarsT* nv) { std::copy_n(nv, noVars(), vars); }
    void setObjs(const VarsT* no) { std::copy_n(no, noObjs(), objs); }

    void setVars(std::initializer_list<VarsT>&& nv) {
        std::copy_n(nv.begin(), noVars(), vars);
    }
    void setObjs(std::initializer_list<ObjsT>&& no) {
        std::copy_n(no.begin(), noObjs(), objs);
    }

    void evaluate() { problem->evaluate(vars, objs); }
    void validate() { problem->validate(vars);       }
    void validateEvaluate() { validate(); evaluate(); }

    bool equals(const Solution& other,
                const double eps=std::numeric_limits<double>::epsilon()) const {
        return equals(other, eps, eps);
    }

    bool equals(const Solution& other, const double vars_eps,
                const double objs_eps) const {
        return std::equal(vars, vars + noVars(), other.vars,
                          [vars_eps](const VarsT& av, const VarsT& bv) {
                              return std::abs(av - bv) < vars_eps; }) &&
               std::equal(objs, objs + noObjs(), other.objs,
                          [objs_eps](const ObjsT& ao, const ObjsT& bo) {
                              return std::abs(ao - bo) < objs_eps; });
    }

    bool dominates(const Solution& b) const {
        bool isworse = false, isequal = true;
        for (int i = 0; i < noObjs() && !isworse; ++i)
        {
            isworse = objs[i] > b.objs[i];
            isequal = (objs[i] == b.objs[i]) && isequal;
        }
        return !isequal && !isworse;
    }

    bool weaklyDominates(const Solution& b) const {
        bool better_or_eq = true;
        for (int i = 0; i < noObjs() && better_or_eq; ++i)
            better_or_eq &= (objs[i] <= b.objs[i]);
        return better_or_eq;
    }

    friend bool operator==(const Solution& a, const Solution& b) {
        return a.equals(b);
    }

    friend bool operator!=(const Solution& a, const Solution& b) {
        return !a.equals(b);
    }

    friend bool operator<(const Solution& a, const Solution& b) {
        return a.dominates(b);
    }

    friend bool operator<=(const Solution& a, const Solution& b) {
        return a.weaklyDominates(b);
    }

    friend std::ostream& operator<<(std::ostream& os, const Solution& s) {
        for (int i = 0; i < s.noVars(); ++i)
            os << s.vars[i] << " ";
        os << ";";
        for (int i = 0; i < s.noObjs(); ++i)
            os << s.objs[i] << " ";
        return os;
    }

    friend std::istream& operator>>(std::istream& is, Solution& s) {
        for (int i = 0; i < s.noVars(); ++i)
            is >> s.vars[i];
        char c;
        do { is >> c; } while (c != ';' && is.good());
        for (int i = 0; i < s.noObjs(); ++i)
            is >> s.objs[i];
        return is;
    };

    friend void swap(Solution& a, Solution& b) noexcept { a.swap(b); };
};

struct AllocSolution : Solution
{
    using vars_container = Eigen::Matrix<VarsT,1,Eigen::Dynamic,Eigen::RowMajor>;
    using objs_container = Eigen::Matrix<ObjsT,1,Eigen::Dynamic,Eigen::RowMajor>;
    vars_container vars_data;
    objs_container objs_data;

    AllocSolution(MOProblem& prob, const int idx=-1) :
        Solution(&prob, nullptr, nullptr, idx),
        vars_data(prob.noVars()), objs_data(prob.noObjs())
    {
        this->vars = vars_data.data();
        this->objs = objs_data.data();
    };

    virtual ~AllocSolution() {};
};

#endif /* SOLUTION_H */
