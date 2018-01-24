#include "nd_archive.h"

bool NdArchive::contain(const Solution& s) const
{
    const int no_vars = s.noVars();
    return std::any_of(solutions.begin(), solutions.begin() + act_size,
                       [&s,no_vars,this](const Solution& a) {
        return euclideanDist(s.vars, a.vars, no_vars) < this->eps;
    });
}

bool NdArchive::add(const Solution& s)
{
    throw_assert(act_size < max_size, "min distance archive is full (reached "
                                      "max size of " << max_size << ')');
    const size_t pos = act_size;
    solutions.resize(++act_size);
    solutions[pos] = s;

    // first solution is always included
    if (nd_set.empty())
    {
        nd_set.insert(pos);
        return true;
    }

    for (auto it = nd_set.begin(); it != nd_set.end();)
    {
        if (*it == pos)
            continue;
        else if (solutions[pos] < solutions[*it]) // erase dominated solutions
            it = nd_set.erase(it);
        else if (solutions[*it] < solutions[pos]) // if is dominated, dont include
            return false;
        else
            it++;
    }

    // is not dominated
    nd_set.insert(pos);
    return true;
}

size_t NdArchive::add(const Population& pop)
{
    size_t no_nd = 0;
    for (const auto& s : pop)
        no_nd += int(add(s));
    return no_nd;
}

