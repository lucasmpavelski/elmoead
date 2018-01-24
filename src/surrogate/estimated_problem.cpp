#include "estimated_problem.h"



AdaptMethod str2AdaptMethod(const std::string &s)
{
    if (s == "NONE") return AdaptMethod::NONE;
    if (s == "OFFLINE_INIT") return AdaptMethod::OFFLINE_INIT;
    if (s == "ONLINE_MSE") return AdaptMethod::ONLINE_MSE;
    throw_assert(false, "adapt_method must be (NONE | OFFLINE_INIT | ONLINE_MSE)");
    return AdaptMethod::NONE;
}
