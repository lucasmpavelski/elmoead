#include "rand.h"

std::mt19937_64 RNG::engine;
std::random_device RNG::true_rand_engine;

long RNG::seed(long s) {
    engine.seed(s);
    srand(s);
    return s;
}
