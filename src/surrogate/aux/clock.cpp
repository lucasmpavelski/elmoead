#include "clock.h"

Clock::time_point Clock::start;

void Clock::begin() {
    start = system_clock::now();
}

Clock::duration Clock::tellTime() {
    return system_clock::now() - start;
}
