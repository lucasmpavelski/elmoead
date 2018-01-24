#include "../surrogate/aux.h"

#include "helpers.h"

TEST(SlidingWindow, Iteration)
{
    SlidingWindow<int> sl(3, 0);

    sl.append(1);
    sl.append(2);
    sl.append(3);
    sl.append(4);
    sl.append(5);

    int expect[3] = {3, 4, 5};
    int i = 0;
    for (const auto& el : sl)
        EXPECT_EQ(el, expect[i++]);
}
