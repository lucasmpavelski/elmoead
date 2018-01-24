#include "../surrogate/plotter.h"

TEST(PlotterTests, plot1d) {
    double data[] = {1, 3, 2};
    Plotter p(data, 3, 1);
};

TEST(PlotterTests, plot2d) {
    double data[] = {1, 3,
                     2, 4,
                     1, 6};
    Plotter p(data, 3, 2);
};

TEST(PlotterTests, plot3d) {
    double data[] = {1, 3, 2,
                     2, 4, 1,
                     1, 6, 5};
    Plotter p(data, 3, 3);
};

TEST(PlotterTests, plot4d) {
    double data[] = {1, 3, 2, 3,
                     2, 4, 1, 1,
                     1, 6, 5, 2};
    Plotter p(data, 3, 4);
};
