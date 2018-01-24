#include "../surrogate/ml.h"
#include "../surrogate/aux.h"
#include "../surrogate/plotter.h"

#include "helpers.h"
#include <cmath>
#include <numeric>
#include <iostream>
using namespace std;

using std::hypot;

#include "Eigen/Core"
#include "Eigen/Dense"

TEST(KMeans, Apply) {
    KMeans::matrix data(2, 8);
    data <<
       0, 9, 1, 9, 1, 0, 8, 8,
       0, 8, 0, 9, 1, 1, 8, 9;
    KMeans km(data.rows(), 2);
    km.apply(data);
    auto near = [](double v, double e) { return v - e < 0.0001; };
    auto cts = km.centers;
    ASSERT_TRUE((near(cts(0,0), 0.5) && near(cts(0,1), 8.5) &&
                 near(cts(1,0), 0.5) && near(cts(1,1), 8.5)) ||
                (near(cts(0,0), 8.5) && near(cts(0,1), 0.5) &&
                 near(cts(1,0), 8.5) && near(cts(1,1), 0.5)));
}

TEST(RegELM, TrainPredict) {
    RNG::seed();
    const int no_ins = 2;
    const int no_outs = 1;
    const int no_hidden = 300;
    const int no_samples = 1000;
    const double C = std::pow(2.0, -0);
    const auto act_fun = SLFN::ELMActFunc::MULTIQUADRIC;
    ematrixd init_x(no_samples, no_ins);
    ematrixd init_y(no_samples, no_outs);
    svectord lb = {-1, -1};
    svectord ub = {1, 1};

    auto sombrero = [](double a, double b) -> double {
        const double scale = 10.0;
        const double h = hypot(a * scale, b * scale);
        return h == 0 ? 1 : sin(h) / h;
    };

    init_x.setRandom(no_samples, 2);
    for (int i = 0; i < no_samples; ++i)
        init_y(i, 0) = sombrero(init_x(i, 0), init_x(i, 1));

    RegELM elm(no_ins, no_hidden, no_outs, act_fun, C, true, true, lb, ub);

    elm.train(init_x, init_y);

    init_x.setRandom(no_samples, 2);
    for (int i = 0; i < no_samples; ++i)
        init_y(i, 0) = sombrero(init_x(i, 0), init_x(i, 1));

    ematrixd pred_y(no_samples, 1);
    elm.predictAll(init_x, pred_y);

    auto mse = (pred_y - init_y).squaredNorm() / no_samples;

    cout << mse << endl;
    EXPECT_TRUE(mse < 0.02);

    //xy_pred.col(2) << init_y;
    //Plotter plt2;
    //plt2.plot(xy_pred.data(), no_samples, 3);
//    ematrixd xy_pred(no_samples, 3);
//    xy_pred << init_x, pred_y;
//    cout << "x1 | x2 | pred_y\n" << xy_pred << endl;
//    Plotter plt;
//    plt.plot(xy_pred.data(), no_samples, 3);
}

TEST(RBF, TrainPredict) {
    using namespace SLFN;
    //ematrixd input_t(2, 8);
    //input_t <<
    //   0, 9, 1, 9, 1, 0, 8, 8,
    //   0, 8, 0, 9, 1, 1, 8, 9;
    //ematrixd input = input_t.transpose();
    //ematrixd output(8, 1);

    //RBF rbf(2, 2, 1);

    //rbf.train(input, output);

    const int no_ins = 2;
    const int no_outs = 1;
    const int no_hidden = 100;
    const int no_samples = 1000;
    const double C = std::pow(2.0, -0);
    const auto act_fun = SLFN::ELMActFunc::MULTIQUADRIC;
    ematrixd init_x(no_samples, no_ins);
    ematrixd init_y(no_samples, no_outs);
    svectord lb = {-1, -1};
    svectord ub = {1, 1};

    auto sombrero = [](double a, double b) -> double {
        const double scale = 10.0;
        const double h = hypot(a * scale, b * scale);
        return h == 0 ? 1 : sin(h) / h;
    };

    init_x.setRandom(no_samples, 2);
    for (int i = 0; i < no_samples; ++i)
        init_y(i, 0) = sombrero(init_x(i, 0), init_x(i, 1));

    RBF rbf(no_ins, no_hidden, no_outs, lb, ub);

    rbf.train(init_x, init_y);

    init_x.setRandom(no_samples, 2);
    for (int i = 0; i < no_samples; ++i)
        init_y(i, 0) = sombrero(init_x(i, 0), init_x(i, 1));

    ematrixd pred_y(no_samples, 1);
    rbf.predictAll(init_x, pred_y);

    auto mse = (pred_y - init_y).squaredNorm() / no_samples;

    cout << mse << endl;
    EXPECT_TRUE(mse < 0.02);

    //ematrixd xy_pred(no_samples, 3);
    //xy_pred << init_x, pred_y;
    //Plotter plt2;
    //plt2.plot(xy_pred.data(), no_samples, 3);

    //cout << "x1 | x2 | pred_y\n" << xy_pred << endl;
    //xy_pred.col(2) << init_y;
    //Plotter plt;
    //plt.plot(xy_pred.data(), no_samples, 3)
}

/*
TEST(ELMTrainPredict, OSELM) {
    const int no_ins = 2;
    const int no_outs = 1;
    const int no_hidden = 10;
    const auto act_fun = ELMActFunc::SIGMOID;
    ematrixd init_x(50, 2);
    init_x.setRandom(50, 2);
    evectord init_y(50);
    svectord inu = {1, 1};
    svectord inl = {-1, -1};
    double res[1];

    auto sombrero = [](double a, double b) -> double {
        const double h = hypot(a, b);
        return sin(h) / h;
    };

    for (int i = 0; i < 50; ++i)
        init_y(i) = sombrero(init_x(i, 0), init_x(i, 1));

    OnSeqELM elm(no_ins, no_hidden, no_outs, act_fun, init_x.data(),
                 init_y.data(), 50, inu, inl, true, true);

    init_x.setRandom(50, 2);
    for (int i = 0; i < 50; ++i)
    {
        res[0] = sombrero(init_x(i, 0), init_x(i, 1));
        elm.train(init_x.row(i).data(), 1, res);
    }

    init_x.setRandom(50, 2);
    for (int i = 0; i < 50; ++i)
    {
        const double h = sombrero(init_x(i, 0), init_x(i, 1));
        elm.predict(init_x.row(i).data(), 1, res);

        cout << h <<  " ~> " << res[0] << "\n";
    }
}
*/
/*
TEST(ELMModel, GaussianFunction) {
    ematrixd w(3, 2);
    w << 1, -1,
         2, -2,
         3, -3;
    evectord b(1);
    b << 4;
    evectord x(3);
    x << 0,
         2,
        -1;
    evectord result = ELMModel::gaussianAF(x, w, b);

    ASSERT_TRUE(result.rows() == 2 && result.cols() == 1);
    EXPECT_DOUBLE_EQ(result(0,0), exp(-4 * 17));
    EXPECT_DOUBLE_EQ(result(1,0), exp(-4 * 21));
};*/
/*
class ElmTests : public testing::Test
{
protected:
    function<void(double*,double*)> test_func;
    int no_ins;
    int no_outs;
    int no_samples;
    int no_tests;
    double* ins;
    double* outs;
    double* test_ins;
    double* test_outs;

    double urand() {
        return RNG::realUniform<double>();
    };

    virtual void SetUp() {
        test_func = sombrero;
        no_ins = 2;
        no_outs = 1;
        no_samples = 50;
        no_tests = 100;

        srand(123);

        ins = new double[no_ins * no_samples];
        outs = new double[no_outs * no_samples];
        for (int i = 0; i < no_samples; ++i)
        {
            double* in = ins + i * no_ins;
            for (int j = 0; j < no_ins; ++j)
                in[j] = urand() * 20 - 10.0;
            test_func(in, outs + i * no_outs);
        }

        test_ins = new double[no_ins * no_tests];
        test_outs = new double[no_outs * no_tests];
        for (int i = 0; i < no_tests; ++i)
        {
            double* in = test_ins + i * no_ins;
            for (int j = 0; j < no_ins; ++j)
                in[j] = urand() * 20 - 10.0;
            test_func(in, test_outs + i * no_outs);
        };
    };

    virtual void TearDown() {
        delete[] ins;
        delete[] outs;
        delete[] test_ins;
        delete[] test_outs;
    };

    static void sombrero(double* x, double* z) {
        const double h = hypot(x[0], x[1]);
        z[0] = sin(h) / h;
    };
};

TEST_F(ElmTests, All) {
    RegELM elm(no_ins, 500, 1, 2);
    elm.train(ins, no_samples, outs);

    double mse = 0.0;
    double* predicts = elm.predict(test_ins, no_tests);
    for (int i = 0; i < no_tests; ++i)
    {

        //for (int j = 0; j < no_ins; ++j) {
        //    cout << test_ins[i * no_ins + j] << " ";
        //}
        cout << test_outs[i] << " ~> " << predicts[i] << endl;

        const double diff = predicts[i] - test_outs[i];
        mse += diff * diff;
        cout << diff << endl;
        //EXPECT_NEAR(predicts[i], test_outs[i], 0.25);
    }
    cout << "mse: " << mse / no_tests << endl;

    delete[] predicts;
}
*/
#include "../surrogate/ml.h"
