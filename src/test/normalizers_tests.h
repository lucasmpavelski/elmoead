#include "../surrogate/ml/normalizers.h"
#include "../surrogate/aux.h"

#include "helpers.h"

TEST(NonLinearNormalizerTests, HasCorrectBounds)
{
    NonLinearNormalizer<double> nln(-1.0, 1.0, 1.0, 2000.0);

    EXPECT_NEAR(nln.normalize(1.0), -1.0, 0.01);
    EXPECT_NEAR(nln.normalize(2000.0), 1.0, 0.01);
    EXPECT_NEAR(nln.denormalize(-1.0), 1.0, 0.01);
    EXPECT_NEAR(nln.denormalize(1.0), 2000.0, 0.01);
}

TEST(NonLinearNormalizerTests, IsInvertible)
{
    double olb = 23.321770;
    double oub = 49.453453;
    NonLinearNormalizer<double> nln(-1.0, 1.0, olb, oub);

    for (int i = 0; i < 1000; i++)
    {
        double x = (double(i) / 1000.0) * (oub - olb) + olb;
        //cout << nln.normalize(x) << endl;
        EXPECT_NEAR(nln.denormalize(nln.normalize(x)), x, 1.0);
    }
}

TEST(NormalizersTests, HasCorrectBounds)
{
    Normalizer<double> nln(-1.0, 1.0, 1.0, 2000.0);

    EXPECT_NEAR(nln.normalize(1.0), -1.0, 0.01);
    EXPECT_NEAR(nln.normalize(2000.0), 1.0, 0.01);
    EXPECT_NEAR(nln.denormalize(-1.0), 1.0, 0.01);
    EXPECT_NEAR(nln.denormalize(1.0), 2000.0, 0.01);
}

TEST(NormalizersTests, IsInvertible)
{
    double olb = 23.321770;
    double oub = 49.453453;
    Normalizer<double> nln(-1.0, 1.0, olb, oub);

    for (int i = 0; i < 1000; i++)
    {
        double x = (double(i) / 1000.0) * (oub - olb) + olb;
        //cout << nln.normalize(x) << endl;
        EXPECT_NEAR(nln.denormalize(nln.normalize(x)), x, 1.0);
    }
}
