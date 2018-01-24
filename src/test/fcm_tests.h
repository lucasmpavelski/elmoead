#include "../surrogate/ml/fcm.h"
#include "../surrogate/aux.h"

#include "helpers.h"


TEST(FCMTests, ALL)
{
    FuzzyCMeans::matrix data(75, 2);
    data << -0.2, 0.2,
            -0.2, 0.1,
            -0.2, 0.0,
            -0.2,-0.1,
            -0.2,-0.2,

            -0.1, 0.2,
            -0.1, 0.1,
            -0.1, 0.0,
            -0.1,-0.1,
            -0.1,-0.2,

             0.0, 0.2,
             0.0, 0.1,
             0.0, 0.0,
             0.0,-0.1,
             0.0,-0.2,

             0.1, 0.2,
             0.1, 0.1,
             0.1, 0.0,
             0.1,-0.1,
             0.1,-0.2,

             0.2, 0.2,
             0.2, 0.1,
             0.2, 0.0,
             0.2,-0.1,
             0.2,-0.2,

             0.8, 1.2,
             0.8, 1.1,
             0.8, 1.0,
             0.8, 0.9,
             0.8, 0.8,

             0.9, 1.2,
             0.9, 1.1,
             0.9, 1.0,
             0.9, 0.9,
             0.9, 0.8,

             1.0, 1.2,
             1.0, 1.1,
             1.0, 1.0,
             1.0, 0.9,
             1.0, 0.8,

             1.1, 1.2,
             1.1, 1.1,
             1.1, 1.0,
             1.1, 0.9,
             1.1, 0.8,

             1.2, 1.2,
             1.2, 1.1,
             1.2, 1.0,
             1.2, 0.9,
             1.2, 0.8,

             1.8, 0.2,
             1.8, 0.1,
             1.8, 0.0,
             1.8,-0.1,
             1.8,-0.2,

             1.9, 0.2,
             1.9, 0.1,
             1.9, 0.0,
             1.9,-0.1,
             1.9,-0.2,

             2.0, 0.2,
             2.0, 0.1,
             2.0, 0.0,
             2.0,-0.1,
             2.0,-0.2,

             2.1, 0.2,
             2.1, 0.1,
             2.1, 0.0,
             2.1,-0.1,
             2.1,-0.2,

             2.2, 0.2,
             2.2, 0.1,
             2.2, 0.0,
             2.2,-0.1,
             2.2,-0.2;


    //cout << data << endl;

    FuzzyCMeans fcm(75, 2, 3, 27);

    fcm.cluster(data);

    FuzzyCMeans::vector x(2);
    x << 1.4, 1.4;


    std::vector<int> r(2);
    fcm.classify(x, r);

    cout << printSeq(r.begin(), r.end()) << endl;

   /* EXPECT_EQ(r(0, 0), 0);
    EXPECT_EQ(r(1, 0), 0);
    EXPECT_EQ(r(2, 0), 0);
    EXPECT_EQ(r(3, 0), 0);
    EXPECT_EQ(r(4, 0), 1);
    EXPECT_EQ(r(5, 0), 1);
    EXPECT_EQ(r(6, 0), 1);
    EXPECT_EQ(r(7, 0), 1);*/
}
