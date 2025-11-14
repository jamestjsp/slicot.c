#include <gtest/gtest.h>
#include <cmath>
#include "slicot.h"

TEST(MB01QD, FullMatrix) {
    const int32_t m = 3, n = 3;
    double a[9] = {
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0
    };
    const double cfrom = 2.0, cto = 4.0;
    int32_t info;
    int32_t nrows_dummy[1];

    mb01qd('G', m, n, 0, 0, cfrom, cto, 0, nrows_dummy, a, m, &info);

    double expected[9] = {
        2.0, 4.0, 6.0,
        8.0, 10.0, 12.0,
        14.0, 16.0, 18.0
    };

    for (int i = 0; i < 9; i++) {
        EXPECT_NEAR(a[i], expected[i], 1e-14);
    }
}

TEST(MB01QD, LowerTriangular) {
    const int32_t m = 3, n = 3;
    double a[9] = {
        1.0, 0.0, 0.0,
        2.0, 3.0, 0.0,
        4.0, 5.0, 6.0
    };
    const double cfrom = 1.0, cto = 2.0;
    int32_t info;
    int32_t nrows_dummy[1];

    mb01qd('L', m, n, 0, 0, cfrom, cto, 0, nrows_dummy, a, m, &info);

    double expected[9] = {
        2.0, 0.0, 0.0,
        4.0, 6.0, 0.0,
        8.0, 10.0, 12.0
    };

    for (int i = 0; i < 9; i++) {
        EXPECT_NEAR(a[i], expected[i], 1e-14);
    }
}

TEST(MB01QD, UpperTriangular) {
    const int32_t m = 3, n = 3;
    double a[9] = {
        1.0, 2.0, 3.0,
        0.0, 4.0, 5.0,
        0.0, 0.0, 6.0
    };
    const double cfrom = 1.0, cto = 0.5;
    int32_t info;
    int32_t nrows_dummy[1];

    mb01qd('U', m, n, 0, 0, cfrom, cto, 0, nrows_dummy, a, m, &info);

    double expected[9] = {
        0.5, 1.0, 1.5,
        0.0, 2.0, 2.5,
        0.0, 0.0, 3.0
    };

    for (int i = 0; i < 9; i++) {
        EXPECT_NEAR(a[i], expected[i], 1e-14);
    }
}

TEST(MB01QD, UpperHessenberg) {
    const int32_t m = 3, n = 3;
    double a[9] = {
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        0.0, 7.0, 8.0
    };
    const double cfrom = 1.0, cto = 3.0;
    int32_t info;
    int32_t nrows_dummy[1];

    mb01qd('H', m, n, 0, 0, cfrom, cto, 0, nrows_dummy, a, m, &info);

    double expected[9] = {
        3.0, 6.0, 9.0,
        12.0, 15.0, 18.0,
        0.0, 21.0, 24.0
    };

    for (int i = 0; i < 9; i++) {
        EXPECT_NEAR(a[i], expected[i], 1e-14);
    }
}

TEST(MB01QD, EmptyMatrix) {
    const int32_t m = 0, n = 3;
    double a[1];
    const double cfrom = 2.0, cto = 4.0;
    int32_t info;
    int32_t nrows_dummy[1];

    mb01qd('G', m, n, 0, 0, cfrom, cto, 0, nrows_dummy, a, 1, &info);
}

TEST(MB01QD, BlockLowerTriangular) {
    const int32_t m = 4, n = 4;
    double a[16] = {
        1.0, 2.0, 0.0, 0.0,
        3.0, 4.0, 0.0, 0.0,
        5.0, 6.0, 7.0, 8.0,
        9.0, 10.0, 11.0, 12.0
    };
    const double cfrom = 1.0, cto = 2.0;
    int32_t info;
    int32_t nrows[2] = {2, 2};

    mb01qd('L', m, n, 0, 0, cfrom, cto, 2, nrows, a, m, &info);

    double expected[16] = {
        2.0, 4.0, 0.0, 0.0,
        6.0, 8.0, 0.0, 0.0,
        10.0, 12.0, 14.0, 16.0,
        18.0, 20.0, 22.0, 24.0
    };

    for (int i = 0; i < 16; i++) {
        EXPECT_NEAR(a[i], expected[i], 1e-14);
    }
}
