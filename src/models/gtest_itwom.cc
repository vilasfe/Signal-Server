#include "gtest/gtest.h"

#include "itm_math.hh"

namespace {
    TEST(TestQERF, QERF) {
        EXPECT_NEAR(itm_math::qerf(-1), 0.8413, 1e-4);
        EXPECT_NEAR(itm_math::qerf(1), 0.1587, 1e-4);
        EXPECT_NEAR(itm_math::qerf(10), 0.0, 1e-4);
    }

    TEST(TestH0F, H0F) {
        // test via area prediction mode
        constexpr double r = 0.4727387221558643;
        constexpr double et = 2.5221451983881447;

        EXPECT_DOUBLE_EQ(itm_math::h0f(r, et), 34.28154717133048);
        EXPECT_NEAR(itm_math::h0f(2, -1), 9.33, 0.01);
        EXPECT_NEAR(itm_math::h0f(2, 6), 18.53, 0.01);
    }
}
