#include "gtest/gtest.h"

#include "itm_math.hh"

namespace {
    TEST(TestQERF, QERF) {
        EXPECT_NEAR(itm_math::qerf(-1), 0.8413, 1e-4);
        EXPECT_NEAR(itm_math::qerf(1), 0.1587, 1e-4);
        EXPECT_NEAR(itm_math::qerf(10), 0.0, 1e-4);
    }
}
