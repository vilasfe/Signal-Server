#ifndef ITM_MATH_HH
#define ITM_MATH_HH

#include <cmath>
#include <numbers>

namespace itm_math {

    // Function to calculate the standard normal complementary CDF
    constexpr auto normalCCDF(double value) -> double {
        return 0.5 * std::erfc(value / std::numbers::sqrt2);
    }

    constexpr auto qerf(const double& z) -> double {
        return normalCCDF(z);
    }

    // Routine for computing the H01 "frequency gain" function described in
    // Eqn (6.13) of "The ITS Irregular Terrain Model, version 1.2.2: The Algorithm"
    // and used in computing troposcatter attenuation.
    constexpr auto h0f(double r, double et) -> double
    {
        constexpr std::array<double, 5> a = { 25.0, 80.0, 177.0, 395.0, 705.0 };
        constexpr std::array<double, 5> b = { 24.0, 45.0, 68.0, 80.0, 105.0 };
        double q = 0.0;
        int it = static_cast<int>(et);

        if (it <= 0) {
            it = 1;
        }

        else if (it >= 5) {
            it = 5;
        }

        else {
            q = et - it;
        }

        /* x=pow(1.0/r,2.0); */

        const double temp = 1.0 / r;
        const double x = temp * temp;

        double h0fv = 4.343 * std::log((a[it - 1] * x + b[it - 1]) * x + 1.0);

        if (q != 0.0) {
            h0fv =
            (1.0 - q) * h0fv + q * 4.343 * std::log((a[it] * x + b[it]) * x + 1.0);
        }

        return h0fv;
    }

};

#endif
