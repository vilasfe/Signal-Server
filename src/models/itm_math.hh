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

    // Returns the attenuation due to a single knife edge - the Fresnel integral (in decibels,
    // Eqn 4.21 of "The ITS Irregular Terrain Model, version 1.2.2: The Algorithm" – see also
    // Eqn 6.1) evaluated for nu equal to the square root of the input argument.
    constexpr auto aknfe(double v2) -> double
    {
        // Trap for 0 value
        if (v2 <= 0) {
            v2 = 0.00001;
        }

        if (v2 < 5.76) {
            return 6.02 + 9.11 * std::sqrt(v2) - 1.27 * v2;
        }

        return 12.953 + 10 * std::log10(v2);
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
