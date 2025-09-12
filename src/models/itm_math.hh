#ifndef ITM_MATH_HH
#define ITM_MATH_HH

#include <algorithm>
#include <cmath>
#include <memory>
#include <numbers>
#include <ranges>
#include <span>

namespace itm_math {

    // Function to calculate the standard normal complementary CDF
    constexpr auto normalCCDF(double value) -> double {
        return 0.5 * std::erfc(value / std::numbers::sqrt2);
    }

    constexpr auto qerf(const double& z) -> double {
        return normalCCDF(z);
    }

    // The inverse of qerf - the solution for x to q = Q(x). The rational approximation
    // is due to Hastings, Jr. (1995) and the maximum error should be 4.5x10^-4.
    constexpr auto qerfi(double q) -> double
    {
        constexpr double c0 = 2.515516698;
        constexpr double c1 = 0.802853;
        constexpr double c2 = 0.010328;
        constexpr double d1 = 1.432788;
        constexpr double d2 = 0.189269;
        constexpr double d3 = 0.001308;

        const double x = 0.5 - q;
        const double t = std::sqrt(-2.0 * std::log(std::max(0.5 - std::abs(x), 0.000001)));
        const double v = t - ((c2 * t + c1) * t + c0) / (((d3 * t + d2) * t + d1) * t + 1.0);

        if (x < 0.0) {
            return -v;
        }

        return v;
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

    // Supporting function for the height gain in the "three radii method" used
    // in the computation of diffractive attenuation, as described in equations (4.20) and
    // (6.2)-(6.7) of "The ITS Irregular Terrain Model, version 1.2.2: The Algorithm" with
    // inputs corresponding to the "x" and "K" parameters of these equations.
    constexpr auto fht(const double &x, const double &pk) -> double
    {
        double fhtv = 0.0;

        if (x < 200.0) {
            const double w = -std::log(pk);

            if (pk < 1.0e-5 || x * w * w * w > 5495.0) {
                fhtv = -117.0;

                if (x > 1.0) {
                    fhtv = 40.0 * std::log10(x) + fhtv;
                }
            } else {
                fhtv = 2.5e-5 * x * x / pk - 8.686 * w - 15.0;
            }
        }

        else {
            fhtv = 0.05751 * x - 10.0 * std::log10(x);

            if (x < 2000.0) {
                const double w = 0.0134 * x * std::exp(-0.005 * x);
                fhtv = (1.0 - w) * fhtv + w * (40.0 * std::log10(x) - 117.0);
            }
        }
        return fhtv;
    }

    // Returns the function F0(D) (Eqn 6.9 of "The ITS Irregular Terrain Model, version 1.2.2:
    // The Algorithm") used in the computation of tropospheric scatter attenuation,
    // with the input D in meters.
    // The inputs and expected answer are based on an original test for Longley-Rice between
    // for Crystal Palace (South London) to Mursley, England (See Stark, 1967).
    constexpr auto ahd(double td) -> double
    {
        int i = 2;
        constexpr std::array<double, 3> a = { 133.4, 104.6, 71.8 };
        constexpr std::array<double, 3> b = { 0.332e-3, 0.212e-3, 0.157e-3 };
        constexpr std::array<double, 3> c = { -4.343, -1.086, 2.171 };

        if (td <= 10e3) {
            i = 0;
        }

        else if (td <= 70e3) {
            i = 1;
        }

        return a[i] + b[i] * td + c[i] * std::log(td);
    }

    // Tests the empirical curve fitting used in the computation of the Vmd, sigma_T-, and
    // sigma_T+ for estimating time variability effects as a function of the climatic region,
    // as described in equations (5.5) through (5.7) of of "The ITS Irregular Terrain Model,
    // version 1.2.2: The Algorithm" and as captured in Figure 10.13 of NBS Technical Note 101.
    constexpr auto curve(double const &c1, double const &c2, double const &x1,
                         double const &x2, double const &x3, double const &de) -> double
    {
        /* return (c1+c2/(1.0+pow((de-x2)/x3,2.0)))*pow(de/x1,2.0)/(1.0+pow(de/x1,2.0)); */
        const double temp1 = (de - x2) / x3;
        double temp2 = de / x1;

        temp2 *= temp2;

        return (c1 + c2 / (1.0 + temp1 * temp1)) * temp2 / (1.0 + temp2);
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

    /* #include <vector>
    #include <numeric> // For std::accumulate

    struct LinearFitResult {
        double slope;
        double intercept;
    };

    LinearFitResult linearLeastSquares(std::span<double> x, std::span<double> y) {
        if (x.size() != y.size() || x.empty()) {
            // Handle error: unequal sizes or empty data
            return {0.0, 0.0};
        }

        const int n = x.size();
        const double sum_x = std::ranges::accumulate(x, 0.0);
        const double sum_y = std::ranges::accumulate(y, 0.0);
        const double sum_xy = std::transform_reduce(x.begin(), x.end(), y.begin(), 0.0);
        const double sum_x2 = std::transform_reduce(x.begin(), x.end(), x.begin(), 0.0);

        const double denominator = n * sum_x2 - sum_x * sum_x;
        if (denominator == 0) {
            // Handle error: perfect vertical line or all x values are the same
            return {0.0, sum_y / n}; // Return horizontal line at mean y
        }

        const double slope = (n * sum_xy - sum_x * sum_y) / denominator;
        const double intercept = (sum_y - slope * sum_x) / n;

        return {slope, intercept};
    } */


    // A linear least squares fit between x1 and x2, to the function described
    // by the array z.
    // Evaluates a least squares fit to an input function z (in the form of a terrain profile
    // having first element the number of profile samples, second element the spacing between
    // them, and third through end elements the profile data) between horizontal locations
    // x1 and x2.  Returns the interpolated heights at location 0 and the end of the
    // profile.
    // TODO: Make this return z0 and zn as a pair or binding
    constexpr void z1sq1(std::span<double> z, const double &x1, const double &x2, double &z0, double &zn)
    {
        /* Used only with ITM 1.2.2 */
        const double xn = z[0];
        double xa = static_cast<int>(std::fdim(x1 / z[1], 0.0)); // index in z for x1 (int held as double)
        double xb = xn - static_cast<int>(std::fdim(xn, x2 / z[1])); // index in z for x2 (int held as double)

        // Handle case where dest is before start by adding 1 to xa and subtracting 1 from xb (with bounds handling)
        if (xb <= xa) {
            xa = std::fdim(xa, 1.0);
            xb = xn - std::fdim(xn, xb + 1.0);
        }

        int ja = static_cast<int>(xa); // index in z for x1 (int)
        const int jb = static_cast<int>(xb); // index in z for x2 (int)
        xa = xb - xa;
        double x = -0.5 * xa; // for some reason, x starts at either 0 or slightly lower
        xb += x; // Shift xb by this amount

        constexpr int ELEV_OFFSET = 2;
        // Now that things are setup, initialize a and b for the form y = ax + b (TODO correct?)
        // Why are these initialized this way? Some sort of online calculations?
        double a = 0.5 * (z[ja + ELEV_OFFSET] + z[jb + ELEV_OFFSET]);
        double b = 0.5 * (z[ja + ELEV_OFFSET] - z[jb + ELEV_OFFSET]) * x;

        // n = jb - ja
        // for(int i = 2; i <= n ; ++i) {}
        // ja; ja <= jb-2; ++ja
        for (; ja <= jb-ELEV_OFFSET; ++ja) {
            x += 1.0;
            a += z[ja + ELEV_OFFSET];
            b = std::fma(z[ja + ELEV_OFFSET], x, b);
        }

        a /= xa;
        b = b * 12.0 / ((xa * xa + 2.0) * xa);

        z0 = a - b * xb;
        zn = a + b * (xn - xb);
    }

    constexpr void z1sq2(std::span<double>z, const double &x1, const double &x2, double &z0, double &zn)
    {
        /* corrected for use with ITWOM */

        const double xn = z[0];
        double xa = static_cast<int>(std::fdim(x1 / z[1], 0.0));
        double xb = xn - static_cast<int>(std::fdim(xn, x2 / z[1]));

        if (xb <= xa) {
            xa = std::fdim(xa, 1.0);
            xb = xn - std::fdim(xn, xb + 1.0);
        }

        const int jb = static_cast<int>(xb);
        xa = (2 * static_cast<int>((xb - xa) / 2))-1;
        double x = -0.5 * (xa + 1);
        xb += x;
        int ja = jb - 1 - static_cast<int>(xa);
        const int n = jb - ja;
        double a = (z[ja + 2] + z[jb + 2]);
        double b = (z[ja + 2] - z[jb + 2]) * x;
        double bn = 2 * (x * x);

        for (int i = 2; i <= n; ++i) {
            ++ja;
            x += 1.0;
            bn += (x * x);
            a += z[ja + 2];
            b += z[ja + 2] * x;
        }

        a /= (xa + 2);
        b = b / bn;
        z0 = a - (b * xb);
        zn = a + (b * (xn - xb));
    }

    // Use the terrain profile pfl1 to find delta h, interdecile range of elevations between
    // point x1 and point x2, as described in Section 48 by Hufford.
    auto d1thx(double pfl[], const double &x1, const double &x2) -> double
    {
        const int np = static_cast<int>(pfl[0]);
        double xa = x1 / pfl[1];
        double xb = x2 / pfl[1];

        if (xb - xa < 2.0) {	// exit out
            return 0.0;
        }

        const int ka = std::clamp(static_cast<int>(0.1 * (xb - xa + 8.0)), 4, 25);
        const int n = 10 * ka - 5;
        const int kb = n - ka + 1;
        const double sn = n - 1;
        auto s = std::make_unique<double[]>(n + 2);
        s[0] = sn;
        s[1] = 1.0;
        xb = (xb - xa) / sn;
        int k = static_cast<int>(xa + 1.0);
        xa -= static_cast<double>(k);

        for (int j = 0; j < n; j++) {
            while (xa > 0.0 && k < np) {
                xa -= 1.0;
                ++k;
            }

            s[j + 2] = pfl[k + 2] + (pfl[k + 2] - pfl[k + 1]) * xa;
            xa += xb;
        }

        itm_math::z1sq1(std::span<double>(s.get(), s[0]+2), 0.0, sn, xa, xb);
        xb = (xb - xa) / sn;

        for (int j = 0; j < n; j++) {
            s[j + 2] -= xa;
            xa += xb;
        }

        // This was qtile, which really just returns the (clamped) i'th element of the sorted range
        // So the sort was moved out here to reduce the number of calls to it
        std::ranges::sort(std::span(s.get() + 2, n-1), std::greater<>());
        double d1thxv = s[2 + std::clamp(ka - 1, 0, static_cast<int>(n-2))] - s[2 + std::clamp(kb - 1, 0, static_cast<int>(n-2))];
        d1thxv /= 1.0 - 0.8 * std::exp(-(x2 - x1) / 50.0e3);

        return d1thxv;
    }

    auto d1thx2(double pfl[], const double &x1, const double &x2) -> double
    {
        const int np = static_cast<int>(pfl[0]);
        double xa = x1 / pfl[1];
        double xb = x2 / pfl[1];
        double d1thx2v = 0.0;

        if (xb - xa < 2.0) {	// exit out
            return d1thx2v;
        }

        const int kmx = std::max(25, static_cast<int>(83350 / (pfl[1])));
        const int ka = std::clamp(static_cast<int>(0.1 * (xb - xa + 8.0)), 4, kmx);
        const int n = 10 * ka - 5;
        const int kb = n - ka + 1;
        const double sn = n - 1;
        auto s = std::make_unique<double[]>(n + 2);
        s[0] = sn;
        s[1] = 1.0;
        xb = (xb - xa) / sn;
        int k = static_cast<int>(xa + 1.0);
        double xc = xa - static_cast<double>(k);

        for (int j = 0; j < n; j++) {
            while (xc > 0.0 && k < np) {
                xc -= 1.0;
                ++k;
            }

            s[j + 2] = pfl[k + 2] + (pfl[k + 2] - pfl[k + 1]) * xc;
            xc = xc + xb;
        }

        itm_math::z1sq2(std::span<double>(s.get(), s[0]+2), 0.0, sn, xa, xb);
        xb = (xb - xa) / sn;

        for (int j = 0; j < n; j++) {
            s[j + 2] -= xa;
            xa = xa + xb;
        }


        std::ranges::sort(std::span(s.get() + 2, n-1), std::greater<>());
        d1thx2v = s[2 + std::clamp(ka - 1, 0, static_cast<int>(n-2))] - s[2 + std::clamp(kb - 1, 0, static_cast<int>(n-2))];
        d1thx2v /= 1.0 - 0.8 * std::exp(-(x2 - x1) / 50.0e3);
        return d1thx2v;
    }

}; // namespace itm_math

#endif
