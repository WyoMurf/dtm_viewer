/*
 * longley_rice.c -- faithful C port of the NTIA/ITS Irregular Terrain Model
 * (Longley-Rice) reference implementation, point-to-point TLS mode only.
 * See include/longley_rice.h for the public entry point and usage notes.
 *
 * Ported from the official C++ source at https://github.com/NTIA/itm
 * (branch master, fetched 2026-09-22), a US federal government work
 * released to the public domain (Title 15 USC 105) with explicit
 * permission to modify/redistribute (see that repo's LICENSE.md). Each of
 * the ~19 functions below corresponds 1:1 to a same-named function in one
 * of that repo's src/ .cpp files -- this is a mechanical translation
 * (C++ references -> C pointers, std::complex<double> -> C99 double
 * complex, std::vector-based order-statistics -> a fixed array + qsort),
 * not a reimplementation. Every constant, formula and branch is kept
 * exactly as the original, including anything that looks like an odd
 * magic number or a strangely specific threshold -- those come from the
 * underlying NTIA technical notes (referenced in the comments below, e.g.
 * "[TN101, Eqn ...]" or "[ERL 79-ITS 67, Eqn ...]") and are not something
 * this port should "clean up".
 *
 * Left out of this port (not needed for point-to-point work with a real
 * digitized terrain profile, per the task this module was written for):
 * area-prediction mode (ITM_AREA_*, InitializeArea) and the
 * confidence/reliability (CR) variability variant (ITM_P2P_CR) -- only
 * ITM_P2P_TLS is exposed.
 */

#include "include/longley_rice.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>

/////////////////////////////////////////////////////////////////////////
// Constants and macros carried over from the original include/itm.h,
// include/Enums.h -- kept private to this translation unit since only the
// subset a caller actually needs (climate/pol/mdvar/error/warning values)
// is re-declared in the public header.

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define DIM(x, y) (((x) > (y)) ? ((x) - (y)) : (0))

#define PI                                      3.1415926535897932384
#define SQRT2                                   sqrt(2)
#define a_0__meter                              6370e3
#define a_9000__meter                           9000e3
#define THIRD                                   (1.0 / 3.0)

#define MODE__P2P                               0
#define MODE__AREA                              1

// List of valid modes of propagation (LongleyRice's *propmode output --
// this port discards it, since the public API doesn't expose
// IntermediateValues, but LongleyRice still needs somewhere to write it)
#define MODE__NOT_SET                           0
#define MODE__LINE_OF_SIGHT                     1
#define MODE__DIFFRACTION                       2
#define MODE__TROPOSCATTER                      3

// Variability.cpp compares against these under their original bare names
// (not the MDVAR__ prefixed ones) -- same numeric values, see itm.h's
// Enums.h, which is why aliasing is exact and lossless.
#define SINGLE_MESSAGE_MODE                     MDVAR__SINGLE_MESSAGE_MODE
#define ACCIDENTAL_MODE                         MDVAR__ACCIDENTAL_MODE
#define MOBILE_MODE                             MDVAR__MOBILE_MODE
#define BROADCAST_MODE                          MDVAR__BROADCAST_MODE

/* Struct of intermediate values LongleyRice/ITM_P2P_TLS_Ex compute along
 * the way to A__db. The original exposes this to callers (ITM_P2P_TLS_Ex);
 * this port keeps it purely internal since the public API here (see
 * longley_rice.h) only exposes the plain ITM_P2P_TLS signature. */
typedef struct {
    double theta_hzn[2];
    double d_hzn__meter[2];
    double h_e__meter[2];
    double N_s;
    double delta_h__meter;
    double A_ref__db;
    double A_fs__db;
    double d__km;
    int mode;
} IntermediateValues;

/////////////////////////////////////////////////////////////////////////
// Forward declarations (mirrors the DLLEXPORT declarations in itm.h)

static double ComputeDeltaH(const double pfl[], double d_start__meter, double d_end__meter);
static double DiffractionLoss(double d__meter, const double d_hzn__meter[2], const double h_e__meter[2], double complex Z_g,
    double a_e__meter, double delta_h__meter, const double h__meter[2], int mode, double theta_los, double d_sML__meter, double f__mhz);
static double FFunction(double td);
static void FindHorizons(const double pfl[], double a_e__meter, const double h__meter[2], double theta_hzn[2], double d_hzn__meter[2]);
static double FreeSpaceLoss(double d__meter, double f__mhz);
static double FresnelIntegral(double v2);
static double H0Curve(int j, double r);
static double H0Function(double r, double eta_s);
static double HeightFunction(double x__km, double K);
static void InitializePointToPoint(double f__mhz, double h_sys__meter, double N_0, int pol, double epsilon,
    double sigma, double complex *Z_g, double *gamma_e, double *N_s);
static double InverseComplementaryCumulativeDistributionFunction(double q);
static double KnifeEdgeDiffraction(double d__meter, double f__mhz, double a_e__meter, double theta_los, const double d_hzn__meter[2]);
static void LinearLeastSquaresFit(const double pfl[], double d_start, double d_end, double *fit_y1, double *fit_y2);
static double LineOfSightLoss(double d__meter, const double h_e__meter[2], double complex Z_g, double delta_h__meter,
    double M_d, double A_d0, double d_sML__meter, double f__mhz);
static int LongleyRice(const double theta_hzn[2], double f__mhz, double complex Z_g, const double d_hzn__meter[2], const double h_e__meter[2],
    double gamma_e, double N_s, double delta_h__meter, const double h__meter[2], double d__meter, int mode, double *A_ref__db,
    long *warnings, int *propmode);
static void QuickPfl(const double pfl[], double gamma_e, const double h__meter[2], double theta_hzn[2], double d_hzn__meter[2],
    double h_e__meter[2], double *delta_h__meter, double *d__meter);
static double SigmaHFunction(double delta_h__meter);
static double SmoothEarthDiffraction(double d__meter, double f__mhz, double a_e__meter, double theta_los,
    const double d_hzn__meter[2], const double h_e__meter[2], double complex Z_g);
static double TerrainRoughness(double d__meter, double delta_h__meter);
static double TroposcatterLoss(double d__meter, const double theta_hzn[2], const double d_hzn__meter[2], const double h_e__meter[2],
    double a_e__meter, double N_s, double f__mhz, double theta_los, double *h0);
static int ValidateInputs(double h_tx__meter, double h_rx__meter, int climate, double time,
    double location, double situation, double N_0, double f__mhz, int pol,
    double epsilon, double sigma, int mdvar, long *warnings);
static double Variability(double time, double location, double situation, const double h_e__meter[2], double delta_h__meter,
    double f__mhz, double d__meter, double A_ref__db, int climate, int mdvar, long *warnings);
static int ITM_P2P_TLS_Ex(const double h_tx__meter, const double h_rx__meter, const double pfl[], const int climate, const double N_0,
    const double f__mhz, const int pol, const double epsilon, const double sigma, const int mdvar, const double time, const double location,
    const double situation, double *A__db, long *warnings, IntermediateValues *interValues);

/////////////////////////////////////////////////////////////////////////

/* Descending-order comparator for qsort(), used by ComputeDeltaH() below
 * in place of the original's std::nth_element(..., std::greater<double>())
 * calls. The original does two partial-order (O(n)) selections rather
 * than a full O(n log n) sort, but for a fixed-size array of at most 245
 * doubles (see ComputeDeltaH's `s[247]`) the two are indistinguishable in
 * practice, and a full descending sort followed by indexing [p10-1] and
 * [p90] produces bit-identical q10/q90 values to what two nth_element
 * calls with the same comparator would leave there -- nth_element is
 * defined by the *order statistic* it leaves at that index, which is the
 * same value a full sort would put there regardless of algorithm. */
static int CompareDoubleDescending(const void *a, const void *b) {
    double da = *(const double *)a, db = *(const double *)b;
    if (da < db) return 1;
    if (da > db) return -1;
    return 0;
}

/*=============================================================================
 |
 |  Description:  Compute the terrain irregularity parameter, delta_h
 |
 |        Input:  pfl[]          - Terrain data
 |                d_start__meter - Distance into the terrain profile to start
 |                                 considering data, in meters
 |                d_end__meter   - Distance into the terrain profile to end
 |                                 considering data, in meters
 |
 |      Outputs:  [None]
 |
 |      Returns:  delta_h__meter - Terrain irregularity parameter, in meters
 |
 *===========================================================================*/
static double ComputeDeltaH(const double pfl[], double d_start__meter, double d_end__meter) {
    double s[247] = { 0 };                      // Temp pfl data array

    const int np = (int)pfl[0];
    double x_start = d_start__meter / pfl[1];   // index to start considering terrain points
    double x_end = d_end__meter / pfl[1];       // index to stop considering terrain points

    // if there are less than 2 terrain points, return delta_h = 0
    if (x_end - x_start < 2.0)
        return 0;

    int p10 = (int)(0.1 * (x_end - x_start + 8.0));
    p10 = MIN(MAX(4, p10), 25);                 // 10% index

    const int n = 10 * p10 - 5;
    int p90 = n - p10;                         // 90% index

    double np_s = n - 1;
    s[0] = np_s;
    s[1] = 1.0;

    x_end = (x_end - x_start) / np_s;
    int i = (int)x_start;
    x_start -= (float)(i + 1.0);

    for (int j = 0; j < n; j++) {
        while (x_start > 0.0 && (i + 1) < np) {
            x_start--;
            i++;
        }

        s[j + 2] = pfl[i + 3] + (pfl[i + 3] - pfl[i + 2]) * x_start;

        x_start += x_end;
    }

    double fit_y1, fit_y2;
    LinearLeastSquaresFit(s, 0.0, np_s, &fit_y1, &fit_y2);

    fit_y2 = (fit_y2 - fit_y1) / np_s;

    double diffs[245];

    // compute the difference between fitted line and actual data
    for (int j = 0; j < n; j++) {
        diffs[j] = s[j + 2] - fit_y1;

        fit_y1 += fit_y2;
    }

    qsort(diffs, (size_t)n, sizeof(double), CompareDoubleDescending);
    const double q10 = diffs[p10 - 1];
    const double q90 = diffs[p90];

    const double delta_h_d__meter = q10 - q90;

    // [ERL 79-ITS 67, Eqn 3], inverted
    const double delta_h__meter = delta_h_d__meter / (1.0 - 0.8 * exp(-(d_end__meter - d_start__meter) / 50e3));

    return delta_h__meter;
}

/*=============================================================================
 |
 |  Description:  Compute the diffraction loss at a specified distance
 |
 |        Input:  d__meter          - Path distance, in meters
 |                d_hzn__meter[2]   - Horizon distances, in meters
 |                h_e__meter[2]     - Effective terminal heights, in meters
 |                Z_g            - Complex ground impedance
 |                a_e__meter     - Effective earth radius, in meters
 |                delta_h__meter - Terrain irregularity parameter, in meters
 |                h__meter[2]       - Terminal heights, in meters
 |                mode           - Area or Point-to-Point mode flag
 |                theta_los      - Angular distance of line-of-sight region
 |                d_sML__meter   - Maximum line-of-sight distance for
 |                                 a smooth earth, in meters
 |                f__mhz         - Frequency, in MHz
 |
 |      Outputs:  [None]
 |
 |      Returns:  A_d__db        - Diffraction loss, in dB
 |
 *===========================================================================*/
static double DiffractionLoss(double d__meter, const double d_hzn__meter[2], const double h_e__meter[2], double complex Z_g,
    double a_e__meter, double delta_h__meter, const double h__meter[2], int mode, double theta_los, double d_sML__meter, double f__mhz) {
    const double A_k__db = KnifeEdgeDiffraction(d__meter, f__mhz, a_e__meter, theta_los, d_hzn__meter);

    const double A_se__db = SmoothEarthDiffraction(d__meter, f__mhz, a_e__meter, theta_los, d_hzn__meter, h_e__meter, Z_g);

    //////////////////
    // Terrain clutter

    // Terrain roughness term, using d_sML__meter, per [ERL 79-ITS 67, page 3-13]
    const double delta_h_dsML__meter = TerrainRoughness(d_sML__meter, delta_h__meter);

    const double sigma_h_d__meter = SigmaHFunction(delta_h_dsML__meter);

    // Clutter factor
    // [ERL 79-ITS 67, Eqn 3.38c]
    const double A_fo__db = MIN(15.0, 5 * log10(1.0 + 1e-5 * h__meter[0] * h__meter[1] * f__mhz * sigma_h_d__meter));

    //////////////////////////////
    // Combined diffraction losses

    // compute the weighting factor in the following calculations

    const double delta_h_d__meter = TerrainRoughness(d__meter, delta_h__meter);

    double q = h__meter[0] * h__meter[1];
    const double qk = h_e__meter[0] * h_e__meter[1] - q;

    // For low antennas with known path parameters, C ~= 10 [ERL 79-ITS 67, page 3-8]
    if (mode == MODE__P2P)
        q += 10.0;

    const double term1 = sqrt(1.0 + qk / q);                              // square root term in [ERL 79-ITS 67, Eqn 3.23]

    const double d_ML__meter = d_hzn__meter[0] + d_hzn__meter[1];         // Maximum line-of-sight distance for actual path
    q = (term1 + (-theta_los * a_e__meter + d_ML__meter) / d__meter) * MIN(delta_h_d__meter * f__mhz / 47.7, 6283.2);

    // weighting factor [ERL 17-ITS 67, Eqn 3.23]
    const double w = 25.1 / (25.1 + sqrt(q));

    const double A_d__db = w * A_se__db + (1.0 - w) * A_k__db + A_fo__db;

    return A_d__db;
}

/*=============================================================================
 |
 |  Description:  Compute the radio horizon's of the terminals
 |
 |        Input:  pfl[]             - Terrain data
 |                a_e__meter        - Effective earth radius, in meters
 |                h__meter[2]       - Terminal structural heights, in meters
 |
 |      Outputs:  theta_hzn[2]      - Terminal radio horizon angle, in radians
 |                d_hzn__meter[2]   - Terminal radio horizon distance, in meters
 |
 |      Returns:  [None]
 |
 *===========================================================================*/
static void FindHorizons(const double pfl[], double a_e__meter, const double h__meter[2], double theta_hzn[2], double d_hzn__meter[2]) {
    const int np = (int)pfl[0];
    const double xi = pfl[1];

    const double d__meter = pfl[0] * pfl[1];

    // compute radials (ignore radius of earth since it cancels out in the later math)
    const double z_tx__meter = pfl[2] + h__meter[0];
    const double z_rx__meter = pfl[np + 2] + h__meter[1];

    // set the terminal horizon angles as if the terminals are line-of-sight
    // [TN101, Eq 6.15]
    theta_hzn[0] = (z_rx__meter - z_tx__meter) / d__meter - d__meter / (2 * a_e__meter);
    theta_hzn[1] = -(z_rx__meter - z_tx__meter) / d__meter - d__meter / (2 * a_e__meter);

    d_hzn__meter[0] = d__meter;
    d_hzn__meter[1] = d__meter;

    double d_tx__meter = 0.0;
    double d_rx__meter = d__meter;

    double theta_tx, theta_rx;

    for (int i = 1; i < np; i++) {
        d_tx__meter = d_tx__meter + xi;
        d_rx__meter = d_rx__meter - xi;

        theta_tx = (pfl[i + 2] - z_tx__meter) / d_tx__meter - d_tx__meter / (2 * a_e__meter);
        theta_rx = -(z_rx__meter - pfl[i + 2]) / d_rx__meter - d_rx__meter / (2 * a_e__meter);

        if (theta_tx > theta_hzn[0]) {
            theta_hzn[0] = theta_tx;
            d_hzn__meter[0] = d_tx__meter;
        }

        if (theta_rx > theta_hzn[1]) {
            theta_hzn[1] = theta_rx;
            d_hzn__meter[1] = d_rx__meter;
        }
    }
}

/*=============================================================================
 |
 |  Description:  Free space basic transmission loss equation
 |
 |        Input:  d__meter       - Path distance, in meters
 |                f__mhz         - Frequency, in MHz
 |
 |      Outputs:  [None]
 |
 |      Returns:  A_fs__db       - Free space basic transmission loss, in dB
 |
 *===========================================================================*/
static double FreeSpaceLoss(double d__meter, double f__mhz) {
    return 32.45 + 20.0 * log10(f__mhz) + 20.0 * log10(d__meter / 1000.0);
}

/*=============================================================================
 |
 |  Description:  Approximate to ideal knife edge diffraction loss
 |
 |        Input:  v2             - v^2 parameter
 |
 |      Outputs:  [None]
 |
 |      Returns:  A(v, 0)        - Loss, in dB
 |
 *===========================================================================*/
static double FresnelIntegral(double v2) {
    // Note: v2  is v^2, so 5.76 is actually comparing v to 2.4

    if (v2 < 5.76)
        return 6.02 + 9.11 * sqrt(v2) - 1.27 * v2;      // [TN101v2, Eqn III.24b] and [ERL 79-ITS 67, Eqn 3.27a & 3.27b]
    else
        return 12.953 + 10 * log10(v2);                 // [TN101v2, Eqn III.24c] and [ERL 79-ITS 67, Eqn 3.27a & 3.27b]
}

/*=============================================================================
 |
 |  Description:  Curve fit helper function to approximate H_0()
 |
 |        Input:  j              - eta_s curve
 |                r              - Input parameter r_1,2
 |
 |      Outputs:  [None]
 |
 |      Returns: H_01(r, j)      - in dB
 |
 *===========================================================================*/
static double H0Curve(int j, double r) {
    // values from [Algorithm, 6.13]
    const double a[] = { 25.0, 80.0, 177.0, 395.0, 705.0 };
    const double b[] = { 24.0, 45.0, 68.0, 80.0, 105.0 };

    return 10 * log10(1 + a[j] * pow(1 / r, 4) + b[j] * pow(1.0 / r, 2));    // related to TN101v2, Eqn III.49, but from [Algorithm, 6.13]
}

/*=============================================================================
 |
 |  Description:  Troposcatter frequency gain function, H_0(), from
 |                [TN101v1, Ch 9.2]
 |
 |        Input:  r              - Input parameter r_1,2
 |                eta_s          - Parameter eta_s
 |
 |      Outputs:  [None]
 |
 |      Returns:  H_0()          - in dB
 |
 *===========================================================================*/
static double H0Function(double r, double eta_s) {
    eta_s = MIN(MAX(eta_s, 1), 5);  // range 1 <= eta_s <= 5

    const int i = (int)eta_s;       // integer part of eta_s
    const double q = eta_s - i;     // decimal part of eta_s

    double result = H0Curve(i - 1, r);

    if (q != 0.0)                   // interpolate with next curve, if needed
        result = (1.0 - q) * result + q * H0Curve(i, r);

    return result;
}

/*=============================================================================
 |
 |  Description:  Initialize parameters for point-to-point mode
 |
 |        Input:  f__mhz            - Frequency, in MHz
 |                h_sys__meter      - Average height of the path above
 |                                    mean sea level, in meters
 |                N_0               - Refractivity, in N-Units
 |                pol               - Polarization
 |                                      + 0 : POLARIZATION__HORIZONTAL
 |                                      + 1 : POLARIZATION__VERTICAL
 |                epsilon           - Relative permittivity
 |                sigma             - Conductivity
 |
 |      Outputs:  Z_g               - Complex ground impedance
 |                gamma_e           - Curvature of the effective earth
 |                N_s               - Surface refractivity, in N-Units
 |
 |      Returns:  [None]
 |
 *===========================================================================*/
static void InitializePointToPoint(double f__mhz, double h_sys__meter, double N_0, int pol,
    double epsilon, double sigma, double complex *Z_g, double *gamma_e, double *N_s) {
    // gamma_a is the curvature of the actual earth, ~1 / 6370 km
    const double gamma_a = 157e-9;

    // scale the refractivity based on the elevation above mean sea level
    if (h_sys__meter == 0.0)
        *N_s = N_0;
    else
        *N_s = N_0 * exp(-h_sys__meter / 9460.0);               // [TN101, Eq 4.3]

    // gamma_e is the curvature of the effective earth
    *gamma_e = gamma_a * (1.0 - 0.04665 * exp(*N_s / 179.3));   // [TN101, Eq 4.4], reworked

    // complex relative permittivity
    const double complex ep_r = epsilon + (18000 * sigma / f__mhz) * I;

    *Z_g = csqrt(ep_r - 1.0);                        // ground impedance (horizontal polarization)

    if (pol == POLARIZATION__VERTICAL)              // adjust for vertical polarization
        *Z_g = *Z_g / ep_r;
}

/*=============================================================================
 |
 |  Description:  This function computes the inverse complementary
 |                cumulative distribution function approximation as
 |                described in Formula 26.2.23 in Abramowitz & Stegun.
 |                This approximation has an error of
 |                abs(epsilon(p)) < 4.5e-4
 |
 |        Input:  q              - Quantile, 0.0 < q < 1.0
 |
 |      Outputs:  [None]
 |
 |      Returns:  Q_q            - Q(q)^-1
 |
 *===========================================================================*/
static double InverseComplementaryCumulativeDistributionFunction(double q) {
    const double C_0 = 2.515516;
    const double C_1 = 0.802853;
    const double C_2 = 0.010328;
    const double D_1 = 1.432788;
    const double D_2 = 0.189269;
    const double D_3 = 0.001308;

    double x = q;
    if (q > 0.5)
        x = 1.0 - x;

    const double T_x = sqrt(-2.0 * log(x));

    const double zeta_x = ((C_2 * T_x + C_1) * T_x + C_0) / (((D_3 * T_x + D_2) * T_x + D_1) * T_x + 1.0);

    double Q_q = T_x - zeta_x;

    if (q > 0.5)
        Q_q = -Q_q;

    return Q_q;
}

/*=============================================================================
 |
 |  Description:  Compute the knife-edge diffraction loss
 |
 |        Input:  d__meter          - Distance of interest, in meters
 |                f__mhz            - Frequency, in MHz
 |                a_e__meter        - Effective earth radius, in meters
 |                theta_los         - Angular distance of line-of-sight region
 |                d_hzn__meter[2]   - Horizon distances, in meters
 |
 |      Outputs:  [None]
 |
 |      Returns:  A_k__db        - Knife-edge diffraction loss, in dB
 |
 *===========================================================================*/
static double KnifeEdgeDiffraction(double d__meter, double f__mhz, double a_e__meter, double theta_los, const double d_hzn__meter[2]) {
    const double d_ML__meter = d_hzn__meter[0] + d_hzn__meter[1];        // Maximum line-of-sight distance for actual path
    const double theta_nlos = d__meter / a_e__meter - theta_los;         // Angular distance of diffraction region [Algorithm, Eqn 4.12]

    const double d_nlos__meter = d__meter - d_ML__meter;                 // Diffraction distance, in meters

    // 1 / (4 pi) = 0.0795775
    // [TN101, Eqn I.7]
    const double v_1 = 0.0795775 * (f__mhz / 47.7) * pow(theta_nlos, 2) * d_hzn__meter[0] * d_nlos__meter / (d_nlos__meter + d_hzn__meter[0]);
    const double v_2 = 0.0795775 * (f__mhz / 47.7) * pow(theta_nlos, 2) * d_hzn__meter[1] * d_nlos__meter / (d_nlos__meter + d_hzn__meter[1]);

    const double A_k__db = FresnelIntegral(v_1) + FresnelIntegral(v_2);  // [TN101, Eqn I.1]

    return A_k__db;
}

/*=============================================================================
 |
 |  Description:  Perform a linear least squares fit to the terrain data
 |
 |        Input:  pfl[2]         - Input data array, in pfl format
 |                d_start        - Start distance
 |                d_end          - End distance
 |
 |      Outputs:  fit_y1         - Fitted y1 value
 |                fit_y2         - Fitted y2 value
 |
 |      Returns:  [None]
 |
 *===========================================================================*/
static void LinearLeastSquaresFit(const double pfl[], double d_start, double d_end, double *fit_y1, double *fit_y2) {
    const int np = (int)pfl[0];

    int i_start = (int)fdim(d_start / pfl[1], 0.0);
    int i_end = np - (int)fdim(np, d_end / pfl[1]);

    if (i_end <= i_start) {
        i_start = (int)fdim(i_start, 1.0);
        i_end = np - (int)fdim(np, i_end + 1.0);
    }

    const double x_length = i_end - i_start;

    double mid_shifted_index = -0.5 * x_length;
    const double mid_shifted_end = i_end + mid_shifted_index;

    double sum_y = 0.5 * (pfl[i_start + 2] + pfl[i_end + 2]);
    double scaled_sum_y = 0.5 * (pfl[i_start + 2] - pfl[i_end + 2]) * mid_shifted_index;

    for (int i = 2; i <= x_length; i++) {
        i_start++;
        mid_shifted_index++;

        sum_y += pfl[i_start + 2];
        scaled_sum_y += pfl[i_start + 2] * mid_shifted_index;
    }

    sum_y = sum_y / x_length;
    scaled_sum_y = scaled_sum_y * 12.0 / ((x_length * x_length + 2.0) * x_length);

    *fit_y1 = sum_y - scaled_sum_y * mid_shifted_end;
    *fit_y2 = sum_y + scaled_sum_y * (np - mid_shifted_end);
}

/*=============================================================================
 |
 |  Description:  Compute the loss in the line-of-sight region
 |
 |        Input:  d__meter          - Path distance, in meters
 |                h_e__meter[2]     - Terminal effective heights, in meters
 |                Z_g               - Complex surface transfer impedance
 |                delta_h__meter    - Terrain irregularity parameter
 |                M_d               - Diffraction slope
 |                A_d0              - Diffraction intercept
 |                d_sML__meter      - Maximum line-of-sight distance for
 |                                    a smooth earth, in meters
 |                f__mhz            - Frequency, in MHz
 |
 |      Outputs:  [None]
 |
 |      Returns:  A_los__db         - Loss, in dB
 |
 *===========================================================================*/
static double LineOfSightLoss(double d__meter, const double h_e__meter[2], double complex Z_g, double delta_h__meter,
    double M_d, double A_d0, double d_sML__meter, double f__mhz) {
    const double delta_h_d__meter = TerrainRoughness(d__meter, delta_h__meter);

    const double sigma_h_d__meter = SigmaHFunction(delta_h_d__meter);

    // wavenumber, k
    const double wn = f__mhz / 47.7;

    // [Algorithm, Eqn 4.46]
    const double sin_psi = (h_e__meter[0] + h_e__meter[1]) / sqrt(pow(d__meter, 2) + pow(h_e__meter[0] + h_e__meter[1], 2));

    // [Algorithm, Eqn 4.47]
    double complex R_e = (sin_psi - Z_g) / (sin_psi + Z_g) * exp(-MIN(10.0, wn * sigma_h_d__meter * sin_psi));

    // q = Magnitude of R_e', [Algorithm, Eqn 4.48]
    const double q = pow(creal(R_e), 2) + pow(cimag(R_e), 2);
    if (q < 0.25 || q < sin_psi)
        R_e = R_e * sqrt(sin_psi / q);

    // phase difference between rays, [Algorithm, Eqn 4.49]
    double delta_phi = wn * 2.0 * h_e__meter[0] * h_e__meter[1] / d__meter;

    // [Algorithm, Eqn 4.50]
    if (delta_phi > PI / 2.0)
        delta_phi = PI - pow(PI / 2.0, 2) / delta_phi;

    // Two-ray attenuation
    const double complex rr = (cos(delta_phi) - sin(delta_phi) * I) + R_e;
    const double A_t__db = -10 * log10(pow(creal(rr), 2) + pow(cimag(rr), 2));

    // Extended diffraction attenuation
    const double A_d__db = M_d * d__meter + A_d0;

    // weighting factor
    const double w = 1 / (1 + f__mhz * delta_h__meter / MAX(10e3, d_sML__meter));

    const double A_los__db = w * A_t__db + (1 - w) * A_d__db;

    return A_los__db;
}

/*=============================================================================
 |
 |  Description:  Compute the reference attenuation, using the
 |                Longley-Rice method
 |
 |        Input:  theta_hzn[2]      - Terminal horizon angles
 |                f__mhz            - Frequency, in MHz
 |                Z_g               - Complex surface transfer impedance
 |                d_hzn__meter[2]   - Terminal horizon distances, in meters
 |                h_e__meter[2]     - Effective terminal heights, in meters
 |                gamma_e           - Curvature of the effective earth
 |                N_s               - Surface refractivity, in N-Units
 |                delta_h__meter    - Terrain irregularity parameter
 |                h__meter[2]       - Terminal structural heights, in meters
 |                d__meter          - Path distance, in meters
 |                mode              - Mode of operation (P2P or Area)
 |
 |      Outputs:  A_ref__db         - Reference attenuation, in dB
 |                warnings          - Warning flags
 |                propmode          - Mode of propagation value
 |
 |      Returns:  error             - Error code
 |
 *===========================================================================*/
static int LongleyRice(const double theta_hzn[2], double f__mhz, double complex Z_g, const double d_hzn__meter[2],
    const double h_e__meter[2], double gamma_e, double N_s, double delta_h__meter, const double h__meter[2],
    double d__meter, int mode, double *A_ref__db, long *warnings, int *propmode) {
    // effective earth radius
    const double a_e__meter = 1 / gamma_e;

    double d_hzn_s__meter[2];
    // Terrestrial smooth earth horizon distance approximation
    for (int i = 0; i < 2; i++)
        d_hzn_s__meter[i] = sqrt(2.0 * h_e__meter[i] * a_e__meter);

    // Maximum line-of-sight distance for smooth earth
    const double d_sML__meter = d_hzn_s__meter[0] + d_hzn_s__meter[1];

    // Maximum line-of-sight distance for actual path
    const double d_ML__meter = d_hzn__meter[0] + d_hzn__meter[1];

    // Angular distance of line-of-sight region
    const double theta_los = -MAX(theta_hzn[0] + theta_hzn[1], -d_ML__meter / a_e__meter);

    // Check validity of small angle approximation
    if (fabs(theta_hzn[0]) > 200e-3)
        *warnings |= WARN__TX_HORIZON_ANGLE;
    if (fabs(theta_hzn[1]) > 200e-3)
        *warnings |= WARN__RX_HORIZON_ANGLE;

    // Checks that the actual horizon distance can't be less than 1/10 of the smooth earth horizon distance
    if (d_hzn__meter[0] < 0.1 * d_hzn_s__meter[0])
        *warnings |= WARN__TX_HORIZON_DISTANCE_1;
    if (d_hzn__meter[1] < 0.1 * d_hzn_s__meter[1])
        *warnings |= WARN__RX_HORIZON_DISTANCE_1;

    // Checks that the actual horizon distance can't be greater than 3 times the smooth earth horizon distance
    if (d_hzn__meter[0] > 3.0 * d_hzn_s__meter[0])
        *warnings |= WARN__TX_HORIZON_DISTANCE_2;
    if (d_hzn__meter[1] > 3.0 * d_hzn_s__meter[1])
        *warnings |= WARN__RX_HORIZON_DISTANCE_2;

    // Check the surface refractivity
    if (N_s < 150)
        return ERROR__SURFACE_REFRACTIVITY_SMALL;
    if (N_s > 400)
        return ERROR__SURFACE_REFRACTIVITY_LARGE;
    if (N_s < 250) // 150 <= N_s < 250
        *warnings |= WARN__SURFACE_REFRACTIVITY;

    // Check effective earth size
    if (a_e__meter < 4000000 || a_e__meter > 13333333)
        return ERROR__EFFECTIVE_EARTH;

    // Check ground impedance
    if (creal(Z_g) <= fabs(cimag(Z_g)))
        return ERROR__GROUND_IMPEDANCE;

    // Select two distances far in the diffraction region
    const double d_3__meter = MAX(d_sML__meter, d_ML__meter + 5.0 * pow(pow(a_e__meter, 2) / f__mhz, 1.0 / 3.0));
    const double d_4__meter = d_3__meter + 10.0 * pow(pow(a_e__meter, 2) / f__mhz, 1.0 / 3.0);

    // Compute the diffraction loss at the two distances
    const double A_3__db = DiffractionLoss(d_3__meter, d_hzn__meter, h_e__meter, Z_g, a_e__meter, delta_h__meter, h__meter, mode, theta_los, d_sML__meter, f__mhz);
    const double A_4__db = DiffractionLoss(d_4__meter, d_hzn__meter, h_e__meter, Z_g, a_e__meter, delta_h__meter, h__meter, mode, theta_los, d_sML__meter, f__mhz);

    // Compute the slope and intercept of the diffraction line
    const double M_d = (A_4__db - A_3__db) / (d_4__meter - d_3__meter);
    const double A_d0__db = A_3__db - M_d * d_3__meter;

    const double d_min__meter = fabs(h_e__meter[0] - h_e__meter[1]) / 200e-3;

    if (d__meter < d_min__meter)
        *warnings |= WARN__PATH_DISTANCE_TOO_SMALL_1;
    if (d__meter < 1e3)
        *warnings |= WARN__PATH_DISTANCE_TOO_SMALL_2;
    if (d__meter > 1000e3)
        *warnings |= WARN__PATH_DISTANCE_TOO_BIG_1;
    if (d__meter > 2000e3)
        *warnings |= WARN__PATH_DISTANCE_TOO_BIG_2;

    // if the path distance is less than the maximum smooth earth line of sight distance...
    if (d__meter < d_sML__meter) {
        // Compute the diffraction loss at the maximum smooth earth line of sight distance
        const double A_sML__db = d_sML__meter * M_d + A_d0__db;

        // [ERL 79-ITS 67, Eqn 3.16a], in meters instead of km and with MIN() part below
        double d_0__meter = 0.04 * f__mhz * h_e__meter[0] * h_e__meter[1];

        double d_1__meter;
        if (A_d0__db >= 0.0) {
            d_0__meter = MIN(d_0__meter, 0.5 * d_ML__meter);                // other part of [ERL 79-ITS 67, Eqn 3.16a]
            d_1__meter = d_0__meter + 0.25 * (d_ML__meter - d_0__meter);    // [ERL 79-ITS 67, Eqn 3.16d]
        } else
            d_1__meter = MAX(-A_d0__db / M_d, 0.25 * d_ML__meter);

        const double A_1__db = LineOfSightLoss(d_1__meter, h_e__meter, Z_g, delta_h__meter, M_d, A_d0__db, d_sML__meter, f__mhz);

        int flag = 0;

        double kHat_1__db_per_meter = 0;
        double kHat_2__db_per_meter = 0;

        if (d_0__meter < d_1__meter) {
            const double A_0__db = LineOfSightLoss(d_0__meter, h_e__meter, Z_g, delta_h__meter, M_d, A_d0__db, d_sML__meter, f__mhz);

            const double q = log(d_sML__meter / d_0__meter);

            // [ERL 79-ITS 67, Eqn 3.20]
            kHat_2__db_per_meter = MAX(0.0, ((d_sML__meter - d_0__meter) * (A_1__db - A_0__db) - (d_1__meter - d_0__meter) * (A_sML__db - A_0__db)) / ((d_sML__meter - d_0__meter) * log(d_1__meter / d_0__meter) - (d_1__meter - d_0__meter) * q));

            flag = A_d0__db > 0.0 || kHat_2__db_per_meter > 0.0;

            if (flag) {
                // [ERL 79-ITS 67, Eqn 3.21]
                kHat_1__db_per_meter = (A_sML__db - A_0__db - kHat_2__db_per_meter * q) / (d_sML__meter - d_0__meter);

                if (kHat_1__db_per_meter < 0.0) {
                    kHat_1__db_per_meter = 0.0;
                    kHat_2__db_per_meter = DIM(A_sML__db, A_0__db) / q;

                    if (kHat_2__db_per_meter == 0.0)
                        kHat_1__db_per_meter = M_d;
                }
            }
        }

        if (!flag) {
            kHat_1__db_per_meter = DIM(A_sML__db, A_1__db) / (d_sML__meter - d_1__meter);
            kHat_2__db_per_meter = 0.0;

            if (kHat_1__db_per_meter == 0.0)
                kHat_1__db_per_meter = M_d;
        }

        const double A_o__db = A_sML__db - kHat_1__db_per_meter * d_sML__meter - kHat_2__db_per_meter * log(d_sML__meter);

        // [ERL 79-ITS 67, Eqn 3.19]
        *A_ref__db = A_o__db + kHat_1__db_per_meter * d__meter + kHat_2__db_per_meter * log(d__meter);
        *propmode = MODE__LINE_OF_SIGHT;
    } else { // this is a trans-horizon path
        // select to points far into the troposcatter region
        const double d_5__meter = d_ML__meter + 200e3;
        const double d_6__meter = d_ML__meter + 400e3;

        // Compute the troposcatter loss at the two distances
        double h0 = -1;
        const double A_6__db = TroposcatterLoss(d_6__meter, theta_hzn, d_hzn__meter, h_e__meter, a_e__meter, N_s, f__mhz, theta_los, &h0);
        const double A_5__db = TroposcatterLoss(d_5__meter, theta_hzn, d_hzn__meter, h_e__meter, a_e__meter, N_s, f__mhz, theta_los, &h0);

        double M_s, A_s0__db, d_x__meter;

        // if we got a reasonable prediction value back...
        if (A_5__db < 1000.0) {
            // Compute the slope of the troposcatter line
            M_s = (A_6__db - A_5__db) / 200e3;

            // Find the diffraction-troposcatter transition distance
            d_x__meter = MAX(MAX(d_sML__meter, d_ML__meter + 1.088 * pow(pow(a_e__meter, 2) / f__mhz, 1.0 / 3.0) * log(f__mhz)), (A_5__db - A_d0__db - M_s * d_5__meter) / (M_d - M_s));

            // Compute the intercept of the troposcatter line
            A_s0__db = (M_d - M_s) * d_x__meter + A_d0__db;
        } else {
            // troposcatter gives no real results - so use diffraction line parameters for tropo line
            M_s = M_d;
            A_s0__db = A_d0__db;
            d_x__meter = 10e6;
        }

        // Determine if its diffraction or troposcatter and compute the loss
        if (d__meter > d_x__meter) {
            *A_ref__db = M_s * d__meter + A_s0__db;
            *propmode = MODE__TROPOSCATTER;
        } else {
            *A_ref__db = M_d * d__meter + A_d0__db;
            *propmode = MODE__DIFFRACTION;
        }
    }

    // Don't allow a negative loss
    *A_ref__db = MAX(*A_ref__db, 0.0);

    return SUCCESS;
}

/*=============================================================================
 |
 |  Description:  Extract parameters from the terrain pfl
 |
 |        Input:  pfl[2]            - Terrain data in pfl format
 |                gamma_e           - Effective earth curvature
 |                h__meter[2]       - Terminal structural heights, in meters
 |
 |      Outputs:  theta_hzn[2]      - Terminal horizon angles
 |                d_hzn__meter[2]   - Terminal horizon distances, in meters
 |                h_e__meter[2]     - Effective terminal heights, in meters
 |                delta_h__meter    - Terrain irregularity parameter
 |                d__meter          - Path distance, in meters
 |
 |      Returns:  [None]
 |
 *===========================================================================*/
static void QuickPfl(const double pfl[], double gamma_e, const double h__meter[2], double theta_hzn[2],
    double d_hzn__meter[2], double h_e__meter[2], double *delta_h__meter, double *d__meter) {
    double fit_tx, fit_rx, q;
    double d_start__meter;
    double d_end__meter;

    *d__meter = pfl[0] * pfl[1];

    const int np = (int)pfl[0];

    const double a_e__meter = 1 / gamma_e;        // effective earth radius

    FindHorizons(pfl, a_e__meter, h__meter, theta_hzn, d_hzn__meter);

    // "In our own work we have sometimes said that consideration of terrain elevations should begin at a point about 15 times the tower height"
    //      - [Hufford, 1982] Page 25
    d_start__meter = MIN(15.0 * h__meter[0], 0.1 * d_hzn__meter[0]);             // take lesser: 10% of horizon distance or 15x terminal height
    d_end__meter = *d__meter - MIN(15.0 * h__meter[1], 0.1 * d_hzn__meter[1]);   // << ditto, but measured from the far end of the link >>

    *delta_h__meter = ComputeDeltaH(pfl, d_start__meter, d_end__meter);

    if (d_hzn__meter[0] + d_hzn__meter[1] > 1.5 * *d__meter) {
        // The combined horizon distance is at least 50% larger than the total path distance
        //  -> so we are well within the line-of-sight range

        LinearLeastSquaresFit(pfl, d_start__meter, d_end__meter, &fit_tx, &fit_rx);

        h_e__meter[0] = h__meter[0] + fdim(pfl[2], fit_tx);
        h_e__meter[1] = h__meter[1] + fdim(pfl[np + 2], fit_rx);

        for (int i = 0; i < 2; i++)
            d_hzn__meter[i] = sqrt(2.0 * h_e__meter[i] * a_e__meter) * exp(-0.07 * sqrt(*delta_h__meter / MAX(h_e__meter[i], 5.0)));

        const double combined_horizons__meter = d_hzn__meter[0] + d_hzn__meter[1];
        if (combined_horizons__meter <= *d__meter) {
            q = pow(*d__meter / combined_horizons__meter, 2);

            for (int i = 0; i < 2; i++) {
                h_e__meter[i] = h_e__meter[i] * q;
                d_hzn__meter[i] = sqrt(2.0 * h_e__meter[i] * a_e__meter) * exp(-0.07 * sqrt(*delta_h__meter / MAX(h_e__meter[i], 5.0)));
            }
        }

        for (int i = 0; i < 2; i++) {
            q = sqrt(2.0 * h_e__meter[i] * a_e__meter);
            theta_hzn[i] = (0.65 * *delta_h__meter * (q / d_hzn__meter[i] - 1.0) - 2.0 * h_e__meter[i]) / q;
        }
    } else {
        double dummy = 0;

        LinearLeastSquaresFit(pfl, d_start__meter, 0.9 * d_hzn__meter[0], &fit_tx, &dummy);
        h_e__meter[0] = h__meter[0] + fdim(pfl[2], fit_tx);

        LinearLeastSquaresFit(pfl, *d__meter - 0.9 * d_hzn__meter[1], d_end__meter, &dummy, &fit_rx);
        h_e__meter[1] = h__meter[1] + fdim(pfl[np + 2], fit_rx);
    }
}

/*=============================================================================
 |
 |  Description:  Compute sigma_h
 |
 |        Input:  delta_h__meter - Terrain irregularity parameter
 |
 |      Outputs:  [None]
 |
 |      Returns:  sigma_h_meter  - sigma_h
 |
 *===========================================================================*/
static double SigmaHFunction(double delta_h__meter) {
    // "RMS deviation of terrain and terrain clutter within the limits of the first Fresnel zone in the dominant reflecting plane"
    // [ERL 79-ITS 67, Eqn 3.6a]
    return 0.78 * delta_h__meter * exp(-0.5 * pow(delta_h__meter, 0.25));
}

/*=============================================================================
 |
 |  Description:  Compute the smooth earth diffraction loss using the
 |                Vogler 3-radii method
 |
 |        Input:  d__meter          - Path distance, in meters
 |                f__mhz            - Frequency, in MHz
 |                a_e__meter        - Effective earth radius, in meters
 |                theta_los         - Angular distance of line-of-sight region
 |                d_hzn__meter[2]   - Horizon distances, in meters
 |                h_e__meter[2]     - Effective terminal heights, in meters
 |                Z_g               - Complex ground impedance
 |
 |      Outputs:  [None]
 |
 |      Returns:  A_r__db           - Smooth-earth diffraction loss, in dB
 |
 *===========================================================================*/
static double SmoothEarthDiffraction(double d__meter, double f__mhz, double a_e__meter, double theta_los,
    const double d_hzn__meter[2], const double h_e__meter[2], double complex Z_g) {
    double a__meter[3];
    double d__km[3];
    double F_x__db[2];
    double K[3];
    double B_0[3];
    double x__km[3];
    double C_0[3];

    const double theta_nlos = d__meter / a_e__meter - theta_los;                    // [Algorithm, Eqn 4.12]
    const double d_ML__meter = d_hzn__meter[0] + d_hzn__meter[1];                   // Maximum line-of-sight distance for actual path

    // compute 3 radii
    a__meter[0] = (d__meter - d_ML__meter) / (d__meter / a_e__meter - theta_los);   // which is a_e__meter when theta_los = d_ML__meter / a_e__meter
    a__meter[1] = 0.5 * pow(d_hzn__meter[0], 2) / h_e__meter[0];                    // Compute the radius of the effective earth for terminal j using[Volger 1964, Eqn 3] re - arranged
    a__meter[2] = 0.5 * pow(d_hzn__meter[1], 2) / h_e__meter[1];                    // Compute the radius of the effective earth for terminal j using[Volger 1964, Eqn 3] re - arranged

    d__km[0] = (a__meter[0] * theta_nlos) / 1000.0;                                 // angular distance of the "diffraction path"
    d__km[1] = d_hzn__meter[0] / 1000.0;
    d__km[2] = d_hzn__meter[1] / 1000.0;

    for (int i = 0; i < 3; i++) {
        // C_0 is the ratio of the 4/3 earth to effective earth (technically Vogler 1964 ratio is 4/3 to effective earth k value), all raised to the (1/3) power.
        // C_0 = (4 / 3k) ^ (1 / 3) [Vogler 1964, Eqn 2]
        C_0[i] = pow((4.0 / 3.0) * a_0__meter / a__meter[i], THIRD);

        // [Vogler 1964, Eqn 6a / 7a]
        K[i] = 0.017778 * C_0[i] * pow(f__mhz, -THIRD) / cabs(Z_g);

        // compute B_0 for each radius
        // [Vogler 1964, Fig 4]
        B_0[i] = 1.607 - K[i];
    }

    // compute x__km for each radius [Vogler 1964, Eqn 2]
    x__km[1] = B_0[1] * pow(C_0[1], 2) * pow(f__mhz, THIRD) * d__km[1];
    x__km[2] = B_0[2] * pow(C_0[2], 2) * pow(f__mhz, THIRD) * d__km[2];
    x__km[0] = B_0[0] * pow(C_0[0], 2) * pow(f__mhz, THIRD) * d__km[0] + x__km[1] + x__km[2];

    // compute height gain functions
    F_x__db[0] = HeightFunction(x__km[1], K[1]);
    F_x__db[1] = HeightFunction(x__km[2], K[2]);

    // compute distance function
    const double G_x__db = 0.05751 * x__km[0] - 10.0 * log10(x__km[0]);             // [TN101, Eqn 8.4] & [Volger 1964, Eqn 13]

    return G_x__db - F_x__db[0] - F_x__db[1] - 20;                                  // [Algorithm, Eqn 4.20] & [Volger 1964]
}

/*=============================================================================
 |
 |  Description:  Height Function, F(x, K) for smooth earth diffraction
 |
 |        Input:  x__km          - Normalized distance, in meters
 |                K              - K value
 |
 |      Outputs:  [None]
 |
 |      Returns:  F(x, K)        - in dB
 |
 *===========================================================================*/
static double HeightFunction(double x__km, double K) {
    double w;
    double result;

    if (x__km < 200.0) {
        w = -log(K);

        if (K < 1e-5 || x__km * pow(w, 3) > 5495.0) {
            result = -117.0;

            if (x__km > 1.0)
                result = 17.372 * log(x__km) + result;
        } else
            result = 2.5e-5 * pow(x__km, 2) / K - 8.686 * w - 15.0;
    } else {
        result = 0.05751 * x__km - 4.343 * log(x__km);

        if (x__km < 2000) {
            w = 0.0134 * x__km * exp(-0.005 * x__km);
            result = (1.0 - w) * result + w * (17.372 * log(x__km) - 117.0);
        }
    }

    return result;
}

/*=============================================================================
 |
 |  Description:  Compute delta_h_d
 |
 |        Input:  d__meter       - Path distance, in meters
 |                delta_h__meter - Terrain irregularity parameter
 |
 |      Outputs:  [None]
 |
 |      Returns:  delta_h_d      - Terrain irregularity of path
 |
 *===========================================================================*/
static double TerrainRoughness(double d__meter, double delta_h__meter) {
    // [ERL 79-ITS 67, Eqn 3], with distance in meters instead of kilometers
    return delta_h__meter * (1.0 - 0.8 * exp(-d__meter / 50e3));
}

/*=============================================================================
 |
 |  Description:  The attenuation function, F(th * d)
 |
 |        Input:  td             - theta * distance
 |
 |      Outputs:  [None]
 |
 |      Returns:  F()            - in dB
 |
 *===========================================================================*/
static double FFunction(double td) {
    // constants from [Algorithm, 6.9]
    const double a[] = { 133.4, 104.6, 71.8 };
    const double b[] = { 0.332e-3, 0.212e-3, 0.157e-3 };
    const double c[] = { -10, -2.5, 5 };

    int i;

    // select the set of values to use
    if (td <= 10e3)         // <= 10 km
        i = 0;
    else if (td <= 70e3)    // 10 km to 70 km
        i = 1;
    else                    // > 70 km
        i = 2;

    const double F_0 = a[i] + b[i] * td + c[i] * log10(td);  // [Algorithm, 6.9]

    return F_0;
}

/*=============================================================================
 |
 |  Description:  Troposcatter loss
 |
 |        Input:  d__meter          - Path distance, in meters
 |                theta_hzn[2]      - Terminal horizon angles
 |                d_hzn__meter[2]   - Terminal horizon distances, in meters
 |                h_e__meter[2]     - Effective terminal heights, in meters
 |                a_e__meter        - Effective earth radius, in meters
 |                N_s               - Surface refractivity, in N-Units
 |                f__mhz            - Frequency, in MHz
 |                theta_los         - Angular distance of LOS region
 |
 |      Outputs:  h0                - H_0() value
 |
 |      Returns:  F()               - in dB
 |
 *===========================================================================*/
static double TroposcatterLoss(double d__meter, const double theta_hzn[2], const double d_hzn__meter[2], const double h_e__meter[2],
    double a_e__meter, double N_s, double f__mhz, double theta_los, double *h0) {
    double H_0;

    // wavenumber, k
    const double wn = f__mhz / 47.7;

    if (*h0 > 15.0)     // short-circuit calculations if already greater than 15 dB
        H_0 = *h0;
    else {
        double ad = d_hzn__meter[0] - d_hzn__meter[1];
        double rr = h_e__meter[1] / h_e__meter[0];

        if (ad < 0.0) {      // ensure correct frame of reference
            ad = -ad;
            rr = 1.0 / rr;
        }

        const double theta = theta_hzn[0] + theta_hzn[1] + d__meter / a_e__meter;  // angular distance, in radians

        // [TN101, Eqn 9.4a]
        const double r_1 = 2.0 * wn * theta * h_e__meter[0];
        const double r_2 = 2.0 * wn * theta * h_e__meter[1];

        if (r_1 < 0.2 && r_2 < 0.2)
            return 1001;                // "If both r_1 and r_2 are less than 0.2 the function A_scat is not defined (or is infinite)" [Algorithm, page 11]

        double s = (d__meter - ad) / (d__meter + ad);       // asymmetry parameter

        // "In all of this, we truncate the values of s and q at 0.1 and 10" [Algorithm, page 16]
        const double q = MIN(MAX(0.1, rr / s), 10.0);       // TN101, Eqn 9.5
        s = MAX(0.1, s);                                    // TN101, Eqn 9.5

        double h_0__meter = (d__meter - ad) * (d__meter + ad) * theta * 0.25 / d__meter;   // height of cross-over, [Algorithm, 4.66] [TN101v1, 9.3b]

        const double Z_0__meter = 1.7556e3;             // Scale height, [Algorithm, 4.67]
        const double Z_1__meter = 8.0e3;                // [Algorithm, 4.67]
        const double eta_s = (h_0__meter / Z_0__meter) * (1.0 + (0.031 - N_s * 2.32e-3 + pow(N_s, 2) * 5.67e-6) * exp(-pow(MIN(1.7, h_0__meter / Z_1__meter), 6)));  // Scattering efficiency factor, eta_s [TN101 Eqn 9.3a]

        const double H_00 = (H0Function(r_1, eta_s) + H0Function(r_2, eta_s)) / 2;                        // First term in TN101v1, Eqn 9.5
        const double Delta_H_0 = MIN(H_00, 6.0 * (0.6 - log10(MAX(eta_s, 1.0))) * log10(s) * log10(q));

        H_0 = H_00 + Delta_H_0;                             // TN101, Eqn 9.5
        H_0 = MAX(H_0, 0.0);                                // "If Delta_H_0 would make H_0 negative, use H_0 = 0" [TN101v1, p9.4]

        if (eta_s < 1.0)    // if <=1, interpolate with the special case of eta_s = 0
            H_0 = eta_s * H_0 + (1.0 - eta_s) * 10 * log10(pow((1.0 + SQRT2 / r_1) * (1.0 + SQRT2 / r_2), 2) * (r_1 + r_2) / (r_1 + r_2 + 2 * SQRT2));

        // "If, at d_5, calculations show that H_0 will exceed 15 dB, they are replaced by the value it has at d_6" [Algorithm, page 12]
        if (H_0 > 15.0 && *h0 >= 0.0)
            H_0 = *h0;
    }

    *h0 = H_0;
    const double th = d__meter / a_e__meter - theta_los;

    const double D_0__meter = 40e3;   // [Algorithm, 6.8]
    const double H__meter = 47.7;     // [Algorithm, 4.63]
    return FFunction(th * d__meter) + 10 * log10(wn * H__meter * pow(th, 4)) - 0.1 * (N_s - 301.0) * exp(-th * d__meter / D_0__meter) + H_0;    // [Algorithm, 4.63]
}

/*=============================================================================
 |
 |  Description:  Perform input parameter validation.  This function only
 |                applies to the set of variables common to both ITM
 |                point-to-point mode and area mode.
 |
 |        Input:  h_tx__meter    - Structural height of the TX, in meters
 |				  h_rx__meter    - Structural height of the RX, in meters
 |                climate        - Radio climate enum
 |                time           - Time percentage, 0 < time < 100
 |                location       - Location percentage, 0 < location < 100
 |                situation      - Situation percentage, 0 < situation < 100
 |                N_0            - Refractivity, in N-Units
 |                f__mhz         - Frequency, in MHz
 |                pol            - Polarization
 |                epsilon        - Relative permittivity
 |                sigma          - Conductivity
 |                mdvar          - Mode of variability
 |
 |      Outputs:  warnings       - Warning messages
 |
 |      Returns:  [None]
 |
 *===========================================================================*/
static int ValidateInputs(double h_tx__meter, double h_rx__meter, int climate, double time,
    double location, double situation, double N_0, double f__mhz, int pol,
    double epsilon, double sigma, int mdvar, long *warnings) {
    if (h_tx__meter < 1.0 || h_tx__meter > 1000.0)
        *warnings |= WARN__TX_TERMINAL_HEIGHT;

    if (h_tx__meter < 0.5 || h_tx__meter > 3000.0)
        return ERROR__TX_TERMINAL_HEIGHT;

    if (h_rx__meter < 1.0 || h_rx__meter > 1000.0)
        *warnings |= WARN__RX_TERMINAL_HEIGHT;

    if (h_rx__meter < 0.5 || h_rx__meter > 3000.0)
        return ERROR__RX_TERMINAL_HEIGHT;

    if (climate != CLIMATE__EQUATORIAL &&
        climate != CLIMATE__CONTINENTAL_SUBTROPICAL &&
        climate != CLIMATE__MARITIME_SUBTROPICAL &&
        climate != CLIMATE__DESERT &&
        climate != CLIMATE__CONTINENTAL_TEMPERATE &&
        climate != CLIMATE__MARITIME_TEMPERATE_OVER_LAND &&
        climate != CLIMATE__MARITIME_TEMPERATE_OVER_SEA)
        return ERROR__INVALID_RADIO_CLIMATE;

    if (N_0 < 250 || N_0 > 400)
        return ERROR__REFRACTIVITY;

    if (f__mhz < 40.0 || f__mhz > 10000.0)
        *warnings |= WARN__FREQUENCY;

    if (f__mhz < 20 || f__mhz > 20000)
        return ERROR__FREQUENCY;

    if (pol != POLARIZATION__HORIZONTAL &&
        pol != POLARIZATION__VERTICAL)
        return ERROR__POLARIZATION;

    if (epsilon < 1)
        return ERROR__EPSILON;

    if (sigma <= 0)
        return ERROR__SIGMA;

    if ((mdvar < 0) ||
        (mdvar > 3 && mdvar < 10) ||
        (mdvar > 13 && mdvar < 20) ||
        (mdvar > 23 && mdvar < 30) ||
        (mdvar > 33))
        return ERROR__MDVAR;

    if (situation <= 0 || situation >= 100)
        return ERROR__INVALID_SITUATION;

    if (time <= 0 || time >= 100)
        return ERROR__INVALID_TIME;

    if (location <= 0 || location >= 100)
        return ERROR__INVALID_LOCATION;

    return SUCCESS;
}

/*=============================================================================
 |
 |  Description:  Curve helper function for TN101v2 Eqn III.69 & III.70
 |
 |        Input:  c1, c2, x1, x2, x3    - Curve fit parameters
 |                d_e__metre            - Effective distance, in meters
 |
 |      Outputs:  [None]
 |
 |      Returns:  Curve value           - in dB
 |
 *===========================================================================*/
static double Curve(double c1, double c2, double x1, double x2, double x3, double d_e__meter) {
    return (c1 + c2 / (1.0 + pow((d_e__meter - x2) / x3, 2))) * (pow(d_e__meter / x1, 2)) / (1.0 + (pow(d_e__meter / x1, 2)));
}

/*=============================================================================
 |
 |  Description:  Compute the variability loss
 |
 |        Input:  time           - Time percentage, 0 < time < 1
 |                location       - Location percentage, 0 < location < 1
 |                situation      - Situation percentage, 0 < situation < 1
 |                h_e__meter[2]  - Effective antenna heights, in meters
 |                delta_h__meter - Terrain irregularity parameter
 |                f__mhz         - Frequency, in MHz
 |                d__meter       - Path distance, in meters
 |                A_ref__db      - Reference attenuation, in dB
 |                climate        - Radio climate enum
 |                mdvar          - Mode of variability
 |
 |      Outputs:  warnings       - Warning flags
 |
 |      Returns:  F()            - in dB
 |
 *===========================================================================*/
static double Variability(double time, double location, double situation, const double h_e__meter[2], double delta_h__meter,
    double f__mhz, double d__meter, double A_ref__db, int climate, int mdvar, long *warnings) {
    // Asymptotic values from TN101, Fig 10.13
    // -> approximate to TN101v2 Eqn III.69 & III.70
    // -> to describe the curves for each climate
    const double all_year[5][7] = {
        {  -9.67,   -0.62,    1.26,   -9.21,   -0.62,   -0.39,      3.15 },
        {  12.7,     9.19,   15.5,     9.05,    9.19,    2.86,   857.9   },
        { 144.9e3, 228.9e3, 262.6e3,  84.1e3, 228.9e3, 141.7e3, 2222.e3  },
        { 190.3e3, 205.2e3, 185.2e3, 101.1e3, 205.2e3, 315.9e3,  164.8e3 },
        { 133.8e3, 143.6e3,  99.8e3,  98.6e3, 143.6e3, 167.4e3,  116.3e3 }
    };

    const double bsm1[] = { 2.13,      2.66,    6.11,     1.98,   2.68,    6.86,    8.51 };
    const double bsm2[] = { 159.5,     7.67,    6.65,    13.11,   7.16,   10.38,  169.8 };
    const double xsm1[] = { 762.2e3, 100.4e3, 138.2e3, 139.1e3,  93.7e3, 187.8e3, 609.8e3 };
    const double xsm2[] = { 123.6e3, 172.5e3, 242.2e3, 132.7e3, 186.8e3, 169.6e3, 119.9e3 };
    const double xsm3[] = { 94.5e3,  136.4e3, 178.6e3, 193.5e3, 133.5e3, 108.9e3, 106.6e3 };

    const double bsp1[] = { 2.11, 6.87, 10.08, 3.68, 4.75, 8.58, 8.43 };
    const double bsp2[] = { 102.3, 15.53, 9.60, 159.3, 8.12, 13.97, 8.19 };
    const double xsp1[] = { 636.9e3, 138.7e3, 165.3e3, 464.4e3, 93.2e3, 216.0e3, 136.2e3 };
    const double xsp2[] = { 134.8e3, 143.7e3, 225.7e3, 93.1e3, 135.9e3, 152.0e3, 188.5e3 };
    const double xsp3[] = { 95.6e3, 98.6e3, 129.7e3, 94.2e3, 113.4e3, 122.7e3, 122.9e3 };

    const double C_D[] = { 1.224, 0.801, 1.380, 1.000, 1.224, 1.518, 1.518 };	// [Algorithm, Table 5.1], C_d
    const double z_D[] = { 1.282, 2.161, 1.282, 20.0, 1.282, 1.282, 1.282 };	// [Algorithm, Table 5.1], z_d

    const double bfm1[] = { 1.0, 1.0, 1.0, 1.0, 0.92, 1.0, 1.0 };
    const double bfm2[] = { 0.0, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0 };
    const double bfm3[] = { 0.0, 0.0, 0.0, 0.0, 1.77, 0.0, 0.0 };

    const double bfp1[] = { 1.0, 0.93, 1.0, 0.93, 0.93, 1.0, 1.0 };
    const double bfp2[] = { 0.0, 0.31, 0.0, 0.19, 0.31, 0.0, 0.0 };
    const double bfp3[] = { 0.0, 2.00, 0.0, 1.79, 2.00, 0.0, 0.0 };

    double z_T = InverseComplementaryCumulativeDistributionFunction(time / 100);
    double z_L = InverseComplementaryCumulativeDistributionFunction(location / 100);
    const double z_S = InverseComplementaryCumulativeDistributionFunction(situation / 100);

    int climate_idx = climate; // Create an internal copy for modification
    climate_idx--; // 0-based indexes

    const double wn = f__mhz / 47.7;

    // compute the effective distance
    const double d_ex__meter = sqrt(2 * a_9000__meter * h_e__meter[0]) + sqrt(2 * a_9000__meter * h_e__meter[1]) + pow((575.7e12 / wn), THIRD);  // [Algorithm, Eqn 5.3]

    double d_e__meter;
    if (d__meter < d_ex__meter)
        d_e__meter = 130e3 * d__meter / d_ex__meter;
    else
        d_e__meter = 130e3 + d__meter - d_ex__meter;

    //////////////////////////////////
    // situation variability calcs

    // if mdvar >= 20, then "Direct situation variability is to be eliminated as it should when
    //                       considering interference problems.  Note that there may still be a
    //                       small residual situation variability" [Hufford, 1982]
    int mdvar_internal = mdvar;  // Create an internal copy to modify
    const int plus20 = mdvar_internal >= 20;
    if (plus20)
        mdvar_internal -= 20;

    double sigma_S;
    if (plus20)
        sigma_S = 0.0;
    else {
        double D__meter = 100e3;                                // Scale distance, D = 100 km
        sigma_S = 5.0 + 3.0 * exp(-d_e__meter / D__meter);      // [Algorithm, Eqn 5.10]
    }

    //
    //////////////////////////////////


    const int plus10 = mdvar_internal >= 10;
    if (plus10)
        mdvar_internal -= 10;

    const double V_med__db = Curve(all_year[0][climate_idx], all_year[1][climate_idx], all_year[2][climate_idx], all_year[3][climate_idx], all_year[4][climate_idx], d_e__meter);

    if (mdvar_internal == SINGLE_MESSAGE_MODE) {
        z_T = z_S;
        z_L = z_S;
    } else if (mdvar_internal == ACCIDENTAL_MODE)
        z_L = z_S;
    else if (mdvar_internal == MOBILE_MODE)
        z_L = z_T;
    // else using Broadcast Mode (no additional operations)

    if (fabs(z_T) > 3.10 || fabs(z_L) > 3.10 || fabs(z_S) > 3.10)
        *warnings |= WARN__EXTREME_VARIABILITIES;

    //////////////////////////////////
    // location variability calcs

    double sigma_L;
    if (plus10)
        sigma_L = 0.0;
    else {
        const double delta_h_d__meter = TerrainRoughness(d__meter, delta_h__meter);

        sigma_L = 10.0 * wn * delta_h_d__meter / (wn * delta_h_d__meter + 13.0);    // Context of [Algorithm, Eqn 5.9]
    }
    const double Y_L = sigma_L * z_L;

    //
    //////////////////////////////////

    //////////////////////////////////
    // time variability calcs

    const double q = log(0.133 * wn);
    const double g_minus = bfm1[climate_idx] + bfm2[climate_idx] / (pow(bfm3[climate_idx] * q, 2) + 1.0);
    const double g_plus = bfp1[climate_idx] + bfp2[climate_idx] / (pow(bfp3[climate_idx] * q, 2) + 1.0);

    const double sigma_T_minus = Curve(bsm1[climate_idx], bsm2[climate_idx], xsm1[climate_idx], xsm2[climate_idx], xsm3[climate_idx], d_e__meter) * g_minus;
    const double sigma_T_plus = Curve(bsp1[climate_idx], bsp2[climate_idx], xsp1[climate_idx], xsp2[climate_idx], xsp3[climate_idx], d_e__meter) * g_plus;

    const double sigma_TD = C_D[climate_idx] * sigma_T_plus;
    const double tgtd = (sigma_T_plus - sigma_TD) * z_D[climate_idx];

    double sigma_T;
    if (z_T < 0.0)
        sigma_T = sigma_T_minus;
    else if (z_T <= z_D[climate_idx])
        sigma_T = sigma_T_plus;
    else
        sigma_T = sigma_TD + tgtd / z_T;
    const double Y_T = sigma_T * z_T;

    //
    /////////////////////////////////

    const double Y_S_temp = pow(sigma_S, 2) + pow(Y_T, 2) / (7.8 + pow(z_S, 2)) + pow(Y_L, 2) / (24.0 + pow(z_S, 2));  // Part of [Algorithm, Eqn 5.11]

    double Y_R, Y_S;
    if (mdvar_internal == SINGLE_MESSAGE_MODE) {
        Y_R = 0.0;
        Y_S = sqrt(pow(sigma_T, 2) + pow(sigma_L, 2) + Y_S_temp) * z_S;
    } else if (mdvar_internal == ACCIDENTAL_MODE) {
        Y_R = Y_T;
        Y_S = sqrt(pow(sigma_L, 2) + Y_S_temp) * z_S;
    } else if (mdvar_internal == MOBILE_MODE) {
        Y_R = sqrt(pow(sigma_T, 2) + pow(sigma_L, 2)) * z_T;
        Y_S = sqrt(Y_S_temp) * z_S;
    } else { // BROADCAST_MODE
        Y_R = Y_T + Y_L;
        Y_S = sqrt(Y_S_temp) * z_S;
    }

    double result = A_ref__db - V_med__db - Y_R - Y_S;

    // [Algorithm, Eqn 52]
    if (result < 0.0)
        result = result * (29.0 - result) / (29.0 - 10.0 * result);

    return result;
}

/*=============================================================================
 |
 |  Description:  The ITS Irregular Terrain Model (ITM).  This function
 |               exposes point-to-point mode functionality, with variability
 |               specified with time/location/situation (TLS).
 |
 |        Input:  h_tx__meter       - Structural height of the TX, in meters
 |                h_rx__meter       - Structural height of the RX, in meters
 |                pfl[2]            - Terrain data, in PFL format
 |                climate           - Radio climate
 |                N_0               - Refractivity, in N-Units
 |                f__mhz            - Frequency, in MHz
 |                pol               - Polarization
 |                epsilon           - Relative permittivity
 |                sigma             - Conductivity
 |                mdvar             - Mode of variability
 |                time              - Time percentage, 0 < time < 100
 |                location          - Location percentage, 0 < location < 100
 |                situation         - Situation percentage, 0 < situation < 100
 |
 |      Outputs:  A__db             - Basic transmission loss, in dB
 |                warnings          - Warning flags
 |                interValues       - Struct of intermediate values
 |
 |      Returns:  error             - Error code
 |
 *===========================================================================*/
static int ITM_P2P_TLS_Ex(const double h_tx__meter, const double h_rx__meter, const double pfl[], const int climate, const double N_0,
    const double f__mhz, const int pol, const double epsilon, const double sigma, const int mdvar, const double time, const double location,
    const double situation, double *A__db, long *warnings, IntermediateValues *interValues) {
    double N_s;                 // Surface refractivity, in N-Units
    double gamma_e;             // Curvature of the effective earth
    double delta_h__meter;      // Terrain irregularity parameter
    double d__meter;            // Path distance, in meters
    double complex Z_g;         // Ground impedance
    double theta_hzn[2];        // Terminal horizon angles
    double d_hzn__meter[2];     // Terminal horizon distances
    double h_e__meter[2];       // Terminal effective heights

    *warnings = NO_WARNINGS;    // Initialize to no warnings

    // initial input validation check - some validation occurs later in calculations
    int rtn = ValidateInputs(h_tx__meter, h_rx__meter, climate, time, location, situation, N_0, f__mhz, pol, epsilon, sigma, mdvar, warnings);
    if (rtn != SUCCESS)
        return rtn;

    interValues->d__km = (pfl[0] * pfl[1]) / 1000;

    const int np = (int)pfl[0];     // number of points in the pfl

    // compute the average path height, ignoring first and last 10%
    const int p10 = (int)(0.1 * np);  // 10% of np
    double h_sys__meter = 0;        // Height of the system above mean sea level

    for (int i = p10; i <= np - p10; i++)
        h_sys__meter += pfl[i + 2];

    h_sys__meter = h_sys__meter / (np - 2 * p10 + 1);

    InitializePointToPoint(f__mhz, h_sys__meter, N_0, pol, epsilon, sigma, &Z_g, &gamma_e, &N_s);

    const double h__meter[2] = { h_tx__meter, h_rx__meter };
    QuickPfl(pfl, gamma_e, h__meter, theta_hzn, d_hzn__meter, h_e__meter, &delta_h__meter, &d__meter);

    // Reference attenuation, in dB
    double A_ref__db = 0;
    int propmode = MODE__NOT_SET;
    rtn = LongleyRice(theta_hzn, f__mhz, Z_g, d_hzn__meter, h_e__meter, gamma_e, N_s, delta_h__meter, h__meter, d__meter, MODE__P2P,
        &A_ref__db, warnings, &propmode);
    if (rtn != SUCCESS)
        return rtn;

    double A_fs__db = FreeSpaceLoss(d__meter, f__mhz);

    *A__db = Variability(time, location, situation, h_e__meter, delta_h__meter, f__mhz, d__meter, A_ref__db, climate, mdvar, warnings) + A_fs__db;

    // Save off intermediate values
    interValues->A_ref__db = A_ref__db;
    interValues->A_fs__db = A_fs__db;
    interValues->delta_h__meter = delta_h__meter;
    interValues->d_hzn__meter[0] = d_hzn__meter[0];
    interValues->d_hzn__meter[1] = d_hzn__meter[1];
    interValues->h_e__meter[0] = h_e__meter[0];
    interValues->h_e__meter[1] = h_e__meter[1];
    interValues->N_s = N_s;
    interValues->theta_hzn[0] = theta_hzn[0];
    interValues->theta_hzn[1] = theta_hzn[1];
    interValues->mode = propmode;

    if (*warnings != NO_WARNINGS)
        return SUCCESS_WITH_WARNINGS;

    return SUCCESS;
}

/*=============================================================================
 |
 |  Description:  The ITS Irregular Terrain Model (ITM).  This function
 |               exposes point-to-point mode functionality, with variability
 |               specified with time/location/situation (TLS). Public entry
 |               point -- see include/longley_rice.h for full parameter docs.
 |
 *===========================================================================*/
int ITM_P2P_TLS(const double h_tx__meter, const double h_rx__meter, const double pfl[], const int climate, const double N_0, const double f__mhz,
    const int pol, const double epsilon, const double sigma, const int mdvar, const double time, const double location, const double situation,
    double *A__db, long *warnings) {
    IntermediateValues interValues;

    return ITM_P2P_TLS_Ex(h_tx__meter, h_rx__meter, pfl, climate, N_0, f__mhz, pol, epsilon, sigma, mdvar,
        time, location, situation, A__db, warnings, &interValues);
}
