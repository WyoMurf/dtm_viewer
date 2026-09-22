#ifndef LONGLEY_RICE_H
#define LONGLEY_RICE_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * A faithful C port of the NTIA/ITS Irregular Terrain Model (Longley-Rice),
 * point-to-point mode only, "TLS" (time/location/situation) variability
 * variant -- i.e. the equivalent of the reference implementation's
 * ITM_P2P_TLS() entry point. Ported mechanically from the official C++
 * reference at https://github.com/NTIA/itm (branch master), which is a US
 * federal government work (Title 15 USC 105) released to the public domain
 * with explicit permission to modify and redistribute -- see that repo's
 * LICENSE.md. Every formula, constant, and branch is kept exactly as the
 * original; this is a translation, not a reimplementation.
 *
 * This is the real multi-effect ITM algorithm (line-of-sight two-ray
 * interference, single/double knife-edge and smooth-earth diffraction,
 * troposcatter, and empirical time/location/situation variability), unlike
 * terrain_viewer.c's `--viewshed --erp` model, which is a much simpler
 * single-knife-edge-diffraction approximation. Area-prediction mode (no
 * explicit terrain profile, just a distance and a delta-h roughness
 * statistic) and the confidence/reliability (CR) variability variant were
 * both left out of this port -- only point-to-point TLS is implemented,
 * since that's the mode useful for surveying with a real digitized terrain
 * profile.
 */

/* ---- Radio climate (`climate` parameter) ---- */
#define CLIMATE__EQUATORIAL                     1
#define CLIMATE__CONTINENTAL_SUBTROPICAL        2
#define CLIMATE__MARITIME_SUBTROPICAL           3
#define CLIMATE__DESERT                         4
#define CLIMATE__CONTINENTAL_TEMPERATE          5
#define CLIMATE__MARITIME_TEMPERATE_OVER_LAND   6
#define CLIMATE__MARITIME_TEMPERATE_OVER_SEA    7

/* ---- Polarization (`pol` parameter) ---- */
#define POLARIZATION__HORIZONTAL                0
#define POLARIZATION__VERTICAL                  1

/* ---- Mode of variability (`mdvar` parameter) ----
 * Base mode selects how time/location/situation variability combine:
 *   MDVAR__SINGLE_MESSAGE_MODE - all three vary together (one prediction)
 *   MDVAR__ACCIDENTAL_MODE     - location and situation vary together
 *   MDVAR__MOBILE_MODE         - time and location vary together
 *   MDVAR__BROADCAST_MODE      - all three vary independently (typical use)
 * Add 10 to eliminate location variability, add 20 to also eliminate
 * situation variability (e.g. mdvar=23 = broadcast mode w/ situation
 * variability suppressed), per the original ITM documentation. */
#define MDVAR__SINGLE_MESSAGE_MODE              0
#define MDVAR__ACCIDENTAL_MODE                  1
#define MDVAR__MOBILE_MODE                      2
#define MDVAR__BROADCAST_MODE                   3

/* ---- Return / error codes ---- */
#define SUCCESS                                 0
#define NO_WARNINGS                             0
#define SUCCESS_WITH_WARNINGS                   1

#define ERROR__TX_TERMINAL_HEIGHT               1000  /* TX terminal height is out of range */
#define ERROR__RX_TERMINAL_HEIGHT               1001  /* RX terminal height is out of range */
#define ERROR__INVALID_RADIO_CLIMATE            1002  /* Invalid value for radio climate */
#define ERROR__INVALID_TIME                     1003  /* Time percentage is out of range */
#define ERROR__INVALID_LOCATION                 1004  /* Location percentage is out of range */
#define ERROR__INVALID_SITUATION                1005  /* Situation percentage is out of range */
#define ERROR__REFRACTIVITY                     1008  /* Refractivity is out of range */
#define ERROR__FREQUENCY                        1009  /* Frequency is out of range */
#define ERROR__POLARIZATION                     1010  /* Invalid value for polarization */
#define ERROR__EPSILON                          1011  /* Epsilon is out of range */
#define ERROR__SIGMA                            1012  /* Sigma is out of range */
#define ERROR__GROUND_IMPEDANCE                 1013  /* Imaginary part of the complex ground impedance is larger than the real part */
#define ERROR__MDVAR                            1014  /* Invalid value for mode of variability */
#define ERROR__EFFECTIVE_EARTH                  1016  /* Internally computed effective earth radius is invalid */
#define ERROR__SURFACE_REFRACTIVITY_SMALL       1021  /* Internally computed surface refractivity value is too small */
#define ERROR__SURFACE_REFRACTIVITY_LARGE       1022  /* Internally computed surface refractivity value is too large */

/* ---- Warning bit flags (OR'd into *warnings on SUCCESS_WITH_WARNINGS) ---- */
#define WARN__TX_TERMINAL_HEIGHT                0x0001  /* TX terminal height is near its limits */
#define WARN__RX_TERMINAL_HEIGHT                0x0002  /* RX terminal height is near its limits */
#define WARN__FREQUENCY                         0x0004  /* Frequency is near its limits */
#define WARN__PATH_DISTANCE_TOO_BIG_1           0x0008  /* Path distance is near its upper limit */
#define WARN__PATH_DISTANCE_TOO_BIG_2           0x0010  /* Path distance is large -- care must be taken with result */
#define WARN__PATH_DISTANCE_TOO_SMALL_1         0x0020  /* Path distance is near its lower limit */
#define WARN__PATH_DISTANCE_TOO_SMALL_2         0x0040  /* Path distance is small -- care must be taken with result */
#define WARN__TX_HORIZON_ANGLE                  0x0080  /* TX horizon angle is large -- small angle approximations could break down */
#define WARN__RX_HORIZON_ANGLE                  0x0100  /* RX horizon angle is large -- small angle approximations could break down */
#define WARN__TX_HORIZON_DISTANCE_1             0x0200  /* TX horizon distance is less than 1/10 of the smooth earth horizon distance */
#define WARN__RX_HORIZON_DISTANCE_1             0x0400  /* RX horizon distance is less than 1/10 of the smooth earth horizon distance */
#define WARN__TX_HORIZON_DISTANCE_2             0x0800  /* TX horizon distance is greater than 3 times the smooth earth horizon distance */
#define WARN__RX_HORIZON_DISTANCE_2             0x1000  /* RX horizon distance is greater than 3 times the smooth earth horizon distance */
#define WARN__EXTREME_VARIABILITIES             0x2000  /* One of the provided variabilities is located far in the tail of its distribution */
#define WARN__SURFACE_REFRACTIVITY              0x4000  /* Internally computed surface refractivity value is small -- care must be taken with result */

/*
 * Point-to-point Longley-Rice prediction, TLS variability variant.
 *
 *   h_tx__meter  - TX structural height above ground, in meters
 *   h_rx__meter  - RX structural height above ground, in meters
 *   pfl[]        - Terrain profile, in the original ITM "PFL" format:
 *                    pfl[0]     = np, the number of intervals (so np+1
 *                                 elevation samples follow)
 *                    pfl[1]     = spacing between samples, in meters
 *                    pfl[2..]   = np+1 elevations, in meters, evenly spaced
 *                                 along the great-circle path from TX to RX
 *   climate      - Radio climate, one of the CLIMATE__* values above
 *   N_0          - Surface refractivity at sea level, in N-units (typically
 *                  301; valid range 250-400)
 *   f__mhz       - Frequency, in MHz (valid range 20-20000, best accuracy
 *                  40-10000)
 *   pol          - Polarization, one of the POLARIZATION__* values above
 *   epsilon      - Relative permittivity of the ground (typically ~15)
 *   sigma        - Conductivity of the ground, in S/m (typically ~0.005)
 *   mdvar        - Mode of variability, one of the MDVAR__* values above
 *                  (optionally +10/+20, see above)
 *   time         - Time percentage, 0 < time < 100
 *   location     - Location percentage, 0 < location < 100
 *   situation    - Situation percentage, 0 < situation < 100
 *
 * Outputs:
 *   *A__db       - Basic transmission loss, in dB (already includes free
 *                  space loss -- this is the total predicted path loss)
 *   *warnings    - Bitmask of WARN__* flags describing any inputs that are
 *                  near or outside the model's comfortable operating range
 *                  (set even on SUCCESS_WITH_WARNINGS; caller decides
 *                  whether to trust the result)
 *
 * Returns SUCCESS (0), SUCCESS_WITH_WARNINGS (1, result is usable but see
 * *warnings), or one of the ERROR__* codes above if an input was invalid
 * enough that no result was computed (*A__db is left untouched in that
 * case).
 */
int ITM_P2P_TLS(const double h_tx__meter, const double h_rx__meter, const double pfl[], const int climate,
    const double N_0, const double f__mhz, const int pol, const double epsilon, const double sigma,
    const int mdvar, const double time, const double location, const double situation,
    double *A__db, long *warnings);

#ifdef __cplusplus
}
#endif

#endif /* LONGLEY_RICE_H */
