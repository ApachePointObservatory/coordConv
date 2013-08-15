#pragma once
#include <cmath>
/*
Physical constants

References:
  U.S. Naval Observatory, the APPLE subroutine library for J2000.0 data
  U.S. Naval Observatory, "The Astronomical Almanac" for 1987
  S. Selby, "Standard Mathematical Tables", 15th ed, 1967, Chemical Rubber Co.
  C. Allen, "Astrophysical Quantities", 1973, Athlone Press (U. of London)
  L. Taff, "Computational Spherical Astronomy", 1981, Wiley
  R. Green, "Spherical Astronomy", 1985, Cambridge
  S. Aoki, et. al. (1983) Astron. Astrophys. 128, 263-267
*/
namespace coordConv {

    const double Pi = std::atan(1.0)*4;
    const double HoursPerDeg = 24.0 / 360.0;    // hours per degree
    const double RadPerDeg = Pi / 180.0;        // radians per degree
    const double KmPerAU = 149597871.0;         // kilometers per astr. unit (Astr. Almanac)
    const double AUPerParsec = 206264.8062470964;   // astronomical units per parsec (the APPLE subr. library)
    const double ArcsecPerDeg = 3600.0;         // arcseconds per degree
    const double SecPerDay = 24.0 * 3600.0;     // seconds per day
    const double DaysPerYear = 365.25;          // days per TT (or TAI...)  year
    const double VLight = 299792.458 / KmPerAU; // (au/sec), based on c = 299792458 m/sec (exact)
    const double DegK_DegC = 273.15;            // deg. Kelvin - deg. C (Allen, exact)
    const double MJD_UnixTime = 40587 * SecPerDay;  // MJD (seconds) - unix time
        // unix time is 0 at 1970-01-01 00:00:00; unix time = JD − 2,440,587.5 (days)
        // MJD is 0 at 1858-11-17 00:00:00;             MJD = JD − 2,400,000.5 (days)
    const double MJDJ2000 = 51544.5;            // Modified Julian Date at epoch J2000.0 noon (days)
    const double AngstromsPerMicron = 1.0e4;
    const double PascalsPerMillibar = 100.0;


    const double TT_TAI = 32.184; // TT - UTC (seconds) (Astr. Almanac)
//     const double BE_J2000 = 2000.001278; // Bessilian equiv. of J2000.0 (Taff)
//     const double JE_B1950 = 1949.999790; // Julian equiv. of B1950 (Taff)
//     const double JY_per_BY = 36524.2198781 / 36525.0;    // Julian years per Besselian (troPical) Years (Taff)
    const double SiderealPerSolar = 1.00273790934;  // Mean sidereal time units per mean solar time unit; in 1987
        // (Astr. Almanac) (defined as sid. days/UT days; equals sid. days/TAI days on ave.)

}
