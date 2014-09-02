#!/usr/bin/env python
from __future__ import absolute_import, division

import unittest
import math
import coordConv

ArcSecPerRev = 3600.0 * 360.0

RadPerArcSec = coordConv.RadPerDeg / 3600.0

def julianEpochFromTAI(tai):
    ttDays = (tai + coordConv.TT_TAI) / float(3600 * 24)
    return 2000.0 + ((ttDays - 51544.5) / 365.25)

def besselianEpochFromTAI(tai):
    ttDays = (tai + coordConv.TT_TAI) / float(3600 * 24)
    return 1900.0 + ( ttDays - 15019.81352 ) / 365.242198781

def eqeqx (tdb):
    """The equation of the equinoxes (IAU 1994).
    
    Inputs:
    - tdb (MJD): TDB (loosely ET) as a Modified Julian Date
    
    Returns eqeqx, where:
    Greenwich apparent sidereal time = Greenwich mean sidereal time + eqeqx

    Based on Pat Wallace's EQEQX, which in turn is based on:
    IAU Resolution C7, Recommendation 3 (1994)
    Capitaine, N. & Gontier, A.-M., Astron. Astrophys., 275, 645-650 (1993)
    
    History:
    Patrick Wallace   Starlink   21 November 1994
    2002-07-11 ROwen  Converted to Python
    """
    # TDB - J2000 in centuries
    t = (tdb - 51544.5) / 36525.0

    # Longitude of the mean ascending node of the lunar orbit on the
    # ecliptic, measured from the mean equinox of date
    om = RadPerArcSec * (450160.280 + ((-5.0 * ArcSecPerRev) - 482890.539 + (7.455 + 0.008 * t) * t) * t)

    # Nutation
    dpsi, deps, eps0 = nutc(tdb)

    # Equation of the equinoxes
    return dpsi * math.cos(eps0) + RadPerArcSec * (0.00264 * math.sin(om) + 0.000063 * math.sin(om * 2))

def nutc (tdb):
    """
    Computes nutation and obliquity using the IAU 1980 theory.
    
    Warning: this is older than the current slalib function and it gives somewhat different results.
    
    Inputs:
     -  tdb: TDB date (mjd)
    
    Returns a tuple of three elements:
     -  nutation in longitude (rad)
     -  nutation in obliquity (rad)
     -  mean obliquity (rad)
     
    Based on Pat Wallace's NUTC.

    References:
    - Final report of the IAU Working Group on Nutation,
      chairman P.K.Seidelmann, 1980.
    - Kaplan,G.H., 1981, USNO circular no. 163, pA3 - 6.
    
    History:
    P.T.Wallace Starlink    September 1987
    2002-07-11 ROwen  Converted to Python
    """

    #  Units of 0.0001 arcsec to radians
    u2r = RadPerArcSec / 1.0e4


    # TDB - J2000 in centuries
    t = (tdb - coordConv.MJDJ2000) / (coordConv.DaysPerYear * 100.0)

    #
    #  fundamental arguments in the fk5 reference system
    #

    #  Mean longitude of the moon minus mean longitude of the moon's perigee
    el = RadPerArcSec * (485866.733 + ((1325 * ArcSecPerRev) + 715922.633  \
        + (31.310 + (0.064 * t)) * t) * t)

    #  Mean longitude of the sun minus mean longitude of the sun's perigee
    elp = RadPerArcSec * (1287099.804 + ((99 * ArcSecPerRev) + 1292581.224  \
        + (-0.577 - 0.012 * t) * t) * t)

    #  Mean longitude of the moon minus mean longitude of the moon's node
    f = RadPerArcSec * (335778.877 + (1342 * ArcSecPerRev + 295263.137  \
        + (-13.257 + 0.011 * t) * t) * t)

    #  Mean elongation of the moon from the sun
    d = RadPerArcSec * (1072261.307 + (1236 * ArcSecPerRev + 1105601.328  \
        + (-6.891 + 0.019 * t) * t) * t)

    #  Longitude of the mean ascending node of the lunar orbit on the
    #   ecliptic, measured from the mean equinox of date
    om = RadPerArcSec * (450160.280 + (-5 * ArcSecPerRev - 482890.539  \
        + (7.455 + 0.008 * t) * t) * t)

    #  Multiples of arguments
    el2 = el + el
    el3 = el2 + el
    elp2 = elp + elp
    f2 = f + f
    f4 = f2 + f2
    d2 = d + d
    d4 = d2 + d2
    om2 = om + om


    #
    #  series for the nutation
    #
    dp = 0
    de = 0

    #  106
    dp = dp + math.sin(elp + d)
    #  105
    dp = dp -  math.sin(f2 + d4 + om2)
    #  104
    dp = dp + math.sin(el2 + d2)
    #  103
    dp = dp -  math.sin(el - f2 + d2)
    #  102
    dp = dp -  math.sin(el + elp - d2 + om)
    #  101
    dp = dp -  math.sin(-elp + f2 + om)
    #  100
    dp = dp -  math.sin(el - f2 - d2)
    #  99
    dp = dp -  math.sin(elp + d2)
    #  98
    dp = dp -  math.sin(f2 - d + om2)
    #  97
    dp = dp -  math.sin(-f2 + om)
    #  96
    dp = dp + math.sin(-el - elp + d2 + om)
    #  95
    dp = dp + math.sin(elp + f2 + om)
    #  94
    dp = dp -  math.sin(el + f2 - d2)
    #  93
    dp = dp + math.sin(el3 + f2 - d2 + om2)
    #  92
    dp = dp + math.sin(f4 - d2 + om2)
    #  91
    dp = dp -  math.sin(el + d2 + om)
    #  90
    dp = dp -  math.sin(el2 + f2 + d2 + om2)
    #  89
    a = el2 + f2 - d2 + om
    dp = dp + math.sin(a)
    de = de -  math.cos(a)
    #  88
    dp = dp + math.sin(el - elp - d2)
    #  87
    dp = dp + math.sin(-el + f4 + om2)
    #  86
    a =  - el2 + f2 + d4 + om2
    dp = dp -  math.sin(a)
    de = de + math.cos(a)
    #  85
    a = el + f2 + d2 + om
    dp = dp -  math.sin(a)
    de = de + math.cos(a)
    #  84
    a = el + elp + f2 - d2 + om2
    dp = dp + math.sin(a)
    de = de -  math.cos(a)
    #  83
    dp = dp -  math.sin(el2 - d4)
    #  82
    a =  - el + f2 + d4 + om2
    dp = dp - 2 * math.sin(a)
    de = de + math.cos(a)
    #  81
    a =  - el2 + f2 + d2 + om2
    dp = dp + math.sin(a)
    de = de -  math.cos(a)
    #  80
    dp = dp -  math.sin(el - d4)
    #  79
    a =  - el + om2
    dp = dp + math.sin(a)
    de = de -  math.cos(a)
    #  78
    a = f2 + d + om2
    dp = dp + 2 * math.sin(a)
    de = de -  math.cos(a)
    #  77
    dp = dp + 2 * math.sin(el3)
    #  76
    a = el + om2
    dp = dp - 2 * math.sin(a)
    de = de + math.cos(a)
    #  75
    a = el2 + om
    dp = dp + 2 * math.sin(a)
    de = de -  math.cos(a)
    #  74
    a =  - el + f2 - d2 + om
    dp = dp - 2 * math.sin(a)
    de = de + math.cos(a)
    #  73
    a = el + elp + f2 + om2
    dp = dp + 2 * math.sin(a)
    de = de -  math.cos(a)
    #  72
    a =  - elp + f2 + d2 + om2
    dp = dp - 3 * math.sin(a)
    de = de + math.cos(a)
    #  71
    a = el3 + f2 + om2
    dp = dp - 3 * math.sin(a)
    de = de + math.cos(a)
    #  70
    a =  - el2 + om
    dp = dp - 2 * math.sin(a)
    de = de + math.cos(a)
    #  69
    a =  - el - elp + f2 + d2 + om2
    dp = dp - 3 * math.sin(a)
    de = de + math.cos(a)
    #  68
    a = el - elp + f2 + om2
    dp = dp - 3 * math.sin(a)
    de = de + math.cos(a)
    #  67
    dp = dp + 3 * math.sin(el + f2)
    #  66
    dp = dp - 3 * math.sin(el + elp)
    #  65
    dp = dp - 4 * math.sin(d)
    #  64
    dp = dp + 4 * math.sin(el - f2)
    #  63
    dp = dp - 4 * math.sin(elp - d2)
    #  62
    a = el2 + f2 + om
    dp = dp - 5 * math.sin(a)
    de = de + 3 * math.cos(a)
    #  61
    dp = dp + 5 * math.sin(el - elp)
    #  60
    a =  - d2 + om
    dp = dp - 5 * math.sin(a)
    de = de + 3 * math.cos(a)
    #  59
    a = el + f2 - d2 + om
    dp = dp + 6 * math.sin(a)
    de = de - 3 * math.cos(a)
    #  58
    a = f2 + d2 + om
    dp = dp - 7 * math.sin(a)
    de = de + 3 * math.cos(a)
    #  57
    a = d2 + om
    dp = dp - 6 * math.sin(a)
    de = de + 3 * math.cos(a)
    #  56
    a = el2 + f2 - d2 + om2
    dp = dp + 6 * math.sin(a)
    de = de - 3 * math.cos(a)
    #  55
    dp = dp + 6 * math.sin(el + d2)
    #  54
    a = el + f2 + d2 + om2
    dp = dp - 8 * math.sin(a)
    de = de + 3 * math.cos(a)
    #  53
    a =  - elp + f2 + om2
    dp = dp - 7 * math.sin(a)
    de = de + 3 * math.cos(a)
    #  52
    a = elp + f2 + om2
    dp = dp + 7 * math.sin(a)
    de = de - 3 * math.cos(a)
    #  51
    dp = dp - 7 * math.sin(el + elp - d2)
    #  50
    a =  - el + f2 + d2 + om
    dp = dp - 10 * math.sin(a)
    de = de + 5 * math.cos(a)
    #  49
    a = el - d2 + om
    dp = dp - 13 * math.sin(a)
    de = de + 7 * math.cos(a)
    #  48
    a =  - el + d2 + om
    dp = dp + 16 * math.sin(a)
    de = de - 8 * math.cos(a)
    #  47
    a =  - el + f2 + om
    dp = dp + 21 * math.sin(a)
    de = de - 10 * math.cos(a)
    #  46
    dp = dp + 26 * math.sin(f2)
    de = de -  math.cos(f2)
    #  45
    a = el2 + f2 + om2
    dp = dp - 31 * math.sin(a)
    de = de + 13 * math.cos(a)
    #  44
    a = el + f2 - d2 + om2
    dp = dp + 29 * math.sin(a)
    de = de - 12 * math.cos(a)
    #  43
    dp = dp + 29 * math.sin(el2)
    de = de -  math.cos(el2)
    #  42
    a = f2 + d2 + om2
    dp = dp - 38 * math.sin(a)
    de = de + 16 * math.cos(a)
    #  41
    a = el + f2 + om
    dp = dp - 51 * math.sin(a)
    de = de + 27 * math.cos(a)
    #  40
    a =  - el + f2 + d2 + om2
    dp = dp - 59 * math.sin(a)
    de = de + 26 * math.cos(a)
    #  39
    a =  - el + om
    dp = dp + (-58 - 0.1 * t) * math.sin(a)
    de = de + 32 * math.cos(a)
    #  38
    a = el + om
    dp = dp + (63 + 0.1 * t) * math.sin(a)
    de = de - 33 * math.cos(a)
    #  37
    dp = dp + 63 * math.sin(d2)
    de = de - 2 * math.cos(d2)
    #  36
    a =  - el + f2 + om2
    dp = dp + 123 * math.sin(a)
    de = de - 53 * math.cos(a)
    #  35
    a = el - d2
    dp = dp - 158 * math.sin(a)
    de = de -  math.cos(a)
    #  34
    a = el + f2 + om2
    dp = dp - 301 * math.sin(a)
    de = de + (129 - 0.1 * t) * math.cos(a)
    #  33
    a = f2 + om
    dp = dp + (-386 - 0.4 * t) * math.sin(a)
    de = de + 200 * math.cos(a)
    #  32
    dp = dp + (712 + 0.1 * t) * math.sin(el)
    de = de - 7 * math.cos(el)
    #  31
    a = f2 + om2
    dp = dp + (-2274 - 0.2 * t) * math.sin(a)
    de = de + (977 - 0.5 * t) * math.cos(a)
    #  30
    dp = dp -  math.sin(elp + f2 - d2)
    #  29
    dp = dp + math.sin(-el + d + om)
    #  28
    dp = dp + math.sin(elp + om2)
    #  27
    dp = dp -  math.sin(elp - f2 + d2)
    #  26
    dp = dp + math.sin(-f2 + d2 + om)
    #  25
    dp = dp + math.sin(el2 + elp - d2)
    #  24
    dp = dp - 4 * math.sin(el - d)
    #  23
    a = elp + f2 - d2 + om
    dp = dp + 4 * math.sin(a)
    de = de - 2 * math.cos(a)
    #  22
    a = el2 - d2 + om
    dp = dp + 4 * math.sin(a)
    de = de - 2 * math.cos(a)
    #  21
    a =  - elp + f2 - d2 + om
    dp = dp - 5 * math.sin(a)
    de = de + 3 * math.cos(a)
    #  20
    a =  - el2 + d2 + om
    dp = dp - 6 * math.sin(a)
    de = de + 3 * math.cos(a)
    #  19
    a =  - elp + om
    dp = dp - 12 * math.sin(a)
    de = de + 6 * math.cos(a)
    #  18
    a = elp2 + f2 - d2 + om2
    dp = dp + (-16 + 0.1 * t) * math.sin(a)
    de = de + 7 * math.cos(a)
    #  17
    a = elp + om
    dp = dp - 15 * math.sin(a)
    de = de + 9 * math.cos(a)
    #  16
    dp = dp + (17 - (0.1 * t)) * math.sin(elp2)
    #  15
    dp = dp - 22 * math.sin(f2 - d2)
    #  14
    a = el2 - d2
    dp = dp + 48 * math.sin(a)
    de = de + math.cos(a)
    #  13
    a = f2 - d2 + om
    dp = dp + (129 + 0.1 * t) * math.sin(a)
    de = de - 70 * math.cos(a)
    #  12
    a =  - elp + f2 - d2 + om2
    dp = dp + (217 - 0.5 * t) * math.sin(a)
    de = de + (-95 + 0.3 * t) * math.cos(a)
    #  11
    a = elp + f2 - d2 + om2
    dp = dp + (-517 + 1.2 * t) * math.sin(a)
    de = de + (224 - 0.6 * t) * math.cos(a)
    #  10
    dp = dp + (1426 - 3.4 * t) * math.sin(elp)
    de = de + (54 - 0.1 * t) * math.cos(elp)
    #  9
    a = f2 - d2 + om2
    dp = dp + (-13187 - 1.6 * t) * math.sin(a)
    de = de + (5736 - 3.1 * t) * math.cos(a)
    #  8
    dp = dp + math.sin(el2 - f2 + om)
    #  7
    a =  - elp2 + f2 - d2 + om
    dp = dp - 2 * math.sin(a)
    de = de + 1 * math.cos(a)
    #  6
    dp = dp - 3 * math.sin(el - elp - d)
    #  5
    a =  - el2 + f2 + om2
    dp = dp - 3 * math.sin(a)
    de = de + 1 * math.cos(a)
    #  4
    dp = dp + 11 * math.sin(el2 - f2)
    #  3
    a =  - el2 + f2 + om
    dp = dp + 46 * math.sin(a)
    de = de - 24 * math.cos(a)
    #  2
    dp = dp + (2062 + 0.2 * t) * math.sin(om2)
    de = de + (-895 + 0.5 * t) * math.cos(om2)
    #  1
    dp = dp + (-171996 - 174.2 * t) * math.sin(om)
    de = de + (92025 + 8.9 * t) * math.cos(om)

    #  Convert results to radians
    dpsi = dp * u2r
    deps = de * u2r

    #  Mean obliquity
    eps0 = RadPerArcSec * (84381.448 +  \
        (-46.8150 + (-0.00059 + (0.001813 * t)) * t) * t)

    return (dpsi, deps, eps0)

def gmstFromUT1(ut1):
    """Convert from universal time (MJD)
    to Greenwich mean sidereal time, in degrees

    Based on Pat Wallace's GMST, whose comments follow:
    
    The IAU 1982 expression (see page S15 of 1984 Astronomical
    Almanac) is used, but rearranged to reduce rounding errors.
    This expression is always described as giving the GMST at
    0 hours UT.  In fact, it gives the difference between the
    GMST and the UT, which happens to equal the GMST (modulo
    24 hours) at 0 hours UT each day.  In this routine, the
    entire UT is used directly as the argument for the
    standard formula, and the fractional part of the UT is
    added separately;  note that the factor 1.0027379... does
    not appear.
    
    See also the routine GMSTA, which delivers better numerical
    precision by accepting the UT date and time as separate arguments.
    
    P.T.Wallace   Starlink   14 September 1995
    
    History:
    2002-20-12 ROwen    Removed an extra + sign that was doing nothing.
                        Thanks to pychecker.
    """
    # convert date to Julian centuries since J2000
    jc= (ut1 - coordConv.MJDJ2000) / 36525.0

    return coordConv.wrapPos (
        (math.fmod(ut1, 1.0) * 360.0) # fraction of day of UT1, in degrees
         + (24110.54841
            + (8640184.812866
               + (0.093104-(6.2e-6*jc))*jc)*jc)*0.0041666666666666666)

def lastFromTAI(tai, site):
    ut1Days = (tai + site.ut1_tai) / coordConv.SecPerDay
    gmst = gmstFromUT1(ut1Days)

    # compute apparent - mean sidereal time, in degrees
    ttDays = (tai + coordConv.TT_TAI) / coordConv.SecPerDay;
    appMinusMean = eqeqx(ttDays) / coordConv.RadPerDeg;

    # compute local apparent sideral time, in degrees, in range [0, 360)
    return coordConv.wrapPos(gmst + site.corrLong + appMinusMean);

class TestTime(unittest.TestCase):
    """Test time functions
    """
    def testJulianEpochFromTAI(self):
        """Test julianEpochFromTAI and inverse
        """
        for tai in (4232.89, 20000.32, 56350.03, 74222.9):
            self.assertAlmostEqual(coordConv.julianEpochFromTAI(tai), julianEpochFromTAI(tai), places=7)
            
            self.assertAlmostEqual(coordConv.taiFromJulianEpoch(coordConv.julianEpochFromTAI(tai)), tai, places=5)

    def testBesselianEpochFromTAI(self):
        """Test besselianEpochFromTAI and inverse
        """
        for tai in (4232.89, 20000.32, 56350.03, 74222.9):
            self.assertAlmostEqual(coordConv.besselianEpochFromTAI(tai), besselianEpochFromTAI(tai), places=7)
            
            self.assertAlmostEqual(coordConv.taiFromBesselianEpoch(coordConv.besselianEpochFromTAI(tai)), tai, places=5)
    
    def testLastFromTAI(self):
        """Test lastFromTAI
        
        This is not a very accurate test because my reference nutc function above is too old.
        """
        site = coordConv.Site(10, 10, 10) # values are not used
        for ut1_tai in (0, -10.5, 23.3, 900.3):
            site.ut1_tai = ut1_tai
            for tai in (4232.89, 20000.32, 56350.03, 74222.9):
                self.assertAlmostEqual(coordConv.lastFromTAI(tai, site), lastFromTAI(tai, site), places=3)

if __name__ == '__main__':
    unittest.main()
