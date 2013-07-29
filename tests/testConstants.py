#!/usr/bin/env python
import time
import datetime
import unittest
import math
import coordConv

class TestConstants(unittest.TestCase):
    def testConstants(self):
        """Test various constants
        """
        self.assertAlmostEqual(coordConv.Pi, math.pi)
        self.assertAlmostEqual(coordConv.HoursPerDeg, 24.0 / 360.0)
        self.assertAlmostEqual(coordConv.RadPerDeg, math.pi / 180.0)
        self.assertAlmostEqual(coordConv.KmPerAU, 149597871)
        self.assertAlmostEqual(coordConv.AUPerParsec, (180 * 3600) / math.pi)
        self.assertAlmostEqual(coordConv.ArcsecPerDeg, 60 * 60)
        self.assertAlmostEqual(coordConv.SecPerDay, 24 * 60 * 60)
        self.assertAlmostEqual(coordConv.DaysPerYear, 365.25)
        self.assertAlmostEqual(coordConv.VLight, 299792.458 / coordConv.KmPerAU) # in au/sec; based on c=299792458 m/sec
        self.assertAlmostEqual(coordConv.DegK_DegC, 273.15)
        unixZero = datetime.datetime(1970, 1, 1, 0, 0, 0)
        mjdZero = datetime.datetime(1858, 11, 17, 0, 0, 0)
        self.assertAlmostEqual(coordConv.MJD_UnixTime, (unixZero - mjdZero).total_seconds())
        self.assertAlmostEqual(coordConv.MJDJ2000, 51544.5)
        self.assertAlmostEqual(coordConv.AngstromsPerMicron, 1.0e4)
        self.assertAlmostEqual(coordConv.PascalsPerMillibar, 100.0)

if __name__ == '__main__':
    unittest.main()
