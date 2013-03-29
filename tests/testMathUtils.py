#!/usr/bin/env python
import unittest
import math
import numpy
import coordConv

class TestMathUtils(unittest.TestCase):
    """Test mathUtils except for the wrap functions (which are tested elsewhere)
    """
    def testRot2D(self):
        """Test rot2D using various obvious values
        """
        invSqrt2 = 1.0 / math.sqrt(2)
        for inArgs, predOut in (
            ((1, 0, 0), (1, 0)),
            ((1, 0, 90), (0, 1)),
            ((1, 0, -90), (0, -1)),
            ((1, 0, 180), (-1, 0)),
            ((1, 0, -180), (-1, 0)),
            ((1, 0, 45), (invSqrt2, invSqrt2)),
            ((1, 0, -45), (invSqrt2, -invSqrt2)),
            ((1.2345e67, 0, -45), (1.2345e67 * invSqrt2, -1.2345e67 * invSqrt2)),
            ((0, 1, 0), (0, 1)),
            ((0, 1, 90), (-1, 0)),
            ((0, 1, -90), (1, 0)),
            ((0, 1, 180), (0, -1)),
            ((0, 1, -180), (0, -1)),
            ((0, 1, 45), (-invSqrt2, invSqrt2)),
            ((0, 1, -45), (invSqrt2, invSqrt2)),
            ((0, 1.2345e67, -45), (1.2345e67 * invSqrt2, 1.2345e67 * invSqrt2)),
        ):
            out = coordConv.rot2D(*inArgs)
            self.assertTrue(numpy.allclose(out, predOut))
    
    def testHypot(self):
        """Test hypot
        """
        for x in (0, 1e-45, 12.3, 3234.34):
            for y in (0, -2e-45, -3.3, 234.32):
                self.assertTrue(numpy.allclose([coordConv.hypot(x, y)], [math.hypot(x, y)]))
    
    def testIsFinite(self):
        """Test isFinite
        """
        for x, isFinite in (
            (0, True),
            (coordConv.DoubleMax, True),
            (-coordConv.DoubleMax, True),
            (numpy.nan, False),
            (numpy.inf, False), 
            (-numpy.inf, False),
        ):
            self.assertEqual(coordConv.isFinite(x), isFinite)
    
    def testTrig(self):
        """Test degrees-based trig functions
        """
        for ang in (0, -1, 31.2, -235, 234324):
            angRad = ang * coordConv.RadPerDeg
            self.assertAlmostEqual(coordConv.sind(ang), math.sin(angRad))
            self.assertAlmostEqual(coordConv.cosd(ang), math.cos(angRad))
            self.assertAlmostEqual(coordConv.tand(ang), math.tan(angRad))
            self.assertAlmostEqual(coordConv.atand(ang) * coordConv.RadPerDeg, math.atan(ang))
            for ang2 in (0, 35, -364):
                self.assertAlmostEqual(coordConv.atan2d(ang, ang2) * coordConv.RadPerDeg, math.atan2(ang, ang2))
        
        for val in (-1, -0.34, 0, 0.67, 1):
            self.assertAlmostEqual(coordConv.asind(val) * coordConv.RadPerDeg, math.asin(val))
            self.assertAlmostEqual(coordConv.acosd(val) * coordConv.RadPerDeg, math.acos(val))
    
    def testPolarFromXY(self):
        """Test polarFromXY and xyFromPolar
        """
        sqrt2 = math.sqrt(2.0)
        for xy, predPol in (
            ((1, 0), (False, 1, 0)),
            ((-1, 0), (False, 1, 180)),
            ((0, 1), (False, 1, 90)),
            ((0, -1), (False, 1, -90)),
            ((1, 1), (False, sqrt2, 45)),
            ((1, -1), (False, sqrt2, -45)),
            ((-1, 1), (False, sqrt2, 135)),
            ((-1, -1), (False, sqrt2, -135)),
            ((-123.45, -123.45), (False, 123.45 * sqrt2, -135)),
            ((0, 0), (True, 0, 0)),
        ):
            pol = coordConv.polarFromXY(*xy)
            self.assertEqual(pol[0], predPol[0])
            self.assertTrue(numpy.allclose(pol[1:], predPol[1:]))
            compXY = coordConv.xyFromPolar(*pol[1:])
            self.assertTrue(numpy.allclose(xy, compXY))

if __name__ == '__main__':
    unittest.main()
