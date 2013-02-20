#!/usr/bin/env python
import unittest
import math
import numpy
from coordConv import Coord, sind, cosd, tand, atan2d, wrapCtr, wrapPos, \
    DoubleEpsilon, MinParallax, AUPerParsec, RadPerDeg, ArcsecPerDeg, SecPerDay, DaysPerYear, KmPerAU

class TestCoord(unittest.TestCase):
    def testConstants(self):
        """Some of these tests rely on the physical constants being correct, so verify that
        """
        self.assertAlmostEqual(AUPerParsec, (180 * 3600) / math.pi)
        self.assertAlmostEqual(RadPerDeg, math.pi / 180.0)
        self.assertAlmostEqual(ArcsecPerDeg, 60 * 60)
        self.assertAlmostEqual(SecPerDay, 24 * 60 * 60)
        self.assertAlmostEqual(DaysPerYear, 365.25)
        self.assertAlmostEqual(KmPerAU, 149597871)
        
    def testBasics(self):
        """Check sph -> vec and back again for points not at the pole
        
        The round trip test relies on the fact that the internal representation of a coord is vector,
        so simply retrieving the sph info again suffices to test a round trip.
        
        Warning: the test of vector velocity is poor because the code is copied from the C++.
        However, two other tests will help:
        - Convert a coordinate system to itself with two different dates of observation
          and check that the expected amount of spatial motion occurs
        - Compare coordinate conversion results to the old TCC
        """
        for equatAng in (0, 71, -123.4):
            for polarAng in (0, -75, -89.999999, 89.999999):
                for parallax in (0, MinParallax / 0.9000001, MinParallax / 0.8999999, 5):
                    for equatPM in (4,): # (0, 4):
                        for polarPM in (-7,): # (0, -7):
                            for radVel in (0, -34):
                                coord = Coord(equatAng, polarAng, parallax, equatPM, polarPM, radVel)
                                self.assertFalse(coord.atPole())
                                self.assertEqual(coord.atInfinity(), parallax < MinParallax / 0.9)

                                # check vector position
                                vecPos = coord.getVecPos()
                                dist = coord.getDist()
                                self.assertAlmostEqual(numpy.linalg.norm(vecPos), dist, 2)
                                predVecPos = (
                                    dist * cosd(polarAng) * cosd(equatAng),
                                    dist * cosd(polarAng) * sind(equatAng),
                                    dist * sind(polarAng),
                                )
                                self.assertTrue(numpy.allclose(predVecPos, vecPos))
                                
                                # check vector velocity
                                vecVel = coord.getVecVel()
                                RadPerYear_per_ArcsecPerCentury = RadPerDeg / (ArcsecPerDeg * 100.0)
                                SecPerYear = SecPerDay * DaysPerYear
                                AUPerYear_per_KmPerSec = SecPerYear / KmPerAU

                                sinEquat = sind(equatAng)
                                cosEquat = cosd(equatAng)
                                sinPolar = sind(polarAng)
                                cosPolar = cosd(polarAng)

                                # change units of proper motion from arcsec/century to au/year
                                pmAUPerYear1 = equatPM * dist * RadPerYear_per_ArcsecPerCentury
                                pmAUPerYear2 = polarPM * dist * RadPerYear_per_ArcsecPerCentury

                                # change units of radial velocity from km/sec to au/year
                                radVelAUPerYear = radVel * AUPerYear_per_KmPerSec
                                
                                predVecVel = (
                                    - (pmAUPerYear2 * sinPolar * cosEquat) - (pmAUPerYear1 * cosPolar * sinEquat) + (radVelAUPerYear * cosPolar * cosEquat),
                                    - (pmAUPerYear2 * sinPolar * sinEquat) + (pmAUPerYear1 * cosPolar * cosEquat) + (radVelAUPerYear * cosPolar * sinEquat),
                                    + (pmAUPerYear2 * cosPolar)                                                   + (radVelAUPerYear * sinPolar),
                                )
                                self.assertTrue(numpy.allclose(predVecVel, vecVel))

                                # check round trip
                                atPole, destEquatAng, destPolarAng = coord.getSphPos()
                                self.assertAlmostEqual(wrapPos(equatAng), destEquatAng)
                                self.assertAlmostEqual(polarAng, destPolarAng)

                                destParallax = coord.getParallax()
                                predParallax = 0 if coord.atInfinity() else parallax
                                self.assertAlmostEqual(predParallax, destParallax)

                                atPole, destEquatPM, destPolarPM = coord.getPM()
                                self.assertAlmostEqual(equatPM, destEquatPM, 6)
                                self.assertAlmostEqual(polarPM, destPolarPM, 6)

                                destRadVel = coord.getRadVel()
                                self.assertAlmostEqual(radVel, destRadVel)

        
    def testAtPole(self):
        """Test atPole computation
        """
        for equatAng in (0, 33.6, -123.4):
            for polarAng in (0, -89.999, 89.999, 89.9999999, -89.9999999):
                coord = Coord(equatAng, polarAng)
                vec = coord.getVecPos()
                fracXYMag = math.hypot(vec[0], vec[1]) / coord.getDist()
                predAtPole = fracXYMag**2 < DoubleEpsilon
                self.assertEqual(predAtPole, coord.atPole())
    
    def testDist(self):
        """Test distance
        """
        for parallax in (0, MinParallax / 0.899999999, MinParallax / 0.900000001, MinParallax, 1, 3.4, 75.3):
            adjParallax = max(parallax, MinParallax)
            predDist = AUPerParsec / adjParallax
            predAtInf = parallax < MinParallax / 0.9
            for coord in ( # test one with space motion, one without
                Coord(43, 23, parallax),
                Coord(-32, 89.99, parallax, 3, 5, 2),
            ):
                self.assertAlmostEqual(coord.getDist(), predDist, 2)
                self.assertEqual(predAtInf, coord.atInfinity())


if __name__ == '__main__':
    unittest.main()
