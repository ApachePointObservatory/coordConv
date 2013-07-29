#!/usr/bin/env python
import unittest
import math
import numpy
from coordConv import Coord, sind, cosd, wrapPos, wrapNear, angSideAng, \
    DoubleEpsilon, MinParallax, AUPerParsec, RadPerDeg, ArcsecPerDeg, SecPerDay, DaysPerYear, KmPerAU

class TestCoord(unittest.TestCase):
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
                                
                                # check vector proper motion
                                vecPM = coord.getVecPM()
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
                                
                                predVecPM = (
                                    - (pmAUPerYear2 * sinPolar * cosEquat) - (pmAUPerYear1 * cosPolar * sinEquat) + (radVelAUPerYear * cosPolar * cosEquat),
                                    - (pmAUPerYear2 * sinPolar * sinEquat) + (pmAUPerYear1 * cosPolar * cosEquat) + (radVelAUPerYear * cosPolar * sinEquat),
                                    + (pmAUPerYear2 * cosPolar)                                                   + (radVelAUPerYear * sinPolar),
                                )
                                self.assertTrue(numpy.allclose(predVecPM, vecPM))

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
                                
                                coordCopy = Coord(coord)
                                self.assertEqual(repr(coord), repr(coordCopy))
        
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

    def testOffsetSmall(self):
        """Test offset, angularSeparation and orientationTo for small offsets not too near the pole
        
        In this regime delta-long = dist along long / cos(lat) is a reasonable approximation
        but I have no simple way to compute toOrient, except if dist very small then toOrient = fromOrient
        """
        for fromPolarAng in (-40.0, 0.43, 36.7):
            cosPolarAng = cosd(fromPolarAng)
            for fromEquatAng in (0, 41.0): # should not matter
                fromCoord = Coord(fromEquatAng, fromPolarAng)
                for fromOrient in (-90, -72, -45.0, -30, 0.01, 12.5, 31, 47, 56, 68, 89):
                    cosFromOrient = cosd(fromOrient)
                    sinFromOrient = sind(fromOrient)
                    for dist in (0, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-5, 0.001, 0.01, 0.13):
                        predEquatAng = fromEquatAng + (dist * cosFromOrient / cosPolarAng)
                        predPolarAng = fromPolarAng + (dist * sinFromOrient)

                        toCoord, toOrient = fromCoord.offset(fromOrient, dist)
                        atPole, toEquatAng, toPolarAng = toCoord.getSphPos()
                        if abs(dist) > 0.1:
                            places = 3
                        else:
                            places = 5
                        if abs(dist) < 0.0001:
                            orientPlaces = 4
                        else:
                            orientPlaces = 7
                        self.assertAlmostEqual(toEquatAng, predEquatAng, places=places)
                        self.assertAlmostEqual(toPolarAng, predPolarAng, places=places)
                        fromCoord = Coord(fromEquatAng, fromPolarAng)
                        toCoord = Coord(toEquatAng, toPolarAng)
                        self.assertAlmostEqual(dist, fromCoord.angularSeparation(toCoord))
                        predFromOrient = fromCoord.orientationTo(toCoord)
                        if numpy.isfinite(predFromOrient):
                            self.assertAlmostEqual(fromOrient, predFromOrient, places=orientPlaces)
                            self.assertAlmostEqual(toOrient, wrapNear(180 + toCoord.orientationTo(fromCoord), toOrient), places=orientPlaces)
                        else:
                            self.assertLess(dist, 1e-7)
                            self.assertAlmostEqual(fromEquatAng, toEquatAng)
                            self.assertAlmostEqual(fromPolarAng, toPolarAng)
                            self.assertAlmostEqual(fromOrient, toOrient)

    def testOffset(self):
        """Test offset, angularSeparation and orientationTo for offsets over a wide range of angles
        """
        for fromPolarAng in (-87.1, -25.5, 0.43, 36.7, 87.0):
            for fromEquatAng in (0, 41.0): # should not matter
                fromCoord = Coord(fromEquatAng, fromPolarAng)
                for fromOrient in (-89.9, -45.0, 0.01, 12.5, 89.0, 90.0):
                    for dist in (0.0001, 0.01, 0.13, 5.73):
                        sideA = 90.0 - fromPolarAng
                        angB = 90.0 - fromOrient
                        sideC = dist
                        unknownAng, angA, sideB, angC = angSideAng(sideA, angB, sideC)
                        self.assertFalse(unknownAng)
                        predEquatAng = fromEquatAng + angC
                        predPolarAng = 90 - sideB
                        predAng = angA - 90
                        toCoord, toOrient = fromCoord.offset(fromOrient, dist)
                        atPole, toEquatAng, toPolarAng = toCoord.getSphPos()
                        places = 7
                        self.assertAlmostEqual(toEquatAng, predEquatAng, places=places)
                        self.assertAlmostEqual(toPolarAng, predPolarAng, places=places)
                        self.assertAlmostEqual(toOrient, predAng, places=places)
                        fromCoord = Coord(fromEquatAng, fromPolarAng)
                        toCoord = Coord(toEquatAng, toPolarAng)
                        self.assertAlmostEqual(dist, fromCoord.angularSeparation(toCoord))
                        self.assertAlmostEqual(fromOrient, fromCoord.orientationTo(toCoord))
                        self.assertAlmostEqual(toOrient, wrapNear(180 + toCoord.orientationTo(fromCoord), toOrient))


if __name__ == '__main__':
    unittest.main()
