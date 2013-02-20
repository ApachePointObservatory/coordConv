#!/usr/bin/env python
import unittest
from coordConv import offsetSph, angSideAng, Coord, cosd, sind

Eps = 2e-16
EpsTest = Eps * 1.001

class TestOffsetSph(unittest.TestCase):
    """Test offsetSph and Coord.angularSeparation
    """
    def testSmallOffset(self):
        """Test small offsets not too near the pole
        
        In this regime delta-long = dist along long / cos(lat) is a reasonable approximation
        but I have no simple way to compute destAng
        """
        for srcLat in (-40.0, 0.43, 36.7):
            cosLat = cosd(srcLat)
            for srcLong in (0, 41.0, -92.0): # should not matter
                for srcAng in (-45.0, 0.01, 12.5):
                    cosSrcAng = cosd(srcAng)
                    sinSrcAng = sind(srcAng)
                    for dist in (0.01, 0.13):
                        predLong = srcLong + (dist * cosSrcAng / cosLat)
                        predLat = srcLat + (dist * sinSrcAng)
                        destLong, destLat, destAng = offsetSph(srcLong, srcLat, srcAng, dist)
                        if abs(dist) > 0.1:
                            places = 3
                        else:
                            places = 5
                        self.assertAlmostEqual(destLong, predLong, places=places)
                        self.assertAlmostEqual(destLat, predLat, places=places)
                        srcCoord = Coord(srcLong, srcLat)
                        destCoord = Coord(destLong, destLat)
                        self.assertAlmostEqual(dist, srcCoord.angularSeparation(destCoord))

    def testWideRange(self):
        """Test over a wide range of angles
        """
        for srcLat in (-87.1, -25.5, 0.43, 36.7, 87.0):
            for srcLong in (0, 41.0,-92.0): # should not matter
                for srcAng in (-89.9, -45.0, 0.01, 12.5, 89.0, 90.0):
                    for dist in (0.01, 0.13, 5.73):
                        sideA = 90.0 - srcLat
                        angB = 90.0 - srcAng
                        sideC = dist
                        unknownAng, angA, sideB, angC = angSideAng(sideA, angB, sideC)
                        self.assertFalse(unknownAng)
                        predLong = srcLong + angC
                        predLat = 90 - sideB
                        predAng = angA - 90
                        destLong, destLat, destAng = offsetSph(srcLong, srcLat, srcAng, dist)
                        places = 7
                        self.assertAlmostEqual(destLong, predLong, places=places)
                        self.assertAlmostEqual(destLat, predLat, places=places)
                        self.assertAlmostEqual(destAng, predAng, places=places)
                        srcCoord = Coord(srcLong, srcLat)
                        destCoord = Coord(destLong, destLat)
                        self.assertAlmostEqual(dist, srcCoord.angularSeparation(destCoord))
        
                        
def processOutput(outputVec):
    return (
        outputVec[0],
        sind(outputVec[1]), cosd(outputVec[2]),
        outputVec[2],
        sind(outputVec[3]), cosd(outputVec[3]),
    )

if __name__ == '__main__':
    unittest.main()
