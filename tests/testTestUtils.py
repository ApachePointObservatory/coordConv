#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

import unittest

import coordConv


class TestTestUtils(unittest.TestCase):
    """Test the contents of testUtils
    """

    def testAnglesAlmostEquals(self):
        for places in (0, 3, 7):
            maxErr = float(10**-places)
            for angle1 in (-9000, -27.6, 271.3):
                for nWrap in (0, -5, 1, 4):
                    wrappedAngle1 = angle1 + (360.0 * nWrap)
                    eps = wrappedAngle1 * 1e-14
                    for deltaAngle in (maxErr - eps, maxErr + eps):
                        angle2 = wrappedAngle1 + deltaAngle
                        if abs(deltaAngle) <= maxErr:
                            coordConv.assertAnglesAlmostEqual(angle1, angle2, places=places)
                        else:
                            with self.assertRaises(AssertionError):
                                coordConv.assertAnglesAlmostEqual(angle1, angle2, places=places)

    def testPVTsAlmostEqual(self):
        for (posPlaces, velPlaces, tPlaces) in (
            (5, 6, 7),
            (6, 7, 5),
        ):
            maxPosErr = 10**-posPlaces
            maxVelErr = 10**-velPlaces
            maxTErr = 10**-tPlaces
            for pvt1 in (
                coordConv.PVT(5, 0.2, 3543),
                coordConv.PVT(386, -230.3, 5923402.22),
            ):
                for nWrap in (0, -3, 1):
                    if nWrap == 0:
                        doWrapList = (False, True)
                    else:
                        doWrapList = (True,)

                    wrappedPVT = pvt1.copy()
                    wrappedPVT.pos += 360.0 * nWrap
                    posEps = wrappedPVT.pos * 1e-14
                    velEps = pvt1.vel * 1e-14
                    tEps = pvt1.t * 1e-14

                    for deltaPos in (maxPosErr - posEps, maxPosErr + posEps):
                        for deltaVel in (maxVelErr - velEps, maxVelErr + velEps):
                            for deltaT in (maxTErr - tEps, maxTErr + tEps):
                                pvt2 = wrappedPVT.copy()
                                pvt2.pos += deltaPos
                                pvt2.vel += deltaVel
                                pvt2.t += deltaT
                                for doWrap in doWrapList:
                                    if abs(deltaPos) <= maxPosErr \
                                            and abs(deltaVel) <= maxVelErr \
                                            and abs(deltaT) <= maxTErr:
                                        coordConv.assertPVTsAlmostEqual(
                                            pvt1, pvt2,
                                            doWrap = doWrap,
                                            posPlaces = posPlaces,
                                            velPlaces = velPlaces,
                                            tPlaces = tPlaces,
                                        )
                                    else:
                                        with self.assertRaises(AssertionError):
                                            coordConv.assertPVTsAlmostEqual(
                                                pvt1, pvt2,
                                                doWrap = doWrap,
                                                posPlaces = posPlaces,
                                                velPlaces = velPlaces,
                                                tPlaces = tPlaces,
                                            )


if __name__ == '__main__':
    unittest.main()
