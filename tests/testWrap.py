#!/usr/bin/env python
import unittest
import numpy
from coordConv import wrapCtr, wrapPos, wrapNear, sind, cosd, PVT

Eps = 2e-14

class TestWrap(unittest.TestCase):
    def testWrapPos(self):
        for wrap in (-1000, -10, -1, 0, 1, 10, 1000):
            for offset in (-360, -180, -90, 0, 90, 180, 270, 360):
                for epsMult in (-10, -4, -2, -1, 0, 1, 2, 4, 10):
                    ang = (wrap * 360) + offset + (epsMult * Eps)
                    pvt = PVT(ang, ang, 35.0) # pick anything for vel and time

                    posAng = wrapPos(ang)
                    self.assertGreaterEqual(posAng, 0.0)
                    self.assertLess(posAng, 360.0)
                    # prove that posAng and ang are the same angle
                    self.assertAlmostEqual(wrapCtr(posAng - ang), 0)
                    self.assertAlmostEqual(sind(posAng), sind(ang))
                    self.assertAlmostEqual(cosd(posAng), cosd(ang))

                    posPvt = wrapPos(pvt)
                    self.assertEqual(posPvt.pos, posAng)
                    self.assertEqual(posPvt.vel, pvt.vel)
                    self.assertEqual(posPvt.t, pvt.t)

                    ctrAng = wrapCtr(ang)
                    self.assertGreaterEqual(ctrAng, -180.0)
                    self.assertLess(ctrAng, 180.0)
                    # prove that ctrAng and ang are the same angle
                    self.assertAlmostEqual(wrapCtr(ctrAng - ang), 0)
                    self.assertAlmostEqual(sind(ctrAng), sind(ang))
                    self.assertAlmostEqual(cosd(ctrAng), cosd(ang))

                    ctrPvt = wrapCtr(pvt)
                    self.assertEqual(ctrPvt.pos, ctrAng)
                    self.assertEqual(ctrPvt.vel, pvt.vel)
                    self.assertEqual(ctrPvt.t, pvt.t)
                    
                    for otherAngBase in (-180, 0, 180, 360):
                        for otherEpsMult in (-10, -4, -2, -1, 0, 1, 2, 4, 10):
                            otherAng = otherAngBase + (otherEpsMult * Eps)
                            nearAng = wrapNear(ang, otherAng)
                            self.assertGreaterEqual(nearAng - otherAng, -180)
                            self.assertLess(nearAng - otherAng, 180)
                            # prove that nearAng and ang are the same angle
                            self.assertAlmostEqual(wrapCtr(nearAng - ang), 0)
                            self.assertAlmostEqual(sind(nearAng), sind(ang))
                            self.assertAlmostEqual(cosd(nearAng), cosd(ang))
                            

if __name__ == '__main__':
    unittest.main()
