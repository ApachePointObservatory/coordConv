#!/usr/bin/env python
from __future__ import absolute_import # division breaks PVT / float

import unittest
import numpy
import coordConv

DeltaT = 0.01

def refGetPos(pvt, t):
    """Reference implementation of PVT.getPos
    """
    return pvt.pos + (pvt.vel * (t - pvt.t))

def refPolarFromXY(x, y, tai):
    """Reference implementation of PVT.polarFromXY
    """
    atPole = False
    rArr = []
    thetaArr = []
    for testTAI in (tai, tai + DeltaT):
        ap, r, theta = coordConv.polarFromXY(x.getPos(testTAI), y.getPos(testTAI))
        rArr.append(r)
        thetaArr.append(theta)
        atPole = atPole or ap
    rPVT = coordConv.PVT()
    rPVT.pos = rArr[0]
    rPVT.vel = (rArr[1] - rArr[0]) / DeltaT
    rPVT.t = tai
    thetaPVT = coordConv.PVT()
    thetaPVT.pos = thetaArr[0]
    thetaPVT.vel = coordConv.wrapCtr(thetaArr[1] - thetaArr[0]) / DeltaT
    thetaPVT.t = tai
    return atPole, rPVT, thetaPVT

class TestPVT(unittest.TestCase):
    def testConstructors(self):
        """Test PVT constructors
        """
        pvt = coordConv.PVT(1.0, 2.0, 3.0)
        self.assertEquals(pvt.pos, 1.0)
        self.assertEquals(pvt.vel, 2.0)
        self.assertEquals(pvt.t, 3.0)
        
        pvtCopy = coordConv.PVT(pvt)
        self.assertEquals(pvt.pos, pvtCopy.pos)
        self.assertEquals(pvt.vel, pvtCopy.vel)
        self.assertEquals(pvt.t,   pvtCopy.t)
    
    def testCopy(self):
        """Test PVT.copy() and PVT.copy(t)
        """
        pvt1 = coordConv.PVT(1.0, -2.0, 3.0)
        pvt2 = pvt1.copy()
        pvt3 = pvt1.copy(5.0)

        pvt1 *= 2 # modify pvt1 and make sure pvt2 and pvt3 are not affected
        self.assertAlmostEqual(pvt1.pos, 2.0)
        self.assertAlmostEqual(pvt1.vel, -4.0)
        self.assertEqual(pvt1.t, 3.0)

        self.assertEqual(pvt2.pos, 1.0)
        self.assertEqual(pvt2.vel, -2.0)
        self.assertEqual(pvt2.t, 3.0)

        self.assertEqual(pvt3.pos, -3.0)
        self.assertEqual(pvt3.vel, -2.0)
        self.assertEqual(pvt3.t, 5.0)

    def testEquality(self):
        """Test operator== and operator!=
        """
        def pvtIter():
            for pos in (-1.1, 0, 1.1, numpy.nan):
                for vel in (-0.1, 0, 0.1, numpy.nan):
                    for t in (0, 1, 2, numpy.nan):
                        yield coordConv.PVT(pos, vel, t)

        for pvt1 in pvtIter():
            for pvt2 in pvtIter():
                if pvt1.pos == pvt2.pos and pvt1.vel == pvt2.vel and pvt1.t == pvt2.t:
                    self.assertTrue(pvt1 == pvt2)
                else:
                    self.assertTrue(pvt1 != pvt2)
                self.assertNotEqual(pvt1 == pvt2, pvt1 != pvt2)
    
    def testAddPVT(self):
        """Test pvt + pvt
        """
        pvt1 = coordConv.PVT(1.0, 2.0, 3.0)
        pvt2 = coordConv.PVT(-1.5, -3.0, 4.0)
        pvt3 = pvt1 + pvt2
        predPos3 = pvt1.pos + refGetPos(pvt2, pvt1.t)
        self.assertAlmostEqual(pvt3.pos, predPos3)
        self.assertAlmostEqual(pvt3.vel, -1.0)
        self.assertAlmostEqual(pvt3.t, 3.0)

    
    def testSubtractPVT(self):
        """Test pvt - pvt
        """
        pvt1 = coordConv.PVT(1.0, 2.0, 3.0)
        pvt2 = coordConv.PVT(-1.5, -3.0, 4.0)
        pvt3 = pvt1 - pvt2
        predPos3 = pvt1.pos - refGetPos(pvt2, pvt1.t)
        self.assertAlmostEqual(pvt3.pos, predPos3)
        self.assertAlmostEqual(pvt3.vel, 5.0)
        self.assertAlmostEqual(pvt3.t, 3.0)

        pvt1 -= pvt2
        self.assertAlmostEqual(pvt1.pos, predPos3)
        self.assertAlmostEqual(pvt1.vel, 5.0)
        self.assertAlmostEqual(pvt1.t, 3.0)
    
    def testUnaryMinus(self):
        """Test -pvt
        """
        for pvt in (
            coordConv.PVT(1.0, 2.0, 3.0),
            coordConv.PVT(-1.5, -3.0, 4.0)
        ):
            negPVT = -pvt
            self.assertEqual(negPVT.t, pvt.t)
            self.assertEqual(negPVT.pos, -pvt.pos)
            self.assertEqual(negPVT.vel, -pvt.vel)
    
    def testAddSubtractScalar(self):
        """Test pvt + scalar, pvt += scalar, pvt - scalar and pvt -= scalar
        """
        for pvt in (
            coordConv.PVT(1.0, 2.0, 3.0),
            coordConv.PVT(-1.5, -3.0, 4.0)
        ):
            for val in (-1.234e97, -25.3, -1.0e-13, 0.0, 1.0e13, 39.6, 1.34e99):
                addPVT = pvt + val
                self.assertAlmostEqual(addPVT.pos, pvt.pos + val)
                self.assertEqual(addPVT.vel, pvt.vel)
                self.assertEqual(addPVT.t, pvt.t)
                
                inPlaceAddPVT = coordConv.PVT(pvt.pos, pvt.vel, pvt.t)
                inPlaceAddPVT += val
                self.assertEqual(inPlaceAddPVT.pos, addPVT.pos)
                self.assertEqual(inPlaceAddPVT.vel, addPVT.vel)
                self.assertEqual(inPlaceAddPVT.t, addPVT.t)

                subPVT = pvt - val
                self.assertAlmostEqual(subPVT.pos, pvt.pos - val)
                self.assertEqual(subPVT.vel, pvt.vel)
                self.assertEqual(subPVT.t, pvt.t)

                inPlaceSubPVT = coordConv.PVT(pvt.pos, pvt.vel, pvt.t)
                inPlaceSubPVT -= val
                self.assertEqual(inPlaceSubPVT.pos, subPVT.pos)
                self.assertEqual(inPlaceSubPVT.vel, subPVT.vel)
                self.assertEqual(inPlaceSubPVT.t, subPVT.t)
    
    def testMultDivScalar(self):
        """Test pvt * scalar, pvt *= scalar, pvt / scalar and pvt /= scalar
        """
        for pvt in (
            coordConv.PVT(1.0, 2.0, 3.0),
            coordConv.PVT(-1.5, -3.0, 4.0)
        ):
            for val in (-1.234e97, -25.3, -1.0e-13, 0.0, 1.0e13, 39.6, 1.34e99):
                multPVT = pvt * val
                self.assertAlmostEqual(multPVT.pos, pvt.pos * val)
                self.assertEqual(multPVT.vel, pvt.vel * val)
                self.assertEqual(multPVT.t, pvt.t)

                inPlaceMultPVT = coordConv.PVT(pvt.pos, pvt.vel, pvt.t)
                inPlaceMultPVT *= val
                self.assertEqual(inPlaceMultPVT.pos, multPVT.pos)
                self.assertEqual(inPlaceMultPVT.vel, multPVT.vel)
                self.assertEqual(inPlaceMultPVT.t, multPVT.t)

                if val != 0.0:
                    divPVT = pvt / val
                    self.assertAlmostEqual(divPVT.pos, pvt.pos / val)
                    self.assertEqual(divPVT.vel, pvt.vel / val)
                    self.assertEqual(divPVT.t, pvt.t)

                    inPlaceDivPVT = coordConv.PVT(pvt.pos, pvt.vel, pvt.t)
                    inPlaceDivPVT /= val
                    self.assertEqual(inPlaceDivPVT.pos, divPVT.pos)
                    self.assertEqual(inPlaceDivPVT.vel, divPVT.vel)
                    self.assertEqual(inPlaceDivPVT.t, divPVT.t)
            
    def testGetPos(self):
        """Test pvt.getPos
        """
        pvt = coordConv.PVT(1.1, 2.2, 12345.0)
        for dt in (1235.5, -123.3):
            newt = pvt.t + dt
            self.assertAlmostEqual(pvt.getPos(newt), refGetPos(pvt, newt))
    
    def testIsValid(self):
        """Test pvt.isfinite
        """
        pvt = coordConv.PVT(1, 2, 3)
        
        self.assertTrue(pvt.isfinite())
        pvt.pos = numpy.nan
        self.assertFalse(pvt.isfinite())
        pvt.pos = 1.0
        self.assertTrue(pvt.isfinite())

        pvt.vel = numpy.nan
        self.assertFalse(pvt.isfinite())
        pvt.vel = 1.0
        self.assertTrue(pvt.isfinite())
        
        pvt.t = numpy.nan
        self.assertFalse(pvt.isfinite())
        pvt.t = 0.0
        self.assertTrue(pvt.isfinite())
    
    def testInvalidate(self):
        """Test pvt.invalidate
        """
        pvt = coordConv.PVT(1, 2, 3)
        self.assertTrue(pvt.isfinite())
        pvt.invalidate()
        self.assertFalse(pvt.isfinite())
        self.assertFalse(numpy.isfinite(pvt.pos))
        self.assertFalse(numpy.isfinite(pvt.vel))
        self.assertFalse(numpy.isfinite(pvt.t))

        pvt2 = coordConv.PVT(-2, -4, 6)
        pvt2.invalidate(5)
        self.assertFalse(pvt2.isfinite())
        self.assertFalse(numpy.isfinite(pvt2.pos))
        self.assertFalse(numpy.isfinite(pvt2.vel))
        self.assertEqual(pvt2.t, 5.0)

    def testRot2D(self):
        """Test rot2D
        """
        def pvtIter():
            for pos in (5, -3):
                for vel in (0.1, 0, -0.3):
                    for tai in (500, 999):
                        yield coordConv.PVT(pos, vel, tai)

        for fromPVTX in pvtIter():
            for fromPVTY in pvtIter():
                for ang in (0, 21, -75.5):
                    for rotTAI in (fromPVTX.t - 200, fromPVTX.t + 5000):
                        toPVTX = coordConv.PVT()
                        toPVTY = coordConv.PVT()
                        coordConv.rot2D(toPVTX, toPVTY, fromPVTX, fromPVTY, ang, rotTAI)

                        for testTAI in (rotTAI, rotTAI + 1010):
                            fromX = fromPVTX.getPos(testTAI)
                            fromY = fromPVTY.getPos(testTAI)
                            predToX, predToY = coordConv.rot2D(fromX, fromY, ang)
                            self.assertAlmostEqual(predToX, toPVTX.getPos(testTAI))
                            self.assertAlmostEqual(predToY, toPVTY.getPos(testTAI))

    def testPolarFromXY(self):
        """Test polarFromXY and xyFromPolar
        """
        for xPos, yPos, predAtPole in (
            ( 1,  0, False),
            (-1,  0, False),
            ( 0,  1, False),
            ( 0, -1, False),
            ( 1,  1, False),
            ( 1, -1, False),
            (-1,  1, False),
            (-1, -1, False),
            (-123.45, -123.45, False),
            (0, 0, True),
        ):
            for xVel in (-1, 0, 1):
                for yVel in (-1, 0, 1):
                    x = coordConv.PVT(xPos, xVel, 10)
                    y = coordConv.PVT(yPos, yVel, 10)
                    r = coordConv.PVT()
                    theta = coordConv.PVT()
                    for endTime in (5, 10, 15):
                        atPole = coordConv.polarFromXY(r, theta, x, y, endTime)
                        refAtPole, refR, refTheta = refPolarFromXY(x, y, endTime)
                        self.assertEqual(atPole, refAtPole)
                        self.assertEqual(r.t, endTime)
                        self.assertEqual(refR.t, endTime)
                        self.assertEqual(theta.t, endTime)
                        self.assertEqual(refTheta.t, endTime)
                        self.assertTrue(numpy.allclose(
                            (r.pos, r.vel, theta.pos, theta.vel),
                            (refR.pos, refR.vel, refTheta.pos, refTheta.vel),
                        ))
                        refX = coordConv.PVT()
                        refY = coordConv.PVT()
                        coordConv.xyFromPolar(refX, refY, r, theta, endTime)
                        self.assertTrue(numpy.allclose(
                            (x.pos, x.vel, y.pos, y.vel),
                            (refX.getPos(10), refX.vel, refY.getPos(10), refY.vel),
                        ))


if __name__ == '__main__':
    unittest.main()
