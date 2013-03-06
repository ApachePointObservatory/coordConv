#!/usr/bin/env python
import unittest
import numpy
from coordConv import PVT

def predPos(pvt, t):
    return pvt.pos + (pvt.vel * (t - pvt.t))

class TestPVT(unittest.TestCase):
    def testConstructors(self):
        """Test PVT constructors
        """
        pvt = PVT(1.0, 2.0, 3.0)
        self.assertEquals(pvt.pos, 1.0)
        self.assertEquals(pvt.vel, 2.0)
        self.assertEquals(pvt.t, 3.0)
    
    def testAddPVT(self):
        """Test pvt + pvt
        """
        pvt1 = PVT(1.0, 2.0, 3.0)
        pvt2 = PVT(-1.5, -3.0, 4.0)
        pvt3 = pvt1 + pvt2
        predPos3 = pvt1.pos + predPos(pvt2, pvt1.t)
        self.assertAlmostEqual(pvt3.pos, predPos3)
        self.assertAlmostEqual(pvt3.vel, -1.0)
        self.assertAlmostEqual(pvt3.t, 3.0)

    
    def testSubtractPVT(self):
        """Test pvt - pvt
        """
        pvt1 = PVT(1.0, 2.0, 3.0)
        pvt2 = PVT(-1.5, -3.0, 4.0)
        pvt3 = pvt1 - pvt2
        predPos3 = pvt1.pos - predPos(pvt2, pvt1.t)
        self.assertAlmostEqual(pvt3.pos, predPos3)
        self.assertAlmostEqual(pvt3.vel, 5.0)
        self.assertAlmostEqual(pvt3.t, 3.0)

        pvt1 -= pvt2
        self.assertAlmostEqual(pvt1.pos, predPos3)
        self.assertAlmostEqual(pvt1.vel, 5.0)
        self.assertAlmostEqual(pvt1.t, 3.0)
    
    def testAddSubtractScalar(self):
        """Test pvt + scalar, pvt += scalar, pvt - scalar and pvt -= scalar
        """
        for pvt in (
            PVT(1.0, 2.0, 3.0),
            PVT(-1.5, -3.0, 4.0)
        ):
            for val in (-1.234e97, -25.3, -1.0e-13, 0.0, 1.0e13, 39.6, 1.34e99):
                addPVT = pvt + val
                self.assertAlmostEqual(addPVT.pos, pvt.pos + val)
                self.assertEqual(addPVT.vel, pvt.vel)
                self.assertEqual(addPVT.t, pvt.t)
                
                inPlaceAddPVT = PVT(pvt.pos, pvt.vel, pvt.t)
                inPlaceAddPVT += val
                self.assertEqual(inPlaceAddPVT.pos, addPVT.pos)
                self.assertEqual(inPlaceAddPVT.vel, addPVT.vel)
                self.assertEqual(inPlaceAddPVT.t, addPVT.t)

                subPVT = pvt - val
                self.assertAlmostEqual(subPVT.pos, pvt.pos - val)
                self.assertEqual(subPVT.vel, pvt.vel)
                self.assertEqual(subPVT.t, pvt.t)

                inPlaceSubPVT = PVT(pvt.pos, pvt.vel, pvt.t)
                inPlaceSubPVT -= val
                self.assertEqual(inPlaceSubPVT.pos, subPVT.pos)
                self.assertEqual(inPlaceSubPVT.vel, subPVT.vel)
                self.assertEqual(inPlaceSubPVT.t, subPVT.t)
    
    def testMultDivScalar(self):
        """Test pvt * scalar, pvt *= scalar, pvt / scalar and pvt /= scalar
        """
        for pvt in (
            PVT(1.0, 2.0, 3.0),
            PVT(-1.5, -3.0, 4.0)
        ):
            for val in (-1.234e97, -25.3, -1.0e-13, 0.0, 1.0e13, 39.6, 1.34e99):
                multPVT = pvt * val
                self.assertAlmostEqual(multPVT.pos, pvt.pos * val)
                self.assertEqual(multPVT.vel, pvt.vel * val)
                self.assertEqual(multPVT.t, pvt.t)

                inPlaceMultPVT = PVT(pvt.pos, pvt.vel, pvt.t)
                inPlaceMultPVT *= val
                self.assertEqual(inPlaceMultPVT.pos, multPVT.pos)
                self.assertEqual(inPlaceMultPVT.vel, multPVT.vel)
                self.assertEqual(inPlaceMultPVT.t, multPVT.t)

                if val != 0.0:
                    divPVT = pvt / val
                    self.assertAlmostEqual(divPVT.pos, pvt.pos / val)
                    self.assertEqual(divPVT.vel, pvt.vel / val)
                    self.assertEqual(divPVT.t, pvt.t)

                    inPlaceDivPVT = PVT(pvt.pos, pvt.vel, pvt.t)
                    inPlaceDivPVT /= val
                    self.assertEqual(inPlaceDivPVT.pos, divPVT.pos)
                    self.assertEqual(inPlaceDivPVT.vel, divPVT.vel)
                    self.assertEqual(inPlaceDivPVT.t, divPVT.t)
            
    def testGetPos(self):
        """Test pvt.getPos
        """
        pvt = PVT(1.1, 2.2, 12345.0)
        for dt in (1235.5, -123.3):
            newt = pvt.t + dt
            self.assertAlmostEqual(pvt.getPos(newt), predPos(pvt, newt))
    
    def testIsValid(self):
        """Test pvt.isValid
        """
        pvt = PVT(1, 2, 3)
        
        self.assertTrue(pvt.isValid())
        pvt.pos = numpy.nan
        self.assertFalse(pvt.isValid())
        pvt.pos = 1.0
        self.assertTrue(pvt.isValid())

        pvt.vel = numpy.nan
        self.assertFalse(pvt.isValid())
        pvt.vel = 1.0
        self.assertTrue(pvt.isValid())
        
        pvt.t = numpy.nan
        self.assertFalse(pvt.isValid())
        pvt.t = 0.0
        self.assertTrue(pvt.isValid())

if __name__ == '__main__':
    unittest.main()
