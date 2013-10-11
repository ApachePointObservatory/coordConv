#!/usr/bin/env python
import unittest
import numpy
import coordConv

class TestCoord(unittest.TestCase):
    def checkPVTCoord(self, pvtCoord):
        coord = pvtCoord.getCoord()
        vel = pvtCoord.getVel()
        orient = pvtCoord.getOrient()
        tai = pvtCoord.getTAI()
            
        atPole, equatAng, polarAng = coord.getSphPos()
        
        if abs(vel) > 1e-15:
            distList = (0, 0.001, 0.01, 0.1, 1, 10, 100, 179)
        else:
            distList = (0,)
        for dist in distList:
            if dist == 0:
                dt = 0
            else:
                dt =  dist / vel
            coord1 = pvtCoord.getCoord(tai + dt)
            measDist = coord.angularSeparation(coord1)
            distErr = (measDist - dist) / max(dist, 1e-7)
            self.assertLess(abs(distErr), 1e-7)
            if dist > 1e-7:
                measOrient = coord.orientationTo(coord1)
                if not numpy.isfinite(measOrient):
                    continue
                self.assertAlmostEqual(coord.orientationTo(coord1), orient)
        
        equatPVT = coordConv.PVT()
        polarPVT = coordConv.PVT()
        atPole = pvtCoord.getSphPVT(equatPVT, polarPVT, tai)
        self.assertEqual(atPole, coord.atPole() or atPole)
        self.assertAlmostEqual(coordConv.wrapCtr(equatPVT.pos - equatAng), 0)
        self.assertAlmostEqual(equatPVT.t, tai)
        self.assertAlmostEqual(polarPVT.pos, polarAng)
        self.assertAlmostEqual(polarPVT.t, tai)

    def testOneCoordConstructor(self):
        """Test PVTCoord(coord, orient, vel, tai) constructor
        """
        for equatAng in (0, 71, -123.4):
            for polarAng in (0, -75, -89.999999, 89.999999):
                coord = coordConv.Coord(equatAng, polarAng)
                for orient in (0, -45, 31.23):
                    for vel in (0, 0.1, 0.23):
                        for tai in (500.5, 10001.3):
                            pvtCoord = coordConv.PVTCoord(coord, orient, vel, tai)
                            self.assertAlmostEqual(pvtCoord.getOrient(), orient)
                            self.assertAlmostEqual(pvtCoord.getVel(), vel)
                            self.assertAlmostEqual(pvtCoord.getTAI(), tai)
                            self.assertAlmostEqual(pvtCoord.getCoord().angularSeparation(coord), 0)
                            self.assertTrue(pvtCoord.isfinite())
                            self.checkPVTCoord(pvtCoord)
    
    def testTwoCoordConstructor(self):
        """Test PVTCoord(coord0, coord1, tai, deltaT, defOrient)
        """
        numSkipped = 0
        numDefOrient = 0
        numNotDefOrient = 0
        for equatAng in (0, 71, -123.4):
            for polarAng in (0, -75, -89.999999, 89.999999):
                coord0 = coordConv.Coord(equatAng, polarAng)
                for orient in (0, -37.5, 90, 128):
                    for tai in (500.5, 10001.3):
                        for deltaT in (1e-7, 1e-5, 1e-3, 0.1, 10):
                            # test that default orientation is used when the coords are too close together
                            pvtCoord = coordConv.PVTCoord(coord0, coord0, tai, deltaT, orient)
                            self.assertAlmostEqual(pvtCoord.getOrient(), orient)
                            self.assertAlmostEqual(pvtCoord.getVel(), 0)
                            self.assertAlmostEqual(pvtCoord.getTAI(), tai)
                            self.assertTrue(pvtCoord.isfinite())
                            self.checkPVTCoord(pvtCoord)
                            
                            # test cases where coord1 may not equal coord0
                            for dist in (0, 0.001, 0.01, 0.1, 1, 10, 100, 179):
                                predVel = dist / float(deltaT)
                                if predVel > 100:
                                    numSkipped+= 1
                                    continue
                                coord1, toOrient = coord0.offset(orient, dist)
                                defOrient = orient + 10
                                pvtCoord = coordConv.PVTCoord(coord0, coord1, tai, deltaT, defOrient)
                                if numpy.isfinite(coord0.orientationTo(coord1)):
                                    numNotDefOrient += 1
                                    self.assertAlmostEqual(pvtCoord.getOrient(), orient)
                                    self.assertAlmostEqual(pvtCoord.getVel(), predVel)
                                else:
                                    numDefOrient += 1
                                    self.assertAlmostEqual(pvtCoord.getOrient(), defOrient)
                                    self.assertAlmostEqual(pvtCoord.getVel(), 0)
                                self.assertAlmostEqual(pvtCoord.getTAI(), tai)
                                self.assertTrue(pvtCoord.isfinite())
                                self.checkPVTCoord(pvtCoord)
                                
                                # make sure deltaT = 0 raises RuntimeError
                                self.assertRaises(RuntimeError, coordConv.PVTCoord, coord0, coord1, tai, 0, orient)
        self.assertGreater(numNotDefOrient, 100)
        self.assertGreater(numDefOrient, 100)

    def testEquality(self):
        """Test operator== and operator!=
        """
        def pvtCoordIter():
            for equatAng in (0, 71, -123.4):
                for polarAng in (0, -75, -89.999999, 89.999999):
                    coord = coordConv.Coord(equatAng, polarAng)
                    for orient in (0, -45, 31.23):
                        for vel in (0, 0.1, 0.23):
                            for tai in (500.5, 10001.3):
                                yield coordConv.PVTCoord(coord, orient, vel, tai)

        for pvtCoord1 in pvtCoordIter():
            for pvtCoord2 in pvtCoordIter():
                if pvtCoord1.getCoord() == pvtCoord2.getCoord() \
                    and pvtCoord1.getOrient() == pvtCoord2.getOrient() \
                    and pvtCoord1.getVel() == pvtCoord2.getVel() \
                    and pvtCoord1.getTAI() == pvtCoord2.getTAI():
                    self.assertTrue(pvtCoord1 == pvtCoord2)
                else:
                    self.assertTrue(pvtCoord1 != pvtCoord2)
                self.assertNotEqual(pvtCoord1 == pvtCoord2, pvtCoord1 != pvtCoord2)

    def testOffset(self):
        """Test PVTCoord.offset
        """
        numAtPole = 0
        numNotAtPole = 0
        toOrientPVT = coordConv.PVT()
        for equatAng in (0, 71, -123.4):
            for polarAng in (0, -75, -89.99999, -90, 89.99999):
                coord = coordConv.Coord(equatAng, polarAng)
                for orient in (0, 31.23):
                    for vel in (0, 0.23):
                        for tai in (500.5, 10001.3):
                            pvtCoord = coordConv.PVTCoord(coord, orient, vel, tai)
                            for offOrient in (0, -72):
                                for offOrientVel in (0, 0.1):
                                    for offOrientTAI in (tai, orient + 10):
                                        offOrientPVT = coordConv.PVT(offOrient, offOrientVel, offOrientTAI)
                                        for offDist in (0, 0.1, 1):
                                            for offDistVel in (0, 0.3):
                                                for offDistTAI in (tai, orient - 5):
                                                    offDistPVT = coordConv.PVT(offDist, offDistVel, offDistTAI)
                                                    if coord.atPole():
                                                        numAtPole += 1
                                                        self.assertRaises(RuntimeError, pvtCoord.offset, toOrientPVT, offOrientPVT, offDistPVT, offTAI)
                                                    else:
                                                        numNotAtPole += 1
                                                        for offTAI in (tai, tai - 1, tai + 2):
                                                            offPVTCoord = pvtCoord.offset(toOrientPVT, offOrientPVT, offDistPVT, offTAI)
                                                            self.assertTrue(offPVTCoord.isfinite())
                                                            offOrientAtOffTAI = offOrientPVT.getPos(offTAI)
                                                            offDistAtOffTAI = offDistPVT.getPos(offTAI)
                                                            toOrientAtOffTAI = toOrientPVT.getPos(offTAI)
                                                            offCoordAtOffTAI = offPVTCoord.getCoord(offTAI)
                                                            coordAtOffTAI = pvtCoord.getCoord(offTAI)
                                                            predOffCoord, predToOrient = coordAtOffTAI.offset(offOrientAtOffTAI, offDistAtOffTAI)
                                                            
                                                            self.assertAlmostEqual(toOrientAtOffTAI, predToOrient)
                                                            self.assertAlmostEqual(offCoordAtOffTAI.angularSeparation(predOffCoord), 0)
        self.assertGreater(numNotAtPole, 100)
        self.assertGreater(numAtPole, 100)


if __name__ == '__main__':
    unittest.main()
