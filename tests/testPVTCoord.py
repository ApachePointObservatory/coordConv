#!/usr/bin/env python
import itertools
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

                            for parallax in (0, 0.012):
                                for equatPM in (0, 0.11):
                                    for polarPM in (0, -0.12):
                                        for radVel in (0, 0.13):
                                            coordWithPM = coordConv.Coord(equatAng, polarAng, parallax, equatPM, polarPM, radVel)
                                            pvtCoord = coordConv.PVTCoord(coordWithPM, orient, vel, tai)
                                            self.assertAlmostEqual(pvtCoord.getOrient(), orient)
                                            self.assertAlmostEqual(pvtCoord.getVel(), vel)
                                            self.assertAlmostEqual(pvtCoord.getTAI(), tai)
                                            self.assertAlmostEqual(pvtCoord.getCoord().angularSeparation(coord), 0)
                                            self.assertTrue(pvtCoord.isfinite())
                                            self.checkPVTCoord(pvtCoord)
                                            coordToCheck = pvtCoord.getCoord()
                                            self.assertAlmostEqual(radVel, coordToCheck.getRadVel())
                                            atPole, checkEquatPM, checkPolarPM = coordToCheck.getPM()
                                            if not atPole:
                                                self.assertAlmostEqual(equatPM, checkEquatPM)
                                                self.assertAlmostEqual(polarPM, checkPolarPM)
                                            self.assertAlmostEqual(parallax, coordToCheck.getParallax())
    
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

    def testTwoPVTConstructors(self):
        """Test both two-PVT constructors:

        PVTCoord(equatPVT, polarPVT, tai, parallax, defOrient)
        PVTCoord(equatPVT, polarPVT, tai, parallax, equatPM, polarPM, radVel, defOrient)
        """
        for equatAng in (0, 71, -123.4):
            for polarAng in (0, -75, -89.99, 89.99):
                for orient in (0, -45, 31.23):
                    sinOrient = coordConv.sind(orient)
                    cosOrient = coordConv.cosd(orient)
                    for vel in (0, 0.1, 0.23):
                        for tai in (500.5, 10001.3):
                            coord0 = coordConv.Coord(equatAng, polarAng)
                            polarVel = vel * sinOrient
                            equatVel = vel * cosOrient / coordConv.cosd(polarAng)
                            equatPVT = coordConv.PVT(equatAng, equatVel, tai)
                            polarPVT = coordConv.PVT(polarAng, polarVel, tai)

                            desOrient = orient if vel > 1e-15 else 0
                            pvtCoord = coordConv.PVTCoord(equatPVT, polarPVT, tai)
                            self.assertAlmostEqual(pvtCoord.getOrient(), desOrient)
                            self.assertAlmostEqual(pvtCoord.getVel(), vel)
                            self.assertAlmostEqual(pvtCoord.getTAI(), tai)
                            self.assertTrue(pvtCoord.isfinite())
                            self.checkPVTCoord(pvtCoord)

                            for parallax in (0, 0.012):
                                for equatPM in (0, 0.11):
                                    for polarPM in (0, -0.12):
                                        for radVel in (0, 0.13):
                                            pvtCoord = coordConv.PVTCoord(equatPVT, polarPVT, tai, parallax, equatPM, polarPM, radVel)
                                            self.assertAlmostEqual(pvtCoord.getOrient(), desOrient)
                                            self.assertAlmostEqual(pvtCoord.getVel(), vel)
                                            self.assertAlmostEqual(pvtCoord.getTAI(), tai)
                                            self.assertTrue(pvtCoord.isfinite())
                                            self.checkPVTCoord(pvtCoord)
                                            coordToCheck = pvtCoord.getCoord()
                                            self.assertAlmostEqual(radVel, coordToCheck.getRadVel())
                                            atPole, checkEquatPM, checkPolarPM = coordToCheck.getPM()
                                            if not atPole:
                                                self.assertAlmostEqual(equatPM, checkEquatPM)
                                                self.assertAlmostEqual(polarPM, checkPolarPM)
                                            self.assertAlmostEqual(parallax, coordToCheck.getParallax())


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
        def pvtCoordIter():
            """Return a sequence of PVTCoord"""
            for equatAng in (0, 71, -123.4):
                for polarAng in (0, -75, -89.99999, -90, 89.99999):
                    coord = coordConv.Coord(equatAng, polarAng)
                    for orient in (0, 31.23):
                        for vel in (0, 0.23):
                            if coord.atPole() and vel != 0:
                                continue # do not try to make invalid PVTCoords
                            for tai in (500.5, 10001.3):
                                yield coordConv.PVTCoord(coord, orient, vel, tai)

        def offPVTIter(pvtCoord):
            """return a sequence of (offOrientPVT, offDistPVT)"""
            tai = pvtCoord.getTAI()
            orient = pvtCoord.getOrient()
            for offOrient in (0, -72):
                for offOrientVel in (0, 0.1):
                    for offOrientTAI in (tai, orient + 10):
                        offOrientPVT = coordConv.PVT(offOrient, offOrientVel, offOrientTAI)
                        for offDist in (0, 0.1, 1):
                            for offDistVel in (0, 0.3):
                                for offDistTAI in (tai, orient - 5):
                                    offDistPVT = coordConv.PVT(offDist, offDistVel, offDistTAI)
                                    yield offOrientPVT, offDistPVT

        numAtPole = 0
        numNotAtPole = 0
        toOrientPVT = coordConv.PVT()
        for pvtCoord in pvtCoordIter():
            tai = pvtCoord.getTAI()
            for offOrientPVT, offDistPVT in offPVTIter(pvtCoord):
                if pvtCoord.getCoord().atPole():
                    numAtPole += 1
                    self.assertRaises(RuntimeError, pvtCoord.offset, toOrientPVT, offOrientPVT, offDistPVT, tai)
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

    def testAngularSeparation(self):
        """Test PVTCoord.angularSeparation and orientationTo
        """
        def pvtCoordIter():
            for equatAng in (0, 71, -123.4):
                for polarAng in (0, -75, -89.99999, -90, 89.99999):
                    coord = coordConv.Coord(equatAng, polarAng)
                    for orient in (0, 31.23):
                        for vel in (0, 0.023):
                            if coord.atPole() and vel != 0:
                                continue # do not try to make invalid PVTCoords
                            for tai in (10002.5, 10001.3):
                                yield coordConv.PVTCoord(coord, orient, vel, tai)

        def refAngularSeparation(pvtCoord0, pvtCoord1, tai):
            """Compute angular separation between pvtCoord0, pvtCoord1 using Coord.angularSeparation
            """
            DeltaT = 0.01
            posList = []
            for tempTAI in (tai, tai + DeltaT):
                coord0 = pvtCoord0.getCoord(tempTAI)
                coord1 = pvtCoord1.getCoord(tempTAI)
                posList.append(coord0.angularSeparation(coord1))
            return makePVTFromPair(posList, tai, DeltaT, True)

        def refOrientTo(pvtCoord0, pvtCoord1, tai):
            """Compute orientation from pvtCoord0 to pvtCoord1 using Coord.orientationTo
            """
            DeltaT = 0.01
            posList = []
            for tempTAI in (tai, tai + DeltaT):
                coord0 = pvtCoord0.getCoord(tempTAI)
                coord1 = pvtCoord1.getCoord(tempTAI)
                posList.append(coord0.orientationTo(coord1))
            return makePVTFromPair(posList, tai, DeltaT, True)

        for pvtCoord0 in pvtCoordIter():
            tai0 = pvtCoord0.getTAI()
            for pvtCoord1 in pvtCoordIter():
                for tai in (tai0 + 5, tai0 - 79.2):
                    angSep01 = pvtCoord0.angularSeparation(pvtCoord1, tai)
                     # compute reference separation, which should be the same for 0->1 and 1->0
                    refAngSep = refAngularSeparation(pvtCoord0, pvtCoord1, tai)
                    self.assertPVTsAlmostEqual(angSep01, refAngSep, isAngle=False)

                    angSep10 = pvtCoord1.angularSeparation(pvtCoord0, tai)
                    self.assertPVTsAlmostEqual(angSep10, refAngSep)

                    if pvtCoord0 == pvtCoord1:
                        self.assertAlmostEqual(angSep01.pos, 0)

                    orientTo01 = pvtCoord0.orientationTo(pvtCoord1, tai)
                    refOrientTo01 = refOrientTo(pvtCoord0, pvtCoord1, tai)
                    if not orientTo01.isfinite():
                        self.assertFalse(refOrientTo01.isfinite())
                    else:
                        self.assertPVTsAlmostEqual(orientTo01, refOrientTo01, isAngle=True)


                    orientTo10 = pvtCoord1.orientationTo(pvtCoord0, tai)
                    refOrientTo10 = refOrientTo(pvtCoord1, pvtCoord0, tai)
                    if not orientTo10.isfinite():
                        self.assertFalse(refOrientTo10.isfinite())
                    else:
                        self.assertPVTsAlmostEqual(orientTo10, refOrientTo10, isAngle=True)

    def testConvertFrom(self):
        """Test a few instances of CoordSys.convertFrom on PVTCoords
        
        This test assumes that CoordSys.convertFrom works on Coords (tested elsewhere)
        """
        site = coordConv.Site(-105.822616, 32.780988, 2788)
        fromCoordSys = coordConv.ICRSCoordSys()
        toCoordSys = coordConv.GalCoordSys()
        cnvTAI = 4889900000.2
        dt = 0.001
        taiPair = (cnvTAI, cnvTAI + dt)

        def pvtCoordIter():
            for equatAng in (0, 71, -123.4):
                for polarAng in (0, -75, -89.9, 89.9):
                    coord = coordConv.Coord(equatAng, polarAng)
                    for orient in (0, 31.23):
                        for vel in (0, 0.023):
                            if coord.atPole() and vel != 0:
                                continue # do not try to make invalid PVTCoords
                            for tai in (4889100000.5, 1000.1):
                                yield coordConv.PVTCoord(coord, orient, vel, tai)

        for fromPVTCoord in pvtCoordIter():
            fromCoordPair = [fromPVTCoord.getCoord(t) for t in taiPair]
            toPVTCoord = toCoordSys.convertFrom(fromCoordSys, fromPVTCoord, site, cnvTAI)
            for fromCoord, tai in itertools.izip(fromCoordPair, taiPair):
                predToCoord = toCoordSys.convertFrom(fromCoordSys, fromCoord, site)
                measToCoord = toPVTCoord.getCoord(tai)
                self.assertAlmostEqual(predToCoord.angularSeparation(measToCoord), 0)

            toDir = coordConv.PVT()
            for fromDirPos in (45, -150):
                for fromDirVel in (0, 0.1):
                    fromDir = coordConv.PVT(fromDirPos, fromDirVel, cnvTAI)
                    toPVTCoord, scaleChange = toCoordSys.convertFrom(toDir, fromCoordSys, fromPVTCoord, fromDir, site, cnvTAI)
                    for fromCoord, tai in itertools.izip(fromCoordPair, taiPair):
                        predToCoord, predToDir, predScaleChange = toCoordSys.convertFrom(fromCoordSys, fromCoord, fromDir.getPos(tai), site)
                        measToCoord = toPVTCoord.getCoord(tai)
                        self.assertAlmostEqual(predToCoord.angularSeparation(measToCoord), 0)
                        self.assertAlmostEqual(toDir.getPos(tai), predToDir)
                        self.assertAlmostEqual(scaleChange, predScaleChange)
                        self.assertAlmostEqual(scaleChange, 1.0) # no scale change because mean to mean conversion


    def assertPVTsAlmostEqual(self, pvt0, pvt1, posDig=7, velDig=7, isAngle=False):
        """Compare two PVTS; both must have the same time
        """
        self.assertEqual(pvt0.t, pvt1.t)
        if isAngle:
            self.assertAlmostEqual(coordConv.wrapCtr(pvt0.pos - pvt1.pos), 0, posDig)
        else:
            self.assertAlmostEqual(pvt0.pos, pvt1.pos, posDig)
        self.assertAlmostEqual(pvt0.vel, pvt1.vel, velDig)

def makePVTFromPair(posPair, tai, deltaT, isAngle):
    pos = posPair[0];
    if (isAngle):
        vel = coordConv.wrapCtr(posPair[1] - posPair[0]) / deltaT
    else:
        vel = (posPair[1] - posPair[0]) / deltaT
    return coordConv.PVT(pos, vel, tai)


if __name__ == '__main__':
    unittest.main()
