#!/usr/bin/env python
from __future__ import absolute_import, division

import itertools
import unittest
import math

import numpy

import coordConv

class TestCoord(unittest.TestCase):
    def testOneCoordConstructor(self):
        """Test PVTCoord(coord, vel, tai) constructor
        """
        for equatAng in (0, 71, -123.4):
            for polarAng in (0, -75, -89.999999, 89.999999):
                coord = coordConv.Coord(equatAng, polarAng)
                for vel in (
                    (0, 0, 0),
                    (10, -23, 1.23),
                ):
                    vel = numpy.array(vel, dtype=numpy.float64)
                    for tai in (500.5, 10001.3):
                        pvtCoord = coordConv.PVTCoord(coord, vel, tai)
                        self.assertEqual(coord, pvtCoord.getCoord())
                        self.assertEqual(coord, pvtCoord.getCoord(tai))
                        self.assertTrue(numpy.allclose(vel, pvtCoord.getVel()))
                        self.assertAlmostEqual(pvtCoord.getTAI(), tai)
                        self.assertTrue(pvtCoord.isfinite())

                        for parallax in (0, 0.012):
                            for equatPM in (0, 0.11):
                                for polarPM in (0, -0.12):
                                    for radVel in (0, 0.13):
                                        coordWithPM = coordConv.Coord(equatAng, polarAng, parallax, equatPM, polarPM, radVel)
                                        pvtCoordWithPM = coordConv.PVTCoord(coordWithPM, vel, tai)
                                        self.assertEqual(coordWithPM, pvtCoordWithPM.getCoord())
                                        self.assertEqual(coordWithPM, pvtCoordWithPM.getCoord(tai))
                                        self.assertTrue(numpy.allclose(vel, pvtCoordWithPM.getVel()))
                                        self.assertAlmostEqual(pvtCoordWithPM.getTAI(), tai)
                                        self.assertTrue(pvtCoordWithPM.isfinite())

                                        self.assertAlmostEqual(pvtCoordWithPM.getCoord().angularSeparation(coord), 0)
                                        self.assertTrue(pvtCoordWithPM.isfinite())
                                        coordToCheck = pvtCoordWithPM.getCoord()
                                        self.assertAlmostEqual(radVel, coordToCheck.getRadVel())
                                        atPole, checkEquatPM, checkPolarPM = coordToCheck.getPM()
                                        if not atPole:
                                            self.assertAlmostEqual(equatPM, checkEquatPM)
                                            self.assertAlmostEqual(polarPM, checkPolarPM)
                                        self.assertAlmostEqual(parallax, coordToCheck.getParallax())

    def testTwoCoordConstructor(self):
        """Test PVTCoord(coord0, coord1, tai, deltaT, defOrient)
        """
        numRan = 0
        for equatAng in (0, 71, -123.4):
            for polarAng in (0, -75, -89.999999, 89.999999):
                coord0 = coordConv.Coord(equatAng, polarAng)
                for orient in (0, -37.5, 90, 128):
                    for tai in (500.5, 10001.3):
                        for deltaT in (1e-7, 1e-5, 1e-3, 0.1, 10):
                            # test that velocity is 0 when the coords are identical
                            pvtCoordA = coordConv.PVTCoord(coord0, coord0, tai, deltaT)
                            self.assertTrue(pvtCoordA.isfinite())
                            self.assertAlmostEqual(pvtCoordA.getTAI(), tai)

                            self.assertTrue(numpy.allclose(pvtCoordA.getVel(), (0,0,0)))
                            self.assertEqual(pvtCoordA.getCoord(), coord0)
                            self.assertEqual(pvtCoordA.getCoord(tai), coord0)

                            # test cases where coord1 may not equal coord0
                            for dist in (0, 0.001, 0.01, 0.1, 1, 10, 100, 179):
                                predVel = dist / float(deltaT)
                                if predVel > 100:
                                    continue
                                numRan += 1
                                coord1, toOrient = coord0.offset(orient, dist)
                                pvtCoordB = coordConv.PVTCoord(coord0, coord1, tai, deltaT)
                                self.assertEqual(pvtCoordB.getCoord(), coord0)
                                self.assertEqual(pvtCoordB.getCoord(tai), coord0)

                                gotCoord1 = pvtCoordB.getCoord(tai + deltaT)
                                self.assertTrue(numpy.allclose(coord1.getVecPos(), gotCoord1.getVecPos()))
                                self.assertTrue(numpy.allclose(coord1.getVecPM(), gotCoord1.getVecPM()))
        self.assertTrue(numRan > 1000)

    def testTwoPVTConstructors(self):
        """Test both two-PVT constructors:

        PVTCoord(equatPVT, polarPVT, tai, distPVT)
        PVTCoord(equatPVT, polarPVT, tai, distPVT, equatPM, polarPM, radVel, defOrient)
        """
        for equatAng in (0, 71, -123.4):
            for polarAng in (0, -75, -89.99, 89.99):
                for orient in (0, -45, 31.23):
                    sinOrient = coordConv.sind(orient)
                    cosOrient = coordConv.cosd(orient)
                    for vel in (0, 0.1, 0.23):
                        for tai in (500.5, 10001.3):
                            # coord0 = coordConv.Coord(equatAng, polarAng)
                            polarVel = vel * sinOrient
                            equatVel = vel * cosOrient / coordConv.cosd(polarAng)
                            equatPVT = coordConv.PVT(equatAng, equatVel, tai)
                            polarPVT = coordConv.PVT(polarAng, polarVel, tai)

                            pvtCoord = coordConv.PVTCoord(equatPVT, polarPVT)
                            self.assertAlmostEqual(pvtCoord.getTAI(), tai)
                            self.assertTrue(pvtCoord.isfinite())

                            gotEquatPVT = coordConv.PVT()
                            gotPolarPVT = coordConv.PVT()
                            pvtCoord.getSphPVT(gotEquatPVT, gotPolarPVT)
                            coordConv.assertPVTsAlmostEqual(equatPVT, gotEquatPVT, doWrap=True)
                            coordConv.assertPVTsAlmostEqual(polarPVT, gotPolarPVT, doWrap=False)

                            for parallax in (0, 0.012):
                                dist = coordConv.distanceFromParallax(parallax)
                                for distVel in (0, 10000):
                                    distPVT = coordConv.PVT(dist, distVel, tai)
                                    for equatPM in (0, 0.11):
                                        for polarPM in (0, -0.12):
                                            for radVel in (0, 0.13):
                                                pvtCoordPM = coordConv.PVTCoord(equatPVT, polarPVT, distPVT, equatPM, polarPM, radVel)
                                                self.assertAlmostEqual(pvtCoordPM.getTAI(), tai)
                                                self.assertTrue(pvtCoordPM.isfinite())

                                                gotEquatPVT = coordConv.PVT()
                                                gotPolarPVT = coordConv.PVT()
                                                pvtCoordPM.getSphPVT(gotEquatPVT, gotPolarPVT)
                                                coordConv.assertPVTsAlmostEqual(equatPVT, gotEquatPVT, doWrap=True)
                                                coordConv.assertPVTsAlmostEqual(polarPVT, gotPolarPVT, doWrap=False)

                                                coordToCheck = pvtCoordPM.getCoord()
                                                self.assertAlmostEqual(radVel, coordToCheck.getRadVel())
                                                atPole, checkEquatPM, checkPolarPM = coordToCheck.getPM()
                                                if not atPole:
                                                    self.assertAlmostEqual(equatPM, checkEquatPM)
                                                    self.assertAlmostEqual(polarPM, checkPolarPM)
                                                self.assertAlmostEqual(parallax, coordToCheck.getParallax())
                                                if parallax != 0:
                                                    self.assertAlmostEqual(dist, coordToCheck.getDistance(), places=5)
                                                    coordConv.assertPVTsAlmostEqual(distPVT, pvtCoordPM.getDistance(),
                                                        velPlaces = 5, posPlaces=5)


    def testEquality(self):
        """Test operator== and operator!=
        """
        def pvtCoordIter():
            for equatAng in (0, 71, -123.4):
                for polarAng in (0, -75, -89.999999, 89.999999):
                    for equatVel in (0, 0.023):
                        for polarVel in (0, -math.copysign(0.012, polarAng)):
                            for tai in (4889100000.5, 1000.1):
                                equatPVT = coordConv.PVT(equatAng, equatVel, tai)
                                polarPVT = coordConv.PVT(polarAng, polarVel, tai)
                                yield coordConv.PVTCoord(equatPVT, polarPVT)

        for pvtCoord1 in pvtCoordIter():
            for pvtCoord2 in pvtCoordIter():
                if pvtCoord1.getCoord() == pvtCoord2.getCoord() \
                    and tuple(pvtCoord1.getVel()) == tuple(pvtCoord2.getVel()) \
                    and pvtCoord1.getTAI() == pvtCoord2.getTAI():
                    self.assertTrue(pvtCoord1 == pvtCoord2)
                else:
                    self.assertTrue(pvtCoord1 != pvtCoord2)
                self.assertNotEqual(pvtCoord1 == pvtCoord2, pvtCoord1 != pvtCoord2)

    def testOffset(self):
        """Test PVTCoord.offset
        """
        def pvtCoordIter():
            for equatAng in (0, 71, -123.4):
                for polarAng in (0, -75, -89.99, 89.99, -90, 90):
                    for equatVel in (0, 0.023):
                        for polarVel in (0, -math.copysign(0.012, polarAng)):
                            for tai in (4889100000.5, 1000.1):
                                equatPVT = coordConv.PVT(equatAng, equatVel, tai)
                                polarPVT = coordConv.PVT(polarAng, polarVel, tai)
                                yield coordConv.PVTCoord(equatPVT, polarPVT)

        def offPVTIter(pvtCoord):
            """return a sequence of (offOrientPVT, offDistPVT)"""
            tai = pvtCoord.getTAI()
            for offOrient in (0, -72):
                for offOrientVel in (0, 0.1):
                    for offOrientTAI in (tai, tai + 10):
                        offOrientPVT = coordConv.PVT(offOrient, offOrientVel, offOrientTAI)
                        for offDist in (0, 0.1, 1):
                            for offDistVel in (0, 0.3):
                                for offDistTAI in (tai, tai - 5):
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
                    self.assertRaises(RuntimeError, pvtCoord.offset, toOrientPVT, offOrientPVT, offDistPVT)
                else:
                    numNotAtPole += 1
                    offPVTCoord = pvtCoord.offset(toOrientPVT, offOrientPVT, offDistPVT)
                    self.assertTrue(offPVTCoord.isfinite())
                    offOrientAtTAI = offOrientPVT.getPos(tai)
                    offDistAtTAI = offDistPVT.getPos(tai)
                    toOrientAtTAI = toOrientPVT.getPos(tai)
                    offCoordAtTAI = offPVTCoord.getCoord(tai)
                    coordAtTAI = pvtCoord.getCoord(tai)
                    predOffCoord, predToOrient = coordAtTAI.offset(offOrientAtTAI, offDistAtTAI)

                    self.assertAlmostEqual(toOrientAtTAI, predToOrient)
                    self.assertAlmostEqual(offCoordAtTAI.angularSeparation(predOffCoord), 0)

        self.assertGreater(numNotAtPole, 100)
        self.assertGreater(numAtPole, 100)

    def testAngularSeparation(self):
        """Test PVTCoord.angularSeparation and orientationTo
        """
        def pvtCoordIter():
            for equatAng in (0, 71, -123.4):
                for polarAng in (0, -75, -89.99999, -90, 89.99999):
                    for equatVel in (0, 0.023):
                        for polarVel in (0, -math.copysign(0.012, polarAng)):
                            for tai in (4889100000.5, 1000.1):
                                equatPVT = coordConv.PVT(equatAng, equatVel, tai)
                                polarPVT = coordConv.PVT(polarAng, polarVel, tai)
                                yield coordConv.PVTCoord(equatPVT, polarPVT)

        def refAngularSeparation(pvtCoord0, pvtCoord1):
            """Compute angular separation between pvtCoord0, pvtCoord1 using Coord.angularSeparation
            """
            tai = pvtCoord0.getTAI()
            DeltaT = 0.01
            posList = []
            for tempTAI in (tai, tai + DeltaT):
                coord0 = pvtCoord0.getCoord(tempTAI)
                coord1 = pvtCoord1.getCoord(tempTAI)
                posList.append(coord0.angularSeparation(coord1))
            return makePVTFromPair(posList, tai, DeltaT, True)

        def refOrientTo(pvtCoord0, pvtCoord1):
            """Compute orientation from pvtCoord0 to pvtCoord1 using Coord.orientationTo
            """
            tai = pvtCoord0.getTAI()
            DeltaT = 0.01
            posList = []
            for tempTAI in (tai, tai + DeltaT):
                coord0 = pvtCoord0.getCoord(tempTAI)
                coord1 = pvtCoord1.getCoord(tempTAI)
                posList.append(coord0.orientationTo(coord1))
            if numpy.all(numpy.isfinite(posList)):
                return makePVTFromPair(posList, tai, DeltaT, True)
            elif numpy.isfinite(posList[0]):
                return coordConv.PVT(posList[0], 0, tai)
            elif numpy.isfinite(posList[1]):
                return coordConv.PVT(posList[1], 0, tai)
            else:
                return coordConv.PVT()

        for pvtCoord0 in pvtCoordIter():
            for pvtCoord1 in pvtCoordIter():
                angSep = pvtCoord0.angularSeparation(pvtCoord1)
                refAngSep = refAngularSeparation(pvtCoord0, pvtCoord1)
                coordConv.assertPVTsAlmostEqual(angSep, refAngSep)

                if pvtCoord0 == pvtCoord1:
                    self.assertAlmostEqual(angSep.pos, 0)
                    self.assertAlmostEqual(angSep.vel, 0)

                orient = pvtCoord0.orientationTo(pvtCoord1)
                refOrient = refOrientTo(pvtCoord0, pvtCoord1)
                if not orient.isfinite():
                    self.assertFalse(refOrient.isfinite())
                else:
                    coordConv.assertPVTsAlmostEqual(orient, refOrient, doWrap=True)

    def testConvertFromVel(self):
        """Test velocity of convertFrom
        """
        taiDate = 4889900000.205
        site = coordConv.Site(-105.822616, 32.780988, 2788)
        icrsCoordSys = coordConv.ICRSCoordSys()
        appTopoCoordSys = coordConv.AppTopoCoordSys()

        # find ICRS coordinate of a sidereal point on the equator along the meridion
        appTopoCoord = coordConv.Coord(0, 90 - site.meanLat)
        icrsCoord = icrsCoordSys.convertFrom(appTopoCoordSys, appTopoCoord, site, taiDate)

        icrsPVTCoord = coordConv.PVTCoord(icrsCoord, icrsCoord, taiDate, 0.001)

        appTopoPVTCoord = appTopoCoordSys.convertFrom(icrsCoordSys, icrsPVTCoord, site)

        equatPVT = coordConv.PVT()
        polarPVT = coordConv.PVT()
        appTopoPVTCoord.getSphPVT(equatPVT, polarPVT)
        self.assertEqual(equatPVT.t, taiDate)
        self.assertEqual(polarPVT.t, taiDate)
        equatSpaceVel = equatPVT.vel * coordConv.cosd(polarPVT.pos)
        self.assertAlmostEqual(polarPVT.vel, 0, places=3)
        self.assertAlmostEqual(equatSpaceVel, -1/240.0, places=3) # 360 deg/day

        # check round trip of scale and orientation
        for fromDir in (0, 45):
            for fromVel in (0, 0.01):
                fromDirPVT = coordConv.PVT(fromDir, fromVel, taiDate)
                toDirPVT = coordConv.PVT()
                fromDir2PVT = coordConv.PVT()
                at2PVTCoord, scaleChange = appTopoCoordSys.convertFrom(toDirPVT, icrsCoordSys, icrsPVTCoord, fromDirPVT, site)
                icrs2PVTCoord, scaleChange2 = icrsCoordSys.convertFrom(fromDir2PVT, appTopoCoordSys, at2PVTCoord, toDirPVT, site)
                self.assertAlmostEqual(scaleChange, 1.0/scaleChange2, places=7)
                coordConv.assertPVTsAlmostEqual(fromDirPVT, fromDir2PVT, doWrap=True, velPlaces=6)

    def testConvertFrom(self):
        """Test a few instances of CoordSys.convertFrom on PVTCoords

        This test assumes that CoordSys.convertFrom works on Coords (tested elsewhere)
        """
        site = coordConv.Site(-105.822616, 32.780988, 2788)
        fromCoordSys = coordConv.ICRSCoordSys()
        toCoordSys = coordConv.GalCoordSys()
        dt = 0.001

        def pvtCoordIter():
            for equatAng in (0, 71, -123.4):
                for polarAng in (0, -75, -89.9, 89.9):
                    for equatVel in (0, 0.023):
                        for polarVel in (0, -math.copysign(0.012, polarAng)):
                            for tai in (4889100000.5, 1000.1):
                                equatPVT = coordConv.PVT(equatAng, equatVel, tai)
                                polarPVT = coordConv.PVT(polarAng, polarVel, tai)
                                yield coordConv.PVTCoord(equatPVT, polarPVT)

        for fromPVTCoord in pvtCoordIter():
            tai0 = fromPVTCoord.getTAI()
            taiPair = [tai0, tai0 + dt]
            fromCoordPair = [fromPVTCoord.getCoord(t) for t in taiPair]
            toPVTCoord = toCoordSys.convertFrom(fromCoordSys, fromPVTCoord, site)
            for fromCoord, tai in itertools.izip(fromCoordPair, taiPair):
                predToCoord = toCoordSys.convertFrom(fromCoordSys, fromCoord, site, tai)
                measToCoord = toPVTCoord.getCoord(tai)
                self.assertAlmostEqual(predToCoord.angularSeparation(measToCoord), 0)

    def testLunarVel(self):
        """Sanity-check lunar tracking velocity

        This checks for an issue we had with tracking close objects: position was right,
        but velocity was not. This was because Apparent Geocentric date was not being updated
        inside the PVTCoord version of CoordSys.convertFrom.
        """
        tai = 4914602887
        dt = 0.1
        site = coordConv.Site(-105.822616, 32.780988, 2788)
        geoCoordSys = coordConv.AppGeoCoordSys()
        topoCoordSys = coordConv.AppTopoCoordSys()
        geoCoord = coordConv.Coord(0, 75, 82505922)
        topoPVTList = []
        for evalTAI in (tai, tai + dt):
            geoPVTCoord = coordConv.PVTCoord(geoCoord, geoCoord, evalTAI, 0.01)
            topoPVTCoord = topoCoordSys.convertFrom(geoCoordSys, geoPVTCoord, site)
            topoPVTPair = [coordConv.PVT() for i in range(2)]
            topoPVTCoord.getSphPVT(topoPVTPair[0], topoPVTPair[1])
            topoPVTList.append(topoPVTPair)
        for i in range(2):
            pvt0 = topoPVTList[0][i]
            pvt1 = topoPVTList[1][i]
            coordConv.assertPVTsAlmostEqual(pvt0.copy(pvt1.t), pvt1)

def makePVTFromPair(posPair, tai, deltaT, isAngle):
    pos = posPair[0]
    if (isAngle):
        vel = coordConv.wrapCtr(posPair[1] - posPair[0]) / deltaT
    else:
        vel = (posPair[1] - posPair[0]) / deltaT
    return coordConv.PVT(pos, vel, tai)


if __name__ == '__main__':
    unittest.main()
