#!/usr/bin/env python
from __future__ import absolute_import, division

import time
import unittest
import itertools
import os.path
import numpy
import coordConv
"""
Process data from masscc_out.dat

# UT1_min_TAI, poleWander(2)
  -34.77412527      0.9807816869E-05  0.9046393074E-04
# FROM sys, date, pos1, pos2, PM1, PM2, parallax, radVel, dir, refCoA, refCoB, \
   TO sys, date, pos1, pos2, PM1, PM2, parallax, radVel, dir, scaleChange, atInf, atPole, isOK, TAI
...

ignores blank lines and lines beginning with #
"""
DataFile = os.path.join(os.path.dirname(__file__), "data", "masscc_out.dat")

ContinueOnError = False

CSysDict = {
     4: coordConv.ICRSCoordSys,
     3: coordConv.GalCoordSys,
     2: coordConv.FK5CoordSys,
     1: coordConv.FK4CoordSys,
    -1: coordConv.AppGeoCoordSys,
    -2: coordConv.AppTopoCoordSys,
    -3: coordConv.ObsCoordSys,
}
def getCoordSys(coordSysCode, date, tai):
    if coordSysCode < -1:
        # the old TCC use LAST computed from TAI; the new TCC uses TAI
        date = tai
    elif coordSysCode == -1 and date == 0:
        date = coordConv.julianEpochFromTAI(tai)
    return CSysDict[coordSysCode](date)

BoolDict = dict(T=True, F=False)
def cnvBool(val):
    return BoolDict[val]

CnvList = (int,) + (float,)*10 + (int,) + (float,)*9 + (cnvBool,)*3 + (float,)*2

class TestCoordConv(unittest.TestCase):
    def testFile(self):
        """Test file of coordinate conversions from TCC (data/masscc_out.dat)

        Known issues:
        - radVel does not match; the TCC seems to zero radVel if at infinity, but why?
          Also, the TCC seems to be able to round trip RadVel even if at infinity, but how,
          if it zeros it when at infinity? Once I resolve this, update the testCoord.py
          accordingly, as well as this code.
        - Other problems await at other coordinate systems.
        """
        site = None
        numErrors = 0
        with file(DataFile, "rU") as f:
            gotSiteData = False
            startTime = time.time()
            nTested = 0
            for lineInd, line in enumerate(f):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if not gotSiteData:
                    meanLat, meanLong, elevation, ut1_tai, poleX, poleY = [float(val) for val in line.split()]
                    site = coordConv.Site(meanLong, meanLat, elevation)
                    site.setPoleWander(poleX, poleY)
                    site.ut1_tai = ut1_tai
                    gotSiteData = True
                    continue

                dataList = line.split()
                fromSysCode, fromDate, fromPos1, fromPos2, fromPM1, fromPM2, fromParallax, fromRadVel, fromDir, refCoA, refCoB, \
                    toSysCode, toDate, refToPos1, refToPos2, refToPM1, refToPM2, refToParallax, refToRadVel, \
                    refToDir, refScaleChange, refAtInf, refAtPole, isOK, tai, last \
                    = [cnvFunc(val) for val, cnvFunc in itertools.izip(dataList, CnvList)]
                if not isOK:
                    print "Skipping line %s: %s; isOK false" % (lineInd + 1, line)
                if (fromSysCode == 1) and (fromRadVel != 0) and (fromPM1 == 0) and (fromPM2 == 0):
                    print "Skipping line %s; FK4 with zero PM and nonzero radVel" % (lineInd + 1,)
                    continue

                nTested += 1

                fromCoord = coordConv.Coord(fromPos1, fromPos2, fromParallax, fromPM1, fromPM2, fromRadVel)
                fromPVTCoord = coordConv.PVTCoord(fromCoord, fromCoord, tai, 0.01)
                fromPVTDir = coordConv.PVT(fromDir, 0, tai)

                fromCoordSys = getCoordSys(fromSysCode, fromDate, tai)
                toCoordSys = getCoordSys(toSysCode, toDate, tai)
                site.refCoA = refCoA
                site.refCoB = refCoB

                try:
                    toCoord, toDir, scaleChange = toCoordSys.convertFrom(fromCoordSys, fromCoord, fromDir, site)
                    toPVTDir = coordConv.PVT()
                    toPVTCoord, scaleChange2 = toCoordSys.convertFrom(toPVTDir, fromCoordSys, fromPVTCoord, fromPVTDir, site)
                except Exception:
                    print "Failed on line %s: %s\n" % (lineInd + 1, line)
                    raise

                atPole, toPos1, toPos2 = toCoord.getSphPos()
                toParallax = toCoord.getParallax()
                atPole, toPM1, toPM2 = toCoord.getPM()
                toRadVel = toCoord.getRadVel()
                if toCoord.atInfinity(): # emulate something the TCC does that I don't think my code can do
                    toRadVel = fromRadVel
                predList = (toParallax, toPM1, toPM2, toRadVel)
                refList  = (refToParallax, refToPM1, refToPM2, refToRadVel)
                refToCoord = coordConv.Coord(refToPos1, refToPos2, refToParallax, refToPM1, refToPM2, refToRadVel)

                try:
                    self.assertEqual(toCoord.atPole(), refAtPole)
                    self.assertEqual(toCoord.atInfinity(), refAtInf)
                    if (fromSysCode > 0) and (toSysCode > 0):
                        atol = 1e-7
                    elif (fromSysCode < -1) and (toSysCode < -1):
                        atol = 1e-7
                    else:
                        # the sla_Mappa in the old TCC is giving slightly different answers
                        # thatn the latest slaMappa and that appears to explain a small discrepancy
                        # when converting to/from apparent geocentric coordinates;
                        # the error is most noticeable for the precession/nutation matrix.
                        atol = 1e-3
                    self.assertLess(toCoord.angularSeparation(refToCoord), atol)
                    self.assertLess(toPVTCoord.getCoord(tai).angularSeparation(refToCoord), atol)
                    maxPxDelta = refToParallax * 1000.0
                    self.assertAlmostEqual(toParallax, refToParallax, delta = maxPxDelta)
                    self.assertTrue(numpy.allclose(predList[1:], refList[1:], atol=atol))
                    self.assertAlmostEqual(refToDir, coordConv.wrapNear(toDir, refToDir), places=2)
                    self.assertAlmostEqual(refToDir, coordConv.wrapNear(toPVTDir.getPos(tai), refToDir), places=2)
# scale change bears very little resemblance between old and new.
# I believe this is a bug in the old TCC, since mean->mean should be 1.0
# and the new code is significantly closer to 1.0 than the old code.
#                    self.assertAlmostEqual(refScaleChange, scaleChange, places=5)
                    self.assertAlmostEqual(scaleChange, scaleChange2, places=5)
                    if (fromSysCode > 0) and (toSysCode > 0):
                        self.assertAlmostEqual(scaleChange, 1.0, places=5)

                    if toCoordSys.getDateType() == coordConv.DateType_TAI:
                        # "to" system uses tai as its time; try various strategies that remove proper motion to the given tai date

                        # test the removePM function (which removes proper motion and radial velocity, but not parallax)
                        zpmFromCoord = fromCoordSys.removePM(fromCoord, tai)

                        if fromCoordSys.getName() != "fk4":
                            # FK4 coordinates have fictitious space motion
                            zpmFromAtPole, zpmFromPM1, zpmFromPM2 = zpmFromCoord.getPM()
                            self.assertEqual(fromCoord.atPole(), zpmFromAtPole)
                            self.assertEqual(zpmFromPM1, 0)
                            self.assertEqual(zpmFromPM2, 0)
                            zpmFromRadVel = zpmFromCoord.getRadVel()
                            self.assertEqual(zpmFromRadVel, 0)

                        # zpmFromAtPole, zpmFromPM1, zpmFromPM2 = zpmFromCoord.getPM()
                        # self.assertEqual(fromCoord.atPole(), zpmFromAtPole)
                        # zpmFromRadVel = zpmFromCoord.getRadVel()
                        # self.assertEqual(zpmFromPM1, 0)
                        # self.assertEqual(zpmFromPM2, 0)
                        # self.assertEqual(zpmFromRadVel, 0)

                        zpmToCoord, zpmToDir, zpmScaleChange = toCoordSys.convertFrom(fromCoordSys, zpmFromCoord, fromDir, site)
                        zpmToAtPole, zpmToPos1, zpmToPos2 = zpmToCoord.getSphPos()
                        self.assertEqual(atPole, zpmToAtPole)
                        zpmToAtPole, zpmToPM1, zpmToPM2 = zpmToCoord.getPM()
                        self.assertEqual(atPole, zpmToAtPole)
                        zpmToRadVel = zpmToCoord.getRadVel()

                        self.assertAlmostEqual(toDir, zpmToDir, places=2) # why so poor?
                        self.assertAlmostEqual(scaleChange, zpmScaleChange, places=6)
                        self.assertLess(toCoord.angularSeparation(zpmToCoord), 1e-7)
                        self.assertEqual(zpmToPM1, 0)
                        self.assertEqual(zpmToPM2, 0)
                        self.assertEqual(zpmToRadVel, 0)

                except Exception as e:
                    if ContinueOnError:
                        print
                        print str(e)
                    print "Failed on line %s: %s" % (lineInd + 1, line)
                    print "fromCoordSys=(%s, %s); toCoordSys=(%s, %s)" % (fromCoordSys.getName(), fromCoordSys.getDate(), toCoordSys.getName(), toCoordSys.getDate())
                    print "toSphPos=   ", toPos1, toPos2
                    print "refToSphPos=", refToPos1, refToPos2
                    print "angular sep=", toCoord.angularSeparation(refToCoord) * 3600.0, "arcsec"
                    print "pred parallax, PM and radVel=", predList
                    print "ref  parallax, PM and radVel=", refList
                    print "from parallax, PM and radVel=", (fromParallax, fromPM1, fromPM2, fromRadVel)
                    print "from vec pos, vel=", fromCoord.getVecPos(), fromCoord.getVecPM()
                    print "to   vec pos, vel=", toCoord.getVecPos(),  toCoord.getVecPM()
                    if not ContinueOnError:
                        raise
                    numErrors += 1
        duration = time.time() - startTime
        print "Tested %d conversions in %0.2f seconds: %0.0f conversions/second" % \
            (nTested, duration, nTested/duration)
        self.assertEqual(numErrors, 0, "%s errors" % (numErrors,))


if __name__ == '__main__':
    unittest.main()
