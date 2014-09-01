#!/usr/bin/env python
from __future__ import absolute_import, division

import unittest

import numpy

import coordConv

MeanSysList = (coordConv.ICRSCoordSys, coordConv.FK5CoordSys, coordConv.FK4CoordSys, coordConv.GalCoordSys)
JulianSysList = (coordConv.ICRSCoordSys, coordConv.FK5CoordSys, coordConv.GalCoordSys, coordConv.AppGeoCoordSys)
MeanNameList = ("icrs", "fk5", "fk4", "gal")
AzAltSysList = (coordConv.AppTopoCoordSys, coordConv.ObsCoordSys)
AzAltNameList = ("apptopo", "obs")
AppSysList = (coordConv.AppGeoCoordSys,) + AzAltSysList
AppNameList = ("appgeo",) + AzAltNameList
FullSysList = MeanSysList + AppSysList + (coordConv.NoneCoordSys,)
FullNameList = MeanNameList + AppNameList + ("none",)

class TestCoordSys(unittest.TestCase):
    """Test some aspects of CoordSys and subclasses
    
    Note that most coordinate conversions are tested by testCoordConvAgainstTCC
    """
    def testIsMean(self):
        """Test isMean
        """
        meanNameSet = set(MeanNameList)
        for cls in FullSysList:
            coordSys = cls(2000)
            if coordSys.getName() in meanNameSet:
                self.assertTrue(coordSys.isMean())
            else:
                self.assertFalse(coordSys.isMean())
    
    def testDateFromTAI(self):
        """Test dateFromTAI
        """
        for cls in JulianSysList:
            coordSys = cls(2000)
            for tai in (4232.89, 20000.32, 56350.03, 74222.9):
                predDate = coordConv.julianEpochFromMJDSec(coordConv.TT_TAI + tai)
                self.assertAlmostEqual(predDate, coordSys.dateFromTAI(tai))
        
        coordSys = coordConv.FK4CoordSys(2000)
        for tai in (4232.89, 20000.32, 56350.03, 74222.9):
            ttDays = (tai + coordConv.TT_TAI) / coordConv.SecPerDay
            predDate = 1900.0 + ((ttDays - 15019.81352 ) / 365.242198781)
            self.assertAlmostEqual(predDate, coordSys.dateFromTAI(tai), places=5)
        
        for cls in AzAltSysList:
            for tai in (4232.89, 20000.32, 56350.03, 74222.9):
                coordSys = cls(2000)
                self.assertAlmostEqual(tai, coordSys.dateFromTAI(tai))

    def testEquality(self):
        """Test operator== and operator!=
        """
        def csysIter():
            for csysClass in FullSysList:
                for date in (2002, 3001):
                    for isCurrent in (False, True):
                        if isCurrent:
                            if csysClass in (coordConv.FK4CoordSys, coordConv.FK5CoordSys):
                                continue
                            csys = csysClass(0)
                            csys.setCurrDate(date)
                            yield csys
                        else:
                            yield csysClass(date)

        for csys1 in csysIter():
            for csys2 in csysIter():
                if (csys1.getName() == csys2.getName()) \
                    and (csys1.getDate() == csys2.getDate()) \
                    and (csys1.isCurrent() == csys2.isCurrent()):
                    self.assertTrue(csys1 == csys2)
                else:
                    self.assertTrue(csys1 != csys2)
                self.assertNotEqual(csys1 == csys2, csys1 != csys2)
    
    def testNullConversion(self):
        """Test round trip conversion
        """
        site = coordConv.Site(-105.822616, 32.780988, 2788)
        site.setPoleWander(1.1e-4, -0.5e-4)
        site.ut1_tai = -2e-8
        site.refCoA =  1.2e-2
        site.refCoB = -1.3e-5

        maxRoundTripErr = 0
        for cls in JulianSysList:
            for date in (1975, 2012):
                coordSys = cls(date)
                for equatAng in (100, -45):
                    for polAng in (-20, 25, 89):
                        for dist in (0, 0.01, 1):
                            for pmRA in (0, -3):
                                for pmDec in (0, 5):
                                    for radVel in (0, 7):
                                        fromCoord = coordConv.Coord(equatAng, polAng, dist, pmRA, pmDec, radVel)
                                        fk5Coord = coordSys.toFK5J2000(fromCoord, site)
                                        toCoord = coordSys.fromFK5J2000(fk5Coord, site)
                                        roundTripErr = toCoord.angularSeparation(fromCoord)
                                        maxRoundTripErr = max(roundTripErr, maxRoundTripErr)
                                        self.assertLess(roundTripErr, 1e-8)
        print "maxRoundTripErr for mean and app. geo. coordinate systems =", maxRoundTripErr, "deg"

        maxRoundTripErr = 0
        for cls in AzAltSysList:
            for date in (4842765000, 4872765000):
                coordSys = cls(date)
                for equatAng in (100, -45):
                    for polAng in (0, 25, 89):
                        for dist in (0, 0.01, 1):
                            fromCoord = coordConv.Coord(equatAng, polAng, dist)
                            fk5Coord = coordSys.toFK5J2000(fromCoord, site)
                            toCoord = coordSys.fromFK5J2000(fk5Coord, site)
                            roundTripErr = toCoord.angularSeparation(fromCoord)
                            maxRoundTripErr = max(roundTripErr, maxRoundTripErr)
                            self.assertLess(roundTripErr, 1e-7)
        print "maxRoundTripErr for app. topo and observed coordinate systems =", maxRoundTripErr, "deg"
    
    def testMakeCoordSys(self):
        """Test makeCoordSys and clone
        """
        for csysName in FullNameList:
            predIsMean = csysName in set(MeanNameList)
            for date in (0, 1000.5, 2001):
                if date == 0 and csysName in ("fk4", "fk5"):
                    self.assertRaises(Exception, coordConv.makeCoordSys, csysName, date)
                else:
                    csys = coordConv.makeCoordSys(csysName, date)
                    self.assertEqual(csys.getDate(), date)
                    self.assertEqual(csys.isMean(), predIsMean)
                    self.assertEqual(csys.getName(), csysName)
                    self.assertEqual(csys.canConvert(), csysName != "none")
                    if date == 0:
                        self.assertTrue(csys.isCurrent())
                    else:
                        self.assertFalse(csys.isCurrent())
                    
                    csysClone = csys.clone()
                    self.assertEqual(csys.getDate(), csysClone.getDate())
                    self.assertEqual(csys.isMean(),  csysClone.isMean())
                    self.assertEqual(csys.isCurrent(), csysClone.isCurrent())
                    self.assertEqual(csys.getName(), csysClone.getName())

    def testsetCurrDate(self):
        """Test setCurrDate method
        """
        currTAI = 123456789.12
        for csysName in FullNameList:
            for date in (0, 1000.5):
                if date == 0 and csysName in ("fk4", "fk5"):
                    self.assertRaises(Exception, coordConv.makeCoordSys, csysName, date)
                else:
                    csys = coordConv.makeCoordSys(csysName, date)
                    if not csys.isCurrent():
                        self.assertRaises(Exception, csys.setCurrDate, currTAI)
                    else:
                        csys.setCurrDate(csys.dateFromTAI(currTAI))
                        self.assertEqual(csys.getDate(), 0)
                        self.assertEqual(csys.getDate(False), csys.dateFromTAI(currTAI))

    def testAppGeoCache(self):
        """Test that AppGeoCoordSys cache is updated when expected
        """
        agSys = coordConv.AppGeoCoordSys()
        self.assertEqual(agSys.getDate(), 0)
        self.assertEqual(agSys.getDate(False), 0)
        self.assertFalse(numpy.isfinite(agSys.getCacheDate()))

        startDate = 2020.12345
        maxCacheAge = agSys.getMaxAge()
        maxDDate = agSys.getMaxDDate()
        nextDate = startDate + (maxCacheAge * 5) # well away from startDate
        predCacheDate = numpy.nan
        prevDate = startDate
        for i, (date, shouldUpdateCache) in enumerate((
            (startDate, True),
            (startDate + maxCacheAge - (maxDDate * 2), False),
            (startDate + maxCacheAge + (maxDDate * 0.01), True),
            (startDate + (maxCacheAge * 2), False),
            (nextDate, True),
            (nextDate + maxCacheAge - (maxDDate * 0.5), False),
            (nextDate + maxCacheAge + (maxDDate * 0.01), False), # within maxDDate of prev. update
            (nextDate + maxCacheAge + (maxDDate * 2), True),
        )):
            try:
                agSys.setCurrDate(date)
                if shouldUpdateCache:
                    predCacheDate = date
                self.assertAlmostEqual(agSys.getCacheDate(), predCacheDate)
                self.assertEqual(agSys.getDate(), 0)
                self.assertAlmostEqual(agSys.getDate(False), date)
            except Exception:
                print "failed on %d when date-prevDate=%0.11f; date-cacheDate=%0.11f, maxCacheAge=%0.11f, maxDDate=%0.11f" % \
                    (i, date-prevDate, date-agSys.getCacheDate(), maxCacheAge, maxDDate)
                raise
            prevDate = date

    def testDefaultConstructor(self):
        """Test default constructor
        """
        for csysClass in FullSysList:
            csys = csysClass()
            predDate = dict(
                fk4 = 1950,
                fk5 = 2000,
            ).get(csys.getName(), 0)
            self.assertEqual(csys.getDate(), predDate)
            self.assertEqual(csys.isCurrent(), predDate == 0)

    def testCopyConstructor(self):
        """Test copy constructor
        """
        for csysClass in FullSysList:
            for date in (1000.5, 2001):
                csys = csysClass(date)
                self.assertEqual(csys.getDate(), date)
                
                csysCopy = csysClass(csys)
                self.assertEqual(csys.getDate(), csysCopy.getDate())
                self.assertEqual(csys.getName(), csysCopy.getName())
                
    def testNoneAndOtherCoordSys(self):
        """Test that conversions to and from NoneCoordSys and OtherCoordSys fail
        """
        for nullSys in (
            coordConv.NoneCoordSys(),
            coordConv.OtherCoordSys("foo"),
        ):
            self.assertFalse(nullSys.canConvert())
            site = coordConv.Site(-105.822616, 32.780988, 2788)
            fromCoord = coordConv.Coord(10, 30)
            for csysName in FullNameList:
                otherSys = coordConv.makeCoordSys(csysName, 2001)
                self.assertRaises(Exception, nullSys.convertFrom, otherSys, fromCoord, site)
                self.assertRaises(Exception, otherSys.convertFrom, nullSys, fromCoord, site)
    
    def testInvalidCoordSys(self):
        """Test that makeCoordSys raises ValueError for invalid coordinate system name
        """
        for csysName in ("foo", "bar", "icrsx", "null"):
            self.assertRaises(ValueError, coordConv.makeCoordSys, csysName, 0)

if __name__ == '__main__':
    unittest.main()
