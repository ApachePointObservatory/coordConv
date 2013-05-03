#!/usr/bin/env python
import unittest
import math
import numpy
import coordConv

class TestCoordSys(unittest.TestCase):
    """Test some aspects of CoordSys and subclasses
    
    Note that most coordinate conversions are tested by testCoordConvAgainstTCC
    """
    def testIsMean(self):
        """Test isMean
        """
        for cls in (coordConv.ICRSCoordSys, coordConv.FK5CoordSys, coordConv.FK4CoordSys, coordConv.GalCoordSys):
            coordSys = cls(2000)
            self.assertTrue(coordSys.isMean())
        
        for cls in (coordConv.AppGeoCoordSys, coordConv.AppTopoCoordSys, coordConv.ObsCoordSys):
            coordSys = cls(2000)
            self.assertFalse(coordSys.isMean())
    
    def testDateFromTAI(self):
        for cls in (coordConv.ICRSCoordSys, coordConv.FK5CoordSys, coordConv.GalCoordSys, coordConv.AppGeoCoordSys):
            coordSys = cls(2000)
            for tai in (4232.89, 20000.32, 56350.03, 74222.9):
                predDate = coordConv.julianEpochFromMJDSec(coordConv.TT_TAI + tai)
                self.assertAlmostEqual(predDate, coordSys.dateFromTAI(tai))
        
        coordSys = coordConv.FK4CoordSys(2000)
        for tai in (4232.89, 20000.32, 56350.03, 74222.9):
            ttDays = (tai + coordConv.TT_TAI) / coordConv.SecPerDay
            predDate = 1900.0 + ((ttDays - 15019.81352 ) / 365.242198781)
            self.assertAlmostEqual(predDate, coordSys.dateFromTAI(tai), places=5)
        
        for cls in (coordConv.AppTopoCoordSys, coordConv.ObsCoordSys):
            for tai in (4232.89, 20000.32, 56350.03, 74222.9):
                coordSys = cls(2000)
                self.assertAlmostEqual(tai, coordSys.dateFromTAI(tai))
    
    def testNullConversion(self):
        site = coordConv.Site(-105.822616, 32.780988, 2788)
        site.setPoleWander(1.1e-4, -0.5e-4)
        site.ut1_tai = -2e-8
        site.refCoA =  1.2e-2
        site.refCoB = -1.3e-5

        maxRoundTripErr = 0
        for cls in (coordConv.ICRSCoordSys, coordConv.FK5CoordSys, coordConv.GalCoordSys, coordConv.AppGeoCoordSys):
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
        for cls in (coordConv.AppTopoCoordSys, coordConv.ObsCoordSys):
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
        """Test makeCoordSys
        """
        for csysName in ("icrs", "fk5", "fk4", "gal", "appgeo", "apptopo", "obs", "none", "mount"):
            predIsMean = csysName in set(("icrs", "fk5", "fk4", "gal"))
            for date in (1000.5, 2001):
                csys = coordConv.makeCoordSys(csysName, date)
                self.assertEqual(csys.getDate(), date)
                self.assertEqual(csys.isMean(), predIsMean)
                self.assertEqual(csys.getName(), csysName)
                
    def testNoneAndMountCoordSys(self):
        """Test that conversions to and from NoneCoordSys and MountCoordSys yield a null result
        """
        for nullSys in (
            coordConv.NoneCoordSys(),
            coordConv.MountCoordSys()
        ):
            site = coordConv.Site(-105.822616, 32.780988, 2788)
            fromCoord = coordConv.Coord(10, 30)
            for csysName in ("icrs", "fk5", "fk4", "gal", "appgeo", "apptopo", "obs", "none", "mount"):
                otherSys = coordConv.makeCoordSys(csysName, 2001)
                toCoord = nullSys.convertFrom(otherSys, fromCoord, site)
                self.assertFalse(toCoord.isfinite())
                
                toCoord = otherSys.convertFrom(nullSys, fromCoord, site)
                self.assertFalse(toCoord.isfinite())
                
                
        

if __name__ == '__main__':
    unittest.main()
