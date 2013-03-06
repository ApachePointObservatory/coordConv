#!/usr/bin/env python
import unittest
import math
import numpy
import coordConv

class TestCoordSys(unittest.TestCase):
    """Test some aspects of CoordSys and subclasses
    
    Note that coordinate conversions are tested by testCoordConvAgainstTCC
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


if __name__ == '__main__':
    unittest.main()
