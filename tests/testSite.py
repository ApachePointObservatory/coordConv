#!/usr/bin/env python
import math
import unittest
import numpy
import coordConv


class TestSite(unittest.TestCase):
    def testBasics(self):
        """Test basic functionality of Site
        """
        for long in (54.3, -134.5):
            for lat in (-23.4, 12.3):
                for elev in (2789, 6543):
                    site = coordConv.Site(long, lat, elev)
                    self.assertEqual(long, site.meanLong)
                    self.assertEqual(lat, site.meanLat)
                    self.assertEqual(elev, site.elev)
                    self.assertAlmostEqual(long, site.corrLong)
                    self.assertAlmostEqual(lat, site.corrLat)
                    site.setPoleWander(0, 0)
                    self.assertAlmostEqual(long, site.corrLong)
                    self.assertAlmostEqual(lat, site.corrLat)
                    for xArcsec in (-0.3, 0, 0.3):
                        x = xArcsec / 3600.0
                        for yArcsec in (-0.3, 0, 0.3):
                            y = yArcsec / 3600.0
                            site.setPoleWander(x, y)
                            # the following is an approximation given in the header of slaPolmo
                            approxCorrLong = long +  (x * coordConv.cosd(long)) - (y * coordConv.sind(long))
                            approxCorrLat  = lat + (((x * coordConv.sind(long)) + (y * coordConv.cosd(long))) * coordConv.tand(lat))
                            self.assertAlmostEqual(approxCorrLong, site.corrLong, 3)
                            self.assertAlmostEqual(approxCorrLat, site.corrLat, 3)
    
    def testError(self):
        """Test that constructing with invalid latitude raises an exception
        """
        for long in (54.3, -134.5):
            for elev in (2789, 6543):
                for baseLat in (-90, 90):
                    for dLat in (-1, -0.00001, 0, 0.00001, 1):
                        lat = baseLat + dLat
                        if -90 <= lat <= 90:
                            # make sure Site can be constructed
                            coordConv.Site(long, lat, elev)
                        else:
                            # make sure Site construction raises an exception
                            self.assertRaises(Exception, coordConv.Site, long, lat, elev)
                    
        

if __name__ == '__main__':
    unittest.main()
