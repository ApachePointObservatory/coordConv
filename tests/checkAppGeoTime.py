#!/usr/bin/env python
from __future__ import absolute_import, division

import coordConv
"""
Determine how much delta-T is allowed before apparent geocentric data is too stale
"""
def runOne(fromCoordSys, coord, site, errDict):
    for toDate in (2000, 2010, 2020):
        refAppGeoCoordSys = coordConv.AppGeoCoordSys(toDate)
        refAppGeoCoord = refAppGeoCoordSys.convertFrom(fromCoordSys, coord, site)
        refFromCoord = fromCoordSys.convertFrom(refAppGeoCoordSys, coord, site)
        for dataAgeSec in (50, 100, 200, 400, 1000):
            for ageMult in (-1, 1):
                dataAgeYears = dataAgeSec * ageMult / (coordConv.SecPerDay * coordConv.DaysPerYear)
                compDate = toDate + dataAgeYears
                appGeoCoordSys = coordConv.AppGeoCoordSys(compDate)
                appGeoCoord = appGeoCoordSys.convertFrom(fromCoordSys, coord, site)
                fromCoord = fromCoordSys.convertFrom(appGeoCoordSys, coord, site)
                angSep = max(refAppGeoCoord.angularSeparation(appGeoCoord), refFromCoord.angularSeparation(fromCoord))
                dictKey = (dataAgeSec, ageMult)
                oldAngSep = errDict.get(dictKey)
                if oldAngSep is None or angSep > oldAngSep:
                    errDict[dictKey] = angSep
    
errDict = {}
site = coordConv.Site(-105.822616, 32.780988, 2.788)
for fromCoordSys in (coordConv.ICRSCoordSys(),):
    for equatAng in (0, 45):
        for polarAng in (-85, 0, 30, 89):
            coord = coordConv.Coord(equatAng, polarAng)
            runOne(fromCoordSys, coord, site, errDict)

print "Results (in arcsec) for fromCoordSys=", fromCoordSys.getName(), fromCoordSys.getDate(), "; coord=", coord.getSphPos()
for ageAndMult in sorted(errDict.keys()):
    print "%s: %0.4f" % (ageAndMult, errDict[ageAndMult] * 3600)

# based on these results, 200 seconds gives an error < 0.001 arcseconds,
# and it is not so sensitive that it needs to be a constructor parameter
# localhost$ tests/checkAppGeoTime.py 
# Results (in arcsec) for fromCoordSys= icrs 2000.0 ; coord= [False, 44.99999999999999, 29.999999999999996]
# (50, -1): 0.0002
# (50, 1): 0.0002
# (100, -1): 0.0004
# (100, 1): 0.0004
# (200, -1): 0.0008
# (200, 1): 0.0008
# (400, -1): 0.0017
# (400, 1): 0.0017
# (1000, -1): 0.0042
# (1000, 1): 0.0042
