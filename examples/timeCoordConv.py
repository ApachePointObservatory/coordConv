#!/usr/bin/env python
import time
import coordConv
"""
Measure the time required to perform coordinate conversions relevant to telescope tracking

Also measure the time required to recompute cached apparent geocentric data.
"""
_TimeTupleJ2000 = (2000, 1, 1, 12, 0, 0, 5, 1, 0)

def utcFromPySec(pySec=None):
    """Returns the UTC (MJD) corresponding to the supplied python time, or now if none.
    """
    global _TimeTupleJ2000

    if pySec == None:
        pySec = time.time()
    
    # python time (in seconds) corresponding to 2000-01-01 00:00:00
    # this is probably constant, but there's some chance
    # that on some computer systems it varies with daylights savings time
    pySecJ2000 = time.mktime(_TimeTupleJ2000) - time.timezone
    
    return coordConv.MJDJ2000 + ((pySec - pySecJ2000) / coordConv.SecPerDay)

def makeSite():
    """Create a realistic site object (for the APO 3.5m)"""
    site = coordConv.Site(-105.822616, 32.780988, 2.788)
    site.setPoleWander(0.89e-5,  0.92e-4)
    site.ut1_tai = -34.782
    site.refCoA =  1.2e-2 
    site.refCoB = -1.3e-5
    return site
    
def timeFK5ToAppTopo(alt, niter):
    """Time FK5 to apparent geocentric at the specified altitude
    
    @param[in] alt  initial altitude (degrees)
    @param[in] niter  number of iterations
    
    The coordinate conversion code presently does not have a special branch
    for zero proper motion, so no attempt is made to provide proper motion.
    
    Use the approximation that TAI = UTC, which is plenty close enough for timing.
    
    Increment TAI by 0.1 seconds per iteration, which is reasonable for a control loop.
    """
    fk5Sys = coordConv.FK5CoordSys(1980)
    currTAI = utcFromPySec(time.time())
    appTopoSys = coordConv.AppTopoCoordSys(currTAI)
    initialAppTopoCoord = coordConv.Coord(120, alt)
    site = makeSite()
    fk5Coord = fk5Sys.convertFrom(appTopoSys, initialAppTopoCoord, site)
    startTime = time.time()
    for i in range(niter):
        tai = currTAI + (i * 0.1)
        appTopoSys.setDate(tai)
        appTopoCoord, dir, scaleChange = appTopoSys.convertFrom(fk5Sys, fk5Coord, 5.0, site)
    duration = time.time() - startTime
    print "FK5 to AppTopo: %0.1f conversions/second (%d conversions in %0.2f sec) at alt=%0.1f" % (niter/duration, niter, duration, alt)
    
def timeAppTopoToFK5(alt, niter):
    """Time apparent topocentric to FK5 at the specified altitude
    
    @param[in] alt  initial altitude (degrees)
    @param[in] niter  number of iterations
    
    The coordinate conversion code presently does not have a special branch
    for zero proper motion, so no attempt is made to provide proper motion.
    
    Use the approximation that TAI = UTC, which is plenty close enough for timing.
    
    Increment TAI by 0.1 seconds per iteration, which is reasonable for a control loop.
    """
    fk5Sys = coordConv.FK5CoordSys(1980)
    currTAI = utcFromPySec(time.time())
    appTopoSys = coordConv.AppTopoCoordSys(currTAI)
    appTopoCoord = coordConv.Coord(120, alt)
    site = makeSite()
    startTime = time.time()
    for i in range(niter):
        tai = currTAI + (i * 0.1)
        appTopoSys.setDate(tai)
        fk5Coord, dir, scaleChange = fk5Sys.convertFrom(appTopoSys, appTopoCoord, 5.0, site)
    duration = time.time() - startTime
    print "AppTopo To FK5: %0.1f conversions/second (%d conversions in %0.2f sec) at alt=%0.1f" % (niter/duration, niter, duration, alt)

def timeAppGeoData(niter):
    """Time computation of apparent geocentric data
    
    Space the time out far enough for a cache miss; I'm using a year
    """
    epoch = 1950
    appGeoSys = coordConv.AppGeoCoordSys(epoch)
    startTime = time.time()
    for i in range(niter):
        epoch += 0.1
        if epoch > 2500:
            epoch = 1950
        appGeoSys.setDate(epoch)
    duration = time.time() - startTime
    print "AppGeo data: %0.1f computations/second (%d computations in %0.2f sec)" % (niter/duration, niter, duration)

if __name__ == '__main__':
    for alt in (0, 5, 45):
        timeFK5ToAppTopo(alt, 10000)
    print
    for alt in (0, 5, 45):
        timeAppTopoToFK5(alt, 10000)
    print
    timeAppGeoData(10000)
