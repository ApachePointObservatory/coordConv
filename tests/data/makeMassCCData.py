#!/usr/bin/env python
"""
Generate data for masscoordconv

Each line should contain:
  fromSysCode fromDate fromPos1 fromPos2 fromPM1 fromPM2 fromParlax fromRadVel fromDir refCoA refCoB toSysCode toDate

Subtleties:
- Do not convert from FK4 with zero PM and nonzero radVel; the TCC treats this as zero space motion,
  but the new conversion code does not (and could not easily do so without some ugly hacks).

Old coordinate system codes
    -9 = "GuideImage",
    -8 = "GuideProbe",
    -7 = "Rotator",
    -6 = "Instrument",

    -5 = "Mount",
    -4 = "Physical",

    -3 = "Observed",
    -2 = "AppTopo",
    -1 = "AppGeo",

    0 = "None",

    1 = "FK4",
    2 = "FK5",
    3 = "Gal",
    4 = "ICRS",
"""
def getDateIter(sysCode):
    """Return a list of dates given a coordinate system code
    """
    if sysCode > 0:
        # mean
        dateList = (1960, 2010)
    elif sysCode == -1:
        # apparent geocentric; 0 is relevant
        dateList = (0, 1999.5, 2015)
    else:
        # apparent topocentric or focal plane; 0 is the only reasonable choice
        dateList = (0,)
    return dateList

def getPosIter(fromSys):
    """Return an iterator over pos1, pos2
    """
    if fromSys > -3:
        return [
            (175.5, -84.3),
            (-10.5, -1.253),
            (60.35, 23.64),
            (-151.234, 84.43),
        ]
    elif fromSys == -3:
        return [
            (175.5, 5.0),
            (60.35, 23.64),
            (-151.234, 84.43),
        ]
    else:
        raise RuntimeError("Unknown fromSys=%s" % (fromSys,))

def getFromDirIter():
    return (-42.34, 65.43)

def getPMPxRVIter(sysCode):
    """Get an iterator over proper motion 1, proper motion 2, parallax and radial velocity given a coordinate system code
    """
    valList = [(0, 0, 0, 0)]
    if fromSysCode > 0:
        # from mean; parallax, etc. are relevant
        valList += [
            (5, -3, 7, 10),
        ]
    return valList

def getRefCoIter(fromSysCode, toSysCode):
    if fromSysCode > -2 and toSysCode > -2:
        # convert between mean, apparent geocentric or apparent topocentric coordinates; refraction coefficients not used
        return [(1.2e-2,  -1.3e-5)]
    else:
        # refraction coefficients will (probably) be used
        return [(1.2e-2,  -1.3e-5), (2.2e-2, -1.7e-5)]

coordSysList = (-3, -2, -1, 1, 2, 3, 4)

for fromSysInd, fromSysCode in enumerate(coordSysList):
    for fromDate in getDateIter(fromSysCode):
        for fromPos1, fromPos2 in getPosIter(fromSysCode):
            for fromPM1, fromPM2, fromParallax, fromRadVel in getPMPxRVIter(fromSysCode):
                for fromDir in getFromDirIter():
                    for toSysCode in coordSysList[fromSysInd:]:
                        for toDate in getDateIter(toSysCode):
                            for refCoA, refCoB in getRefCoIter(fromSysCode, toSysCode):
                                print fromSysCode, fromDate, fromPos1, fromPos2, \
                                    fromPM1, fromPM2, fromParallax, fromRadVel, fromDir, refCoA, refCoB, toSysCode, toDate
