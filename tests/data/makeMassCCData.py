#!/usr/bin/env python
"""
Generate data for masscoordconv

Each line should contain:
  fromSysCode fromDate fromPos1 fromPos2 fromPM1 fromPM2 fromParlax fromRadVel fromDir refCoA refCoB toSysCode toDate

Subtleties:
- Do not convert from FK4 with zero PM and nonzero radVel; the TCC treats this as zero space motion,
  but the new conversion code does not (and could not easily do so without some ugly hacks).
- Do not try multiple values of fromDir yet; wait until we have code that handles it
"""

coordSysList = (4, 3, 2, 1, -1, -2, -3)
for fromSysCode in coordSysList:
    if fromSysCode > 0:
        # from mean
        fromDateList = (1960, 2010)
        fromPM1List = (0, 5)
        fromPM2List = (0, -3)
        fromParlaxList = (0, 7)
        fromRadVelList = (0, 10)
    else:
        # from apparent; use current date (since the TCC is not accurate otherwise)
        if fromSysCode == -1:
            fromDateList = (0, 1999.5, 2012)
        else:
            fromDateList = [0]
        fromPM1List = [0]
        fromPM2List = [0]
        fromParlaxList = [0]
        fromRadVelList = [0]
    for fromDate in fromDateList:
        for toSysCode in coordSysList:
            if toSysCode > 0:
                # to mean
                toDateList = (1940, 2010)
            elif toSysCode == -1:
                # to apparent; use current date (since the TCC is not accurate otherwise)
                toDateList = (0, 2010, 2012.5)
            else:
                # apparent topocentric or observed
                toDateList = [0]

            if (fromSysCode < 0) or (toSysCode < 0):
                # from or to apparent; refCo used
                refCoAList = [ 1.2e-2,  2e-2]
                refCoBList = [-1.3e-5, -3e-5]
            else:
                # mean to mean; refCo ignored
                refCoAList = [ 1.2e-2]
                refCoBList = [-1.3e-5]
            for toDate in toDateList:
                for fromPos1 in fromPM1List:
                    for fromPos2 in fromPM2List:
                        for fromPM1 in fromPM1List:
                            for fromPM2 in fromPM2List:
                                if (fromSysCode == 1) and (fromPM1 == 0) and (fromPM2 == 0):
                                    fromRadVelList = [0] # see Subtleties above
                                for fromParlax in fromParlaxList:
                                    for fromRadVel in fromRadVelList:
                                        for fromDir in (0,): # (0, 45, -90) use a list once the new code supports dir and scale
                                            for refCoA in refCoAList:
                                                for refCoB in refCoBList:
                                                    print fromSysCode, fromDate, fromPos1, fromPos2, fromPM1, fromPM2, \
                                                        fromParlax, fromRadVel, fromDir, refCoA, refCoB, toSysCode, toDate
