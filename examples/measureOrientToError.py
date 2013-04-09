#!/usr/bin/env python
import math
import numpy
import coordConv

def measureOrientationToErr():
    """Measure error in Coord.orientationTo for tiny offsets
    
    For best results modify Coord::orientationTo to not short-circuit tiny values of sinVal, cosVal
    """
    errDict = {} # dict of distArcSec: (max orient error, max(sinVal, cosVal))
    for fromPolarAng in (-89, -72, -45.0, -30, 0.01, 12.5, 31, 47, 56, 68, 89):
        for fromEquatAng in (0, 41.0): # should not matter
            fromCoord = coordConv.Coord(fromEquatAng, fromPolarAng)
            for fromOrient in (-90, -72, -45.0, -30, 0.01, 0, 12.5, 31, 45, 56, 68, 89):
                # offset by such small distances that toOrient=fromOrient
                # by 1e-6 the error is starting to go up, indicating that the approximation is failing
                for distArcSec in (5e-5, 1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 1e-3, 1e-2):
                    toCoord, dumToOrient = fromCoord.offset(fromOrient, distArcSec / 3600)
                    toOrient = coordConv.wrapCtr(180 + toCoord.orientationTo(fromCoord))

                    # expect fromOrient = toOrient
                    newErrArcSec = abs(toOrient - fromOrient) * 3600
                    oldErrArcSec = errDict.get(distArcSec, (0, 0))[0]
                    if newErrArcSec > oldErrArcSec:
                        toU = toCoord.getVecPos()
                        toU /= numpy.linalg.norm(toU)
                        fromU = fromCoord.getVecPos()
                        fromU /= numpy.linalg.norm(fromU)
                        sinVal = (toU[1] * fromU[0]) - (toU[0] * fromU[1])
                        cosVal = (toU[2] * ((fromU[0] * fromU[0]) + (fromU[1] * fromU[1]))) \
                                  - (fromU[2] * ((toU[0] * fromU[0]) + (toU[1] * fromU[1])))
                        errDict[distArcSec] = (newErrArcSec, max(sinVal, cosVal))
    distKeys = sorted(errDict.keys())
    print "Note that by 0.01 arcsec the error creeps up, indicating that the approximation toOrient=fromOrient is failing"
    for distArcSec in distKeys:
        errArcSec, maxSinCos = errDict[distArcSec]
        print "distArcSec=%0.3g arcsec, maxErr=%0.3g arcsec, maxSin/Cos=%0.3g" % (distArcSec, errArcSec, maxSinCos)

if __name__ == "__main__":
    measureOrientationToErr()
