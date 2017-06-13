from __future__ import absolute_import, division, print_function

from .coordConvLib import wrapCtr


def assertPVTsAlmostEqual(pvt1, pvt2, doWrap=False, posPlaces=7, velPlaces=7, tPlaces=7):
    """Assert that two PVTs are almost equal

    @param[in] pvt1  first coordConv.PVT
    @param[in] pvt2  second coordConv.PVT
    @param[in] doWrap  if True then wrap position difference to [-180, 180] before comparing
        (ignored for velocity and time comparison).
    @param[in] posPlaces  number of decimal places for position difference
    @param[in] velPlaces  number of decimal places for velocity difference
    @param[in] tPlaces  number of decimal places for time difference

    For the places arguments the associated value must be less than 10^-places
    """
    if abs(pvt1.t - pvt2.t) > 10.0**-tPlaces:
        raise AssertionError("%s.t != %s.t to within %s places" % (pvt1, pvt2, tPlaces))
    if abs(pvt1.vel - pvt2.vel) > 10.0**-velPlaces:
        raise AssertionError("%s.vel != %s.vel to within %s places" % (pvt1, pvt2, velPlaces))
    if doWrap:
        if abs(wrapCtr(pvt1.pos - pvt2.pos)) > 10.0**-posPlaces:
            raise AssertionError("%s.pos != %s.pos to within %s places (wrapped)" % (pvt1, pvt2, posPlaces))
    else:
        if abs(pvt1.pos - pvt2.pos) > 10.0**-posPlaces:
            raise AssertionError("%s.pos != %s.pos to within %s places" % (pvt1, pvt2, posPlaces))


def assertAnglesAlmostEqual(angle1, angle2, places=7):
    """Assert that two PVTs are almost equal

    @param[in] angle1  first angle
    @param[in] angle  second angle
    @param[in] places  number of decimal places for position difference
    @throw AssertionError if abs(angle1 - angle2) > 10.0**-places
    """
    if abs(wrapCtr(angle1 - angle2)) > 10.0**-places:
        raise AssertionError("%r != %r to within %s places (wrapped)" % (angle1, angle2, places))
