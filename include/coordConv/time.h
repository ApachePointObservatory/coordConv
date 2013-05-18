#pragma once

#include "coordConv/site.h"

namespace coordConv {

    /**
    Compute local apparent sidereal time from TAI
    
    @param[in] tai: universal time (MJD, seconds)
    @param[in] ut1_tai: UT1-TAI (seconds)
    @param[in] longitude: longitude east (degrees)
    @return Local mean sidereal time, in degrees, in range [0, 360)
    */
    double lastFromTAI(double tai, Site const &site);
    
    /**
    Convert MJD seconds to Julian epoch
    */
    double julianEpochFromMJDSec(double mjdSec);
    
    /**
    Convert Julian epoch to MJD seconds
    */
    double mjdSecFromJulianEpoch(double julianEpoch);

}
