#pragma once

#include "coordConv/site.h"

namespace coordConv {

    /**
    Compute local apparent sidereal time from TAI
    
    @param[in] tai  universal time (MJD, seconds)
    @param[in] site  site information (a coordConv::Site)
        read fields: ut1_tai, longitude
    @return Local mean sidereal time, in degrees, in range [0, 360)
    */
    double lastFromTAI(double tai, Site const &site);
    
    /**
    Convert MJD seconds to Julian epoch

    @param[in] mjdSec  MJD date (sec)
    */
    double julianEpochFromMJDSec(double mjdSec);
    
    /**
    Convert Julian epoch to MJD seconds

    @param[in] julianEpoch  date as a Julian epoch (years)
    */
    double mjdSecFromJulianEpoch(double julianEpoch);

}
