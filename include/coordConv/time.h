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
    Convert TAI (MJD seconds) to Julian epoch

    @param[in] mjdSec  MJD date (sec)
    */
    double julianEpochFromTAI(double mjdSec);
    
    /**
    Convert Julian epoch to TAI (MJD seconds)

    @param[in] julianEpoch  date as a Julian epoch (years)
    */
    double taiFromJulianEpoch(double julianEpoch);

    /**
    Convert TAI (MJD seconds) to Besselian epoch

    @param[in] mjdSec  MJD date (sec)
    */
    double besselianEpochFromTAI(double mjdSec);
    
    /**
    Convert Besselian epoch to TAI (MJD seconds)

    @param[in] besselianEpoch  date as a Besselian epoch (years)
    */
    double taiFromBesselianEpoch(double besselianEpoch);

}
