#pragma once

#include "Eigen/Dense"

#include "coordConv/coord.h"
#include "coordConv/site.h"

namespace coordConv {

    /**
    Convert apparent geocentric coordinates to apparent topocentric at a specified date

    @param[in] appGeoCoord  apparent geocentric coord at the specified TAI date
    @param[in] site  site information
    @param[in] tai  TAI date (MJD, sec)
    @return position in apparent geocentric coordinates at the specified TAI date
    */
    Coord appTopoFromAppGeo(
        Coord const &appGeoCoord,
        Site const &site,
        double tai
    );

}
