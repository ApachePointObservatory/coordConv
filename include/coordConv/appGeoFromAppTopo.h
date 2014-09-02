#pragma once

#include "Eigen/Dense"

#include "coordConv/coord.h"
#include "coordConv/site.h"

namespace coordConv {

    /**
    Convert apparent geocentric coordinates to apparent topocentric at a specified date

    @param[in] appTopoCoord  apparent topocentric coord at the specified TAI date
    @param[in] site  site information
    @param[in] tai  TAI date (MJD, sec)
    @return position in apparent topocentric coordinates at the specified TAI date
    */
    Coord appGeoFromAppTopo(
        Coord const &appTopoCoord,
        Site const &site,
        double tai
    );

}
