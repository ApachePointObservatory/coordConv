#pragma once

#include "Eigen/Dense"

namespace coordConv {

    /**
    Convert Az/Alt position to HA/Dec

    @param[out] haDec  cartesian -HA, Dec (same units as azAlt)
    @param[in] azAlt  cartesian Az/Alt (any units)
        Sign convention: (1,0,0) is south, (0,1,0) is east.
    @param[in] lat  latitude (degrees)
    */
    void haDecFromAzAlt(
        Eigen::Vector3d &haDec,
        Eigen::Vector3d const &azAlt,
        double lat
    );

}
