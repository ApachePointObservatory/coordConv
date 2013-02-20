#pragma once

#include "Eigen/Dense"

namespace coordConv {

    /**
    Convert HA/Dec position to Alt/Az.

    @param[out] azAlt cartesian Az/Alt (same units as haDec);
        Sign convention: (1,0,0) is south, (0,1,0) is east.
    @param[in] haDec cartesian -HA, Dec (any units)
    @param[in] lat latitude (degrees)
    */
    void azAltFromHADec (
        Eigen::Vector3d &azAlt,
        Eigen::Vector3d const &haDec,
        double lat
    );

}
