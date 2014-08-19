#include "coordConv/mathUtils.h"
#include "coordConv/haDecFromAzAlt.h"

namespace coordConv {

    /**
    Compute HA/Dec position from alt/az position

    @param[out] haDec  cartesian -HA, Dec (same units as azAlt)
    @param[in] lat  latitude (degrees)
    @param[in] azAlt  cartesian Az/Alt (any units)
        Sign convention: (1,0,0) is south, (0,1,0) is east.
    */
    void haDecFromAzAlt(
        Eigen::Vector3d &haDec,
        Eigen::Vector3d const &azAlt,
        double lat
    ) {
        double const sinLat = sind(lat);
        double const cosLat = cosd(lat);

        haDec <<
            + (sinLat * azAlt(0)) + (cosLat * azAlt(2)),
            +  azAlt(1),
            - (cosLat * azAlt(0)) + (sinLat * azAlt(2));
    }

}
