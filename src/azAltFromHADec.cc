#include "coordConv/mathUtils.h"
#include "coordConv/azAltFromHADec.h"

namespace coordConv {

    void azAltFromHADec (
        Eigen::Vector3d &azAlt,
        Eigen::Vector3d const &haDec,
        double lat
    ) {
        const double sinLat = sind(lat);
        const double cosLat = cosd(lat);

        azAlt(0) = + sinLat * haDec(0) - cosLat * haDec(2);
        azAlt(1) = + haDec(1);
        azAlt(2) = + cosLat * haDec(0) + sinLat * haDec(2);
    }

}
