#pragma once

#include "Eigen/Dense"

#include "coordConv/coord.h"
#include "coordConv/site.h"

namespace coordConv {

    /**
    Convert apparent topocentric coordinates to observed (refracted apparent topocentric)

    @param[in] appTopoCoord  apparent topocentric coord
    @param[in] site  site information; refCoA and refCoB are read
    @return position in observed coordinates
    */
    Coord obsFromAppTopo(
        Coord const &appTopoCoord,
        Site const &site
    );

}
