#pragma once

#include "Eigen/Dense"

#include "coordConv/coord.h"
#include "coordConv/site.h"

namespace coordConv {

    /**
    Convert observed coordinates (refracted apparent topocentric) to apparent topocentric coordinates

    @param[in] obsCoord  observed (refracted apparent topocentric) coord
    @param[in] site  site information; refCoA and refCoB are read
    @return position in observed coordinates
    */
    Coord appTopoFromObs(
        Coord const &obsCoord,
        Site const &site
    );

}
