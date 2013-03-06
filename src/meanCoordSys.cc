#include "coordConv/physConst.h"
#include "coordConv/time.h"
#include "coordConv/coordSys.h"

namespace coordConv {

    MeanCoordSys::MeanCoordSys(std::string const &name, double date)
    :
        CoordSys(name, true, date)
    {}
    
    // this works for all mean systems that use Julian epoch (all but FK4)
    double MeanCoordSys::dateFromTAI(double tai) const {
        return julianEpochFromMJDSec(tai + TT_TAI);
    }

    Coord MeanCoordSys::removePM(Coord const &coord, double tai) const {
        double epoch = dateFromTAI(tai);

        Eigen::Vector3d corrPos = coord.getVecPos() + ((epoch - this->_date) * coord.getVecVel());
        return Coord(corrPos);
    }

}
