#include <stdexcept>
#include "Eigen/Dense"
#include "coordConv/physConst.h"
#include "coordConv/time.h"
#include "coordConv/coordSys.h"

namespace coordConv {

    MeanCoordSys::MeanCoordSys(std::string const &name, double date, DateTypeEnum dateType)
    :
        CoordSys(name, date, dateType, true, true)
    {}
    
    // this works for all mean systems that use Julian epoch (all but FK4)
    double MeanCoordSys::dateFromTAI(double tai) const {
        return julianEpochFromMJDSec(tai + TT_TAI);
    }

    Coord MeanCoordSys::removePM(Coord const &coord, double tai) const {
        // convert to this coord at tai date, zero velocity and convert back;
        // this is fancier than just adding vecPM to vecPos, but handles fictitious proper motion (e.g. FK4).
        if ((coord.getVecPos().array() == 0.0).all()) {
            // no proper motion to correct; return the coord unchanged
            return coord;
        }

        double dateAtTAI = dateFromTAI(tai);
        CoordSys::Ptr coordSysAtTAIPtr = this->clone(dateAtTAI);
        Site site(10, 10, 10); // values are irrelevant for mean to mean coordinate conversions
        Coord coordAtTAI = coordSysAtTAIPtr->convertFrom(*this, coord, site);
        Eigen::Vector3d posAtTAI = coordAtTAI.getVecPos();
        Coord zpmCoordAtTAI(posAtTAI);
        return this->convertFrom(*coordSysAtTAIPtr, zpmCoordAtTAI, site);
    }

}
