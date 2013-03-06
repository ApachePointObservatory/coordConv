#include <stdexcept>
#include <sstream>
#include "coordConv/mathUtils.h"
#include "coordConv/coordSys.h"

namespace coordConv {

    ICRSCoordSys::ICRSCoordSys(double date)
    :
        MeanCoordSys("icrs", date)
    {};

    boost::shared_ptr<CoordSys> ICRSCoordSys::clone() const {
        return clone(getDate());
    }

    boost::shared_ptr<CoordSys> ICRSCoordSys::clone(double date) const {
        return boost::shared_ptr<CoordSys>(new ICRSCoordSys(date));
    };

    Coord ICRSCoordSys::fromFK5J2000(Coord const &coord, Site const &site) const {
        // use the excellent approximation that ICRS = ICRS J2000
        double const fromDate = 2000.0;
        double const toDate = this->_date;
        
        Eigen::Vector3d fromPos = coord.getVecPos();
        Eigen::Vector3d fromVel = coord.getVecVel();
        
        // correct for velocity (proper motion and radial velocity)
        Eigen::Vector3d toPos = fromPos + (fromVel * (toDate - fromDate));
        
        return Coord(toPos, fromVel);
    };
    
    Coord ICRSCoordSys::toFK5J2000(Coord const &coord, Site const &site) const {
        // use the excellent approximation that ICRS = ICRS J2000
        double const fromDate = this->_date;
        double const toDate = 2000.0;
        
        Eigen::Vector3d fromPos = coord.getVecPos();
        Eigen::Vector3d fromVel = coord.getVecVel();
        
        // correct for velocity (proper motion and radial velocity)
        Eigen::Vector3d toPos = fromPos + (fromVel * (toDate - fromDate));
        
        return Coord(toPos, fromVel);
    };
    
}
