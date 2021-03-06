#include <stdexcept>
#include <sstream>
#include "coordConv/mathUtils.h"
#include "coordConv/coordSys.h"

namespace coordConv {

    ICRSCoordSys::ICRSCoordSys(double date)
    :
        MeanCoordSys("icrs", date)
    {};

    // note: I tried only defining clone() in CoordSys but SWIG would not recognize clone()
    // in any inherited coordinate system.
    CoordSys::Ptr ICRSCoordSys::clone() const {
        return clone(getDate());
    }

    CoordSys::Ptr ICRSCoordSys::clone(double date) const {
        return CoordSys::Ptr(new ICRSCoordSys(date));
    };

    Coord ICRSCoordSys::fromFK5J2000(Coord const &coord, Site const &site) const {
        // use the excellent approximation that ICRS = ICRS J2000
        double const fromDate = 2000.0;
        double const toDate = this->_date;
        
        Eigen::Vector3d fromPos = coord.getVecPos();
        Eigen::Vector3d fromPM = coord.getVecPM();
        
        // correct for velocity (proper motion and radial velocity)
        Eigen::Vector3d toPos = fromPos + (fromPM * (toDate - fromDate));
        
        return Coord(toPos, fromPM);
    };
    
    Coord ICRSCoordSys::toFK5J2000(Coord const &coord, Site const &site) const {
        // use the excellent approximation that ICRS = ICRS J2000
        double const fromDate = this->_date;
        double const toDate = 2000.0;
        
        Eigen::Vector3d fromPos = coord.getVecPos();
        Eigen::Vector3d fromPM = coord.getVecPM();
        
        // correct for velocity (proper motion and radial velocity)
        Eigen::Vector3d toPos = fromPos + (fromPM * (toDate - fromDate));
        
        return Coord(toPos, fromPM);
    };

    std::string ICRSCoordSys::__repr__() const {
        std::ostringstream os;
        os << "ICRSCoordSys(" << getDate() << ")";
        return os.str();
    }
    
}
