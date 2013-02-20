#include "coordConv/coordSys.h"

namespace {

    // WARNING: this data is transposed from how it will be in the Eigen matrix
    static double const _FromICRSCArr[9] = {
        -5.487553969571600e-02, +4.941094533056070e-01, -8.676661358478490e-01,
        -8.734371079953150e-01, -4.448295894318790e-01, -1.980763861308200e-01,
        -4.838349858369940e-01, +7.469822518105100e-01, +4.559837957210930e-01
    };
    static const Eigen::Map<const Eigen::Matrix3d> _FromICRSRotMat(_FromICRSCArr);

}

namespace coordConv {

    GalCoordSys::GalCoordSys(double date)
    :
        CoordSys("gal", date)
    {}

    boost::shared_ptr<CoordSys> GalCoordSys::clone() const {
        return clone(getDate());
    }

    boost::shared_ptr<CoordSys> GalCoordSys::clone(double date) const {
        return boost::shared_ptr<CoordSys>(new GalCoordSys(date));
    };

    Coord GalCoordSys::fromICRS(Coord const &coord, Site const &site) const {
        double const fromDate = 2000.0;
        double const toDate = this->_date;
        
        // adjust space velocity
        Eigen::Vector3d icrsPos = coord.getVecPos();
        Eigen::Vector3d icrsVel = coord.getVecVel();
        Eigen::Vector3d adjIcrsPos = icrsPos + (icrsVel * (toDate - fromDate));
    
        // convert position to galactic coordinates
        Eigen::Vector3d galPos, galVel;
        galPos =  _FromICRSRotMat * adjIcrsPos;
        galVel =  _FromICRSRotMat * icrsVel;
        
        return Coord(galPos, galVel);
    }
    
    Coord GalCoordSys::toICRS(Coord const &coord, Site const &site) const {
        double const fromDate = this->_date;
        double const toDate = 2000.0;

        // correct for velocity (proper motion and radial velocity)
        Eigen::Vector3d galPos = coord.getVecPos();
        Eigen::Vector3d galVel = coord.getVecVel();
        Eigen::Vector3d adjGalPos = galPos + (galVel * (toDate - fromDate));

        // convert position to J2000 coordinates (by INVERSE rotation)
        Eigen::Vector3d icrsPos, icrsVel;
        icrsPos = _FromICRSRotMat.transpose() * adjGalPos;
        icrsVel = _FromICRSRotMat.transpose() * galVel;
        
        return Coord(icrsPos, icrsVel);
    }

}
