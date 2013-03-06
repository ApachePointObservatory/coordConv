#include "coordConv/coordSys.h"

namespace {

    // WARNING: this data is transposed from how it will be in the Eigen matrix
    static double const _fromFK5J2000CArr[9] = {
        -5.487553969571600e-02, +4.941094533056070e-01, -8.676661358478490e-01,
        -8.734371079953150e-01, -4.448295894318790e-01, -1.980763861308200e-01,
        -4.838349858369940e-01, +7.469822518105100e-01, +4.559837957210930e-01
    };
    static const Eigen::Map<const Eigen::Matrix3d> _fromFK5J2000RotMat(_fromFK5J2000CArr);

}

namespace coordConv {

    GalCoordSys::GalCoordSys(double date)
    :
        MeanCoordSys("gal", date)
    {}

    boost::shared_ptr<CoordSys> GalCoordSys::clone() const {
        return clone(getDate());
    }

    boost::shared_ptr<CoordSys> GalCoordSys::clone(double date) const {
        return boost::shared_ptr<CoordSys>(new GalCoordSys(date));
    };

    Coord GalCoordSys::fromFK5J2000(Coord const &coord, Site const &site) const {
        double const fromDate = 2000.0;
        double const toDate = this->_date;
        
        // adjust space velocity
        Eigen::Vector3d fk5J2000Pos = coord.getVecPos();
        Eigen::Vector3d fk5J2000Vel = coord.getVecVel();
        Eigen::Vector3d adjIcrsPos = fk5J2000Pos + (fk5J2000Vel * (toDate - fromDate));
    
        // convert position to galactic coordinates
        Eigen::Vector3d galPos, galVel;
        galPos =  _fromFK5J2000RotMat * adjIcrsPos;
        galVel =  _fromFK5J2000RotMat * fk5J2000Vel;
        
        return Coord(galPos, galVel);
    }
    
    Coord GalCoordSys::toFK5J2000(Coord const &coord, Site const &site) const {
        double const fromDate = this->_date;
        double const toDate = 2000.0;

        // correct for velocity (proper motion and radial velocity)
        Eigen::Vector3d galPos = coord.getVecPos();
        Eigen::Vector3d galVel = coord.getVecVel();
        Eigen::Vector3d adjGalPos = galPos + (galVel * (toDate - fromDate));

        // convert position to J2000 coordinates (by INVERSE rotation)
        Eigen::Vector3d fk5J2000Pos, fk5J2000Vel;
        fk5J2000Pos = _fromFK5J2000RotMat.transpose() * adjGalPos;
        fk5J2000Vel = _fromFK5J2000RotMat.transpose() * galVel;
        
        return Coord(fk5J2000Pos, fk5J2000Vel);
    }

}
