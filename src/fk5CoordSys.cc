#include <stdexcept>
#include <sstream>
#include "slalib.h"
#include "coordConv/mathUtils.h"
#include "coordConv/coordSys.h"

namespace coordConv {

    FK5CoordSys::FK5CoordSys(double date)
    :
        MeanCoordSys("fk5", date),
        _to2000PrecMat()
    {
        setDate(date);
    };
    
    void FK5CoordSys::setDate(double date) {
        this->_date = date;
        if (isFinite(date)) {
            double precMatCArr[3][3];
            slaPrec(date, 2000.0, precMatCArr);
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    _to2000PrecMat(i,j) = precMatCArr[i][j];
                }
            }
        }
    }

    CoordSys::Ptr FK5CoordSys::clone() const {
        return clone(getDate());
    }

    CoordSys::Ptr FK5CoordSys::clone(double date) const {
        return CoordSys::Ptr(new FK5CoordSys(date));
    };

    Coord FK5CoordSys::fromFK5J2000(Coord const &coord, Site const &site) const {
        // use the excellent approximation that ICRS = FK5 J2000
        double const fromDate = 2000.0;
        double const toDate = this->_date;
        
        Eigen::Vector3d fk5J2000Pos = coord.getVecPos();
        Eigen::Vector3d fk5J2000PM = coord.getVecPM();
        
        // correct for velocity (proper motion and radial velocity)
        Eigen::Vector3d tempPos = fk5J2000Pos + (fk5J2000PM * (toDate - fromDate));

        // precess position and velocity
        Eigen::Vector3d fk5Pos, fk5Vel;
        fk5Pos = _to2000PrecMat.transpose() * tempPos;
        fk5Vel = _to2000PrecMat.transpose() * fk5J2000PM;
        
        return Coord(fk5Pos, fk5Vel);
    };
     
    Coord FK5CoordSys::toFK5J2000(Coord const &coord, Site const &site) const {
        // use the excellent approximation that ICRS = FK5 J2000
        double const fromDate = this->_date;
        double const toDate = 2000.0;
        
        Eigen::Vector3d fk5Pos = coord.getVecPos();
        Eigen::Vector3d fk5Vel = coord.getVecPM();
        
        // correct for velocity (proper motion and radial velocity)
        Eigen::Vector3d tempPos = fk5Pos + (fk5Vel * (toDate - fromDate));

        // precess position and velocity
        Eigen::Vector3d fk5J2000Pos, fk5J2000PM;
        fk5J2000Pos = _to2000PrecMat * tempPos;
        fk5J2000PM = _to2000PrecMat * fk5Vel;
        
        return Coord(fk5J2000Pos, fk5J2000PM);
    };

}
