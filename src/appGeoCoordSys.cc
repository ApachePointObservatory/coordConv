#include <stdexcept>
#include <sstream>
#include "slalib.h"
#include "coordConv/mathUtils.h"
#include "coordConv/physConst.h"
#include "coordConv/coordSys.h"

namespace coordConv {

    AppGeoCoordSys::AppGeoCoordSys(double date, double maxAge)
    :
        CoordSys("appgeo", date),
        _maxAge(maxAge),
        _cachedDate(DoubleNaN)
    {
        setDate(date);
    };

    boost::shared_ptr<CoordSys> AppGeoCoordSys::clone() const {
        return clone(getDate());
    }

    boost::shared_ptr<CoordSys> AppGeoCoordSys::clone(double date) const {
        return boost::shared_ptr<CoordSys>(new AppGeoCoordSys(date));
    };
    
    void AppGeoCoordSys::setDate(double date) {
        this->_date = date;
        if (isFinite(date)) {
            if (std::abs(date - _cachedDate) < _maxAge) {
                return;
            }
            double tdbDays = slaEpj2d(date);
            double amprms[21];
            slaMappa(2000.0, tdbDays, amprms);
        
            _pmSpan = amprms[0];
            _gravRad = amprms[7];
            _gammaI = amprms[11];
            for (int i = 0; i < 3; ++i) {
                _bcPos(i) = amprms[1+i];
                _hcDir(i) = amprms[4+i];
                _bcBeta(i) = amprms[8+i];
                for (int j = 0; j < 3; ++j) {
                    _pnMat(i,j) = amprms[12+(i*3)+j];
                }
            }
            _cachedDate = date;
        }
    }

    Coord AppGeoCoordSys::fromICRS(Coord const &coord, Site const &site) const {
        
        Eigen::Vector3d icrsPos = coord.getVecPos();
        Eigen::Vector3d icrsVel = coord.getVecVel();

        // correct for velocity and Earth's offset from the barycenter
        Eigen::Vector3d pos1 = icrsPos + (icrsVel * _pmSpan) - _bcPos;

        // here is where the correction for sun's gravity belongs
        Eigen::Vector3d pos2 = pos1;

        // correct for annual aberration
        double pos2Mag = pos2.norm();
        double dot2 = pos2.dot(_bcBeta) / pos2Mag;
        double vfac = pos2Mag * (1.0 + dot2 / (1.0 + _gammaI)); // the presence of pos2Mag is due to light travel time from the target
        Eigen::Vector3d pos3 = ((_gammaI * pos2) + (vfac * _bcBeta)) / (1.0 + dot2);

        // correct position for precession and nutation
        Eigen::Vector3d appGeoPos = _pnMat * pos3;
        return Coord(appGeoPos);
    };

    /**
    Perform the inverse transform of fromICRS.
    
    Unfortunately, some of the equations (e.g. annual aberration) are not invertable,
    so they have been solved by iteration. To make the code easier to follow, the symbols used here
    are identical to those used in fromICRS.

    The convergence criterion is set by the magic numbers MaxIter and Accuracy.

    References:
      cnv_J2AppGeo*
      ABERAT*, an APPLE (J2000) subroutine; U.S. Naval Observatory
      P.T. Wallace, slaMAPQK (a SLALIB subroutine); Starlink, RGO
      P.T. Wallace, "Proposals for Keck Tel. Point. Algorithms," 1986 (unpub.)
      "The Astronomical Almanac" for 1978, U.S. Naval Observatory
      *these use physical units instead of direction cosines
    */    
    Coord AppGeoCoordSys::toICRS(Coord const &coord, Site const &site) const {
        Eigen::Vector3d appGeoPos = coord.getVecPos();

        /// if the number of iterations exceeds "MaxIter" before converging, raise an exception
        const int MaxIter = 20;
        /// if all three components of (P this iter - P last iter) / |P|
        /// are less than "Accuracy", then the iteration has converged.
        const double Accuracy = 1.0e-10;
        
        // compute constants needed to check iteration
        const double approxMagP = appGeoPos.norm();
        const double allowedErr = Accuracy * approxMagP;

        // correct position for nutation and precession
        Eigen::Vector3d pos3 = _pnMat.transpose() * appGeoPos;

        // iterate to correct for annual aberration
        int iter = 0;
        double maxErr = approxMagP;
        Eigen::Vector3d pos2 = pos3;
        while (maxErr > allowedErr) {
            iter += 1;
            if (iter > MaxIter) {
                std::ostringstream os;
                os << "aberration correction failed to converge in " << MaxIter <<
                    " iterations; error = " << maxErr << " > " << allowedErr << " allowed";
                throw std::runtime_error(os.str());
            }
            
            double p2Mag = pos2.norm();
            double dot2 = pos2.dot(_bcBeta) / p2Mag;
            double fac = p2Mag * (1.0 + (dot2 / (1.0 + _gammaI)));
            Eigen::Vector3d oldP2 = pos2; // ??? is a copy required?
            pos2 = (((1.0 + dot2) * pos3) - (fac * _bcBeta)) / _gammaI;
            maxErr = (pos2 - oldP2).array().abs().maxCoeff();
        }

        // here is where the (iterative) correction for sun's gravity belongs
        Eigen::Vector3d pos1 = pos2;

        // correct for Earth's offset from the barycenter
        Eigen::Vector3d icrsPos = pos1 + _bcPos;
        
        return Coord(icrsPos);
    }

}
