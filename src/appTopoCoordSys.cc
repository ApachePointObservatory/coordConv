#include <stdexcept>
#include <sstream>
#include "coordConv/time.h"
#include "coordConv/azAltFromHADec.h"
#include "coordConv/haDecFromAzAlt.h"
#include "coordConv/mathUtils.h"
#include "coordConv/coordSys.h"

namespace coordConv {

    AppTopoCoordSys::AppTopoCoordSys(double date)
    :
        ApparentCoordSys("apptopo", date)
    {
        setDate(date);
    };
    
    void AppTopoCoordSys::_setDate(double date) const {
        if (date > 0) {
            _appGeoCoordSys.setCurrTAI(date);
        }
        this->_date = date;
    }
    
    CoordSys::Ptr AppTopoCoordSys::clone() const {
        return clone(getDate());
    }

    CoordSys::Ptr AppTopoCoordSys::clone(double date) const {
        return CoordSys::Ptr(new AppTopoCoordSys(date));
    };

    Coord AppTopoCoordSys::fromFK5J2000(Coord const &coord, Site const &site) const {
        Coord appGeoCoord = _appGeoCoordSys.fromFK5J2000(coord, site);
        return fromAppGeo(appGeoCoord, site);
    };

    Coord AppTopoCoordSys::toFK5J2000(Coord const &coord, Site const &site) const {
        Coord appGeoCoord = toAppGeo(coord, site);
        return _appGeoCoordSys.toFK5J2000(appGeoCoord, site);
    };

    Coord AppTopoCoordSys::fromAppGeo(Coord const &coord, Site const &site) const {
        Eigen::Vector3d appGeoPos = coord.getVecPos();

        double const last = lastFromTAI(this->_date, site);
        double const sinLAST = sind(last);
        double const cosLAST = cosd(last);

        // rotate position and offset from RA/Dec to (-HA)/Dec
        Eigen::Vector3d posA;
        posA <<
            + (cosLAST * appGeoPos(0)) + (sinLAST * appGeoPos(1)),
            - (sinLAST * appGeoPos(0)) + (cosLAST * appGeoPos(1)),
               appGeoPos(2);

        // correct position for diurnal parallax
        Eigen::Vector3d posB = posA - site.pos;

        // correct position for diurnal aberration; follows Pat Wallace's slaAOPQK
        Eigen::Vector3d posC;
        double bMag = posB.norm();
        double diurAbScaleCorr = 1.0 - (site.diurAbMag * (posB(1) / bMag));
        posC(0) =  posB(0) * diurAbScaleCorr;
        posC(1) = (posB(1) + (site.diurAbMag * bMag)) * diurAbScaleCorr;
        posC(2) =  posB(2) * diurAbScaleCorr;

        // rotate position from -HA/Dec to alt/az; use latitude corrected for pole wander
        Eigen::Vector3d appTopoPos;
        azAltFromHADec(appTopoPos, posC, site.corrLat);
        
        return Coord(appTopoPos);
    }

    Coord AppTopoCoordSys::toAppGeo(Coord const &coord, Site const &site) const {
        // note: the variable names are identical to appTopoFromAppGeo.cc;
        // hence the reversed order of intermediate position names and the use of bMag instead of cMag

        double const last = lastFromTAI(this->_date, site);
        double const sinLAST = sind(last);
        double const cosLAST = cosd(last);
        
        Eigen::Vector3d appTopoPos = coord.getVecPos();

        // rotate position from alt/az to -HA/Dec; use latitude corrected for pole wander
        Eigen::Vector3d posC;
        haDecFromAzAlt(posC, appTopoPos, site.corrLat);

        // remove correction for diurnal aberration
        // following Pat Wallace's slaOAPQK, the same equation is used
        // as applying the correction, but the sign of diurAbMag is reversed
        Eigen::Vector3d posB;
        double bMag = posC.norm();
        double diurAbScaleCorr = 1.0 + (site.diurAbMag * (posC(1) / bMag));
        posB(0) =  posC(0) * diurAbScaleCorr;
        posB(1) = (posC(1) - (site.diurAbMag * bMag)) * diurAbScaleCorr;
        posB(2) =  posC(2) * diurAbScaleCorr;

        // correct position for diurnal parallax (needed for planets, not stars)
        Eigen::Vector3d posA = posB + site.pos;

        // rotate position from -HA/Dec to RA/Dec
        Eigen::Vector3d appGeoPos;
        appGeoPos(0) = + (cosLAST * posA(0)) - (sinLAST * posA(1));
        appGeoPos(1) = + (sinLAST * posA(0)) + (cosLAST * posA(1));
        appGeoPos(2) =    posA(2);
        return Coord(appGeoPos);
    }

    std::string AppTopoCoordSys::__repr__() const {
        std::ostringstream os;
        os << "AppTopoCoordSys(" << getDate() << ")";
        return os.str();
    }

}
