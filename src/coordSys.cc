#include <sstream>
#include <stdexcept>
#include <vector>
#include "boost/make_shared.hpp"
#include "coordConv/coordSys.h"

namespace {
    const double DeltaT = 0.01;

    /**
    * Return a copy of a CoordSys, setting the date if coordSys is apparent, else leaving the date unchanged
    *
    * @param[in] coordSys  the CoordSys to copy
    * @param[in] taiDate  TAI (MJD sec) date of copy, if coordSys is current, else ignored
    * @return a clone of coordSys
    * @throw std::runtime_error if coordSys is current and taiDate is 0
    */
    inline coordConv::CoordSys::Ptr copyCoordSys(coordConv::CoordSys const &coordSys, double taiDate) {
        if (coordSys.isCurrent()) {
            if (taiDate == 0) {
                std::ostringstream os;
                os << "coordSys " << coordSys << " is current and taiDate=0";
                throw std::runtime_error(os.str());
            }
            return coordSys.clone(coordSys.dateFromTAI(taiDate));
        } else {
            return coordSys.clone();
        }
    }
}

namespace coordConv {
    
    Coord CoordSys::convertFrom(CoordSys const &fromCoordSys, Coord const &fromCoord, Site const &site, double tai) const {
        CoordSys::ConstPtr fromSysPtr = copyCoordSys(fromCoordSys, tai);
        CoordSys::ConstPtr toSysPtr = copyCoordSys(*this, tai);
        Coord icrsCoord = fromSysPtr->toFK5J2000(fromCoord, site);
        return toSysPtr->fromFK5J2000(icrsCoord, site);
    }

    PVTCoord CoordSys::convertFrom(CoordSys const &fromCoordSys, PVTCoord const &fromPVTCoord, Site const &site) const {
        std::vector<Coord> coordArr;
        double const tai = fromPVTCoord.getTAI();
        for (int i = 0; i < 2; ++i) {
            double evalTAI = tai + (i * DeltaT);
            coordArr.push_back(convertFrom(fromCoordSys, fromPVTCoord.getCoord(evalTAI), site, evalTAI));
        }
        return PVTCoord(coordArr[0], coordArr[1], tai, DeltaT);
    }

    Coord CoordSys::convertFrom(double &toDir, double &scaleChange, CoordSys const &fromCoordSys, Coord const &fromCoord, double fromDir, Site const &site, double tai) const {
        double const OffsetLength = 1e-3;
        Coord toCoord = convertFrom(fromCoordSys, fromCoord, site, tai);
        double dumDir;
        Coord offFromCoord = fromCoord.offset(dumDir, fromDir, OffsetLength);
        Coord offToCoord = convertFrom(fromCoordSys, offFromCoord, site, tai);
        scaleChange = toCoord.angularSeparation(offToCoord) / OffsetLength;
        toDir = toCoord.orientationTo(offToCoord);
        return toCoord;
    }

    PVTCoord CoordSys::convertFrom(PVT &toDir, double &scaleChange, CoordSys const &fromCoordSys, PVTCoord const &fromPVTCoord, PVT const &fromDir, Site const &site) const {
        std::vector<Coord> coordArr;
        double toDirPair[2], scaleChangePair[2];
        double const tai = fromPVTCoord.getTAI();
        for (int i = 0; i < 2; ++i) {
            double evalTAI = tai + (i * DeltaT);
            coordArr.push_back(convertFrom(toDirPair[i], scaleChangePair[i],
                fromCoordSys, fromPVTCoord.getCoord(evalTAI), fromDir.getPos(evalTAI), site, evalTAI));
        }
        scaleChange = scaleChangePair[0];
        toDir.setFromPair(toDirPair, tai, DeltaT, true);
        return PVTCoord(coordArr[0], coordArr[1], tai, DeltaT);
    }

    PVTCoord CoordSys::removePM(PVTCoord const &pvtCoord) {
        std::vector<Coord> coordArr;
        double const tai = pvtCoord.getTAI();
        for (int i = 0; i < 2; ++i) {
            double evalTAI = tai + (i * DeltaT);
            coordArr.push_back(removePM(pvtCoord.getCoord(evalTAI), evalTAI));
        }
        return PVTCoord(coordArr[0], coordArr[1], tai, DeltaT);
    }

    CoordSys::Ptr makeCoordSys(std::string const &name, double date) {
        if (name == "icrs") {
            return boost::make_shared<coordConv::ICRSCoordSys>(date);
        } else if (name == "fk5") {
            return boost::make_shared<coordConv::FK5CoordSys>(date);
        } else if (name == "fk4") {
            return boost::make_shared<coordConv::FK4CoordSys>(date);
        } else if (name == "gal") {
            return boost::make_shared<coordConv::GalCoordSys>(date);
        } else if (name == "appgeo") {
            return boost::make_shared<coordConv::AppGeoCoordSys>(date);
        } else if (name == "apptopo") {
            return boost::make_shared<coordConv::AppTopoCoordSys>(date);
        } else if (name == "obs") {
            return boost::make_shared<coordConv::ObsCoordSys>(date);
        } else if (name == "none") {
            return boost::make_shared<coordConv::NoneCoordSys>(date);
        } else {
            std::ostringstream os;
            os << "Unknown coordinate system name: " << name;
            throw std::invalid_argument(os.str());
        }
    }

    std::ostream &operator<<(std::ostream &os, CoordSys const &coordSys) {
        // use overloaded __repr__
        os << coordSys.__repr__();
        return os;
    }

}
