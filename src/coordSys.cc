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
    * @param[in] taiDate  TAI (MJD sec) date to use if CoordSys is apparent
    * @return a clone of coordSys; if coordinate system date type is TAI (apparent topocentric or observed)
        then set date of clone to taiDate; otherwise preserve the date.
    */
    inline coordConv::CoordSys::Ptr copyCoordSys(coordConv::CoordSys const &coordSys, double taiDate) {
        if (coordSys.getDateType() == coordConv::DateType_TAI) {
            return coordSys.clone(taiDate);
        } else {
            return coordSys.clone();
        }
    }
}

namespace coordConv {
    
    Coord CoordSys::convertFrom(CoordSys const &fromCoordSys, Coord const &fromCoord, Site const &site) const {
        Coord icrsCoord = fromCoordSys.toFK5J2000(fromCoord, site);
        return fromFK5J2000(icrsCoord, site);
    }

    PVTCoord CoordSys::convertFrom(CoordSys const &fromCoordSys, PVTCoord const &fromPVTCoord, Site const &site) const {
        double const tai0 = fromPVTCoord.getTAI();
        double const tai1 = tai0 + DeltaT;
        CoordSys::ConstPtr fromSys0Ptr = copyCoordSys(fromCoordSys, tai0);
        CoordSys::ConstPtr fromSys1Ptr = copyCoordSys(fromCoordSys, tai1);
        CoordSys::ConstPtr toSys0Ptr = copyCoordSys(*this, tai0);
        CoordSys::ConstPtr toSys1Ptr = copyCoordSys(*this, tai1);
        Coord toCoord0 = toSys0Ptr->convertFrom(*fromSys0Ptr, fromPVTCoord.getCoord(tai0), site);
        Coord toCoord1 = toSys1Ptr->convertFrom(*fromSys1Ptr, fromPVTCoord.getCoord(tai1), site);
        return PVTCoord(toCoord0, toCoord1, tai0, DeltaT);
    }

    Coord CoordSys::convertFrom(double &toDir, double &scaleChange, CoordSys const &fromCoordSys, Coord const &fromCoord, double fromDir, Site const &site) const {
        double const OffsetLength = 1e-3;
        Coord toCoord = convertFrom(fromCoordSys, fromCoord, site);
        double dumDir;
        Coord offFromCoord = fromCoord.offset(dumDir, fromDir, OffsetLength);
        Coord offToCoord = convertFrom(fromCoordSys, offFromCoord, site);
        scaleChange = toCoord.angularSeparation(offToCoord) / OffsetLength;
        toDir = toCoord.orientationTo(offToCoord);
        return toCoord;
    }

    PVTCoord CoordSys::convertFrom(PVT &toDir, double &scaleChange, CoordSys const &fromCoordSys, PVTCoord const &fromPVTCoord, PVT const &fromDir, Site const &site) const {
        double const tai0 = fromPVTCoord.getTAI();
        double const tai1 = tai0 + DeltaT;
        double toDirPair[2], dumScaleCh;
        CoordSys::ConstPtr fromSys0Ptr = copyCoordSys(fromCoordSys, tai0);
        CoordSys::ConstPtr fromSys1Ptr = copyCoordSys(fromCoordSys, tai1);
        CoordSys::ConstPtr toSys0Ptr = copyCoordSys(*this, tai0);
        CoordSys::ConstPtr toSys1Ptr = copyCoordSys(*this, tai1);
        Coord coord0 = toSys0Ptr->convertFrom(toDirPair[0], scaleChange, *fromSys0Ptr, fromPVTCoord.getCoord(tai0), fromDir.getPos(tai0), site);
        Coord coord1 = toSys1Ptr->convertFrom(toDirPair[1], dumScaleCh,  *fromSys1Ptr, fromPVTCoord.getCoord(tai1), fromDir.getPos(tai1), site);
        toDir.setFromPair(toDirPair, tai0, DeltaT, true);
        return PVTCoord(coord0, coord1, tai0, DeltaT);
    }

    PVTCoord CoordSys::removePM(PVTCoord const &pvtCoord) {
        std::vector<Coord> coordArr;
        double const tai = pvtCoord.getTAI();
        for (int i = 0; i < 2; ++i) {
            double tempTAI = tai + (i * DeltaT);
            coordArr.push_back(removePM(pvtCoord.getCoord(tempTAI), tempTAI));
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
