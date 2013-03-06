#include <sstream>
#include <stdexcept>
#include "boost/make_shared.hpp"
#include "coordConv/coordSys.h"

namespace {
    const double DeltaT = 0.01;
}

namespace coordConv {
    
    Coord CoordSys::convertFrom(CoordSys const &fromCoordSys, Coord const &fromCoord, Site const &site) const {
        Coord icrsCoord = fromCoordSys.toFK5J2000(fromCoord, site);
        return fromFK5J2000(icrsCoord, site);
    }

    PVTCoord CoordSys::convertFrom(CoordSys const &fromCoordSys, PVTCoord const &fromPVTCoord, Site const &site) const {
        double const tai0 = fromPVTCoord.getTAI();
        double const tai1 = tai0 + DeltaT;
        Coord coord0 = this->convertFrom(fromCoordSys, fromPVTCoord.getCoord(tai0), site);
        Coord coord1 = this->convertFrom(fromCoordSys, fromPVTCoord.getCoord(tai1), site);
        return PVTCoord(coord0, coord1, tai0, DeltaT);
    }

    Coord CoordSys::convertFrom(double &toDir, double &scaleChange, CoordSys const &fromCoordSys, Coord const &fromCoord, double fromDir, Site const &site) const {
        Coord toCoord = convertFrom(fromCoordSys, fromCoord, site);

        double const FromDist = 1e-3;
        double dumDir;
        Coord offFromCoord = fromCoord.offset(dumDir, fromDir, FromDist);
        Coord offToCoord = convertFrom(fromCoordSys, offFromCoord, site);
        scaleChange = toCoord.angularSeparation(offToCoord) / FromDist;
        toDir = toCoord.angleTo(offToCoord);
        return toCoord;
    }

    PVTCoord CoordSys::convertFrom(PVT &toDir, double &scaleChange, CoordSys const &fromCoordSys, PVTCoord const &fromPVTCoord, PVT const &fromDir, Site const &site) const {
        double const tai0 = fromPVTCoord.getTAI();
        double const tai1 = tai0 + DeltaT;
        double toDirPair[2], dumScaleCh;
        Coord coord0 = this->convertFrom(toDirPair[0], scaleChange, fromCoordSys, fromPVTCoord.getCoord(tai0), fromDir.getPos(tai0), site);
        Coord coord1 = this->convertFrom(toDirPair[1], dumScaleCh,  fromCoordSys, fromPVTCoord.getCoord(tai1), fromDir.getPos(tai1), site);
        toDir.setFromAnglePair(toDirPair, tai0, DeltaT);
        return PVTCoord(coord0, coord1, tai0, DeltaT);
    }
    
    std::string CoordSys::asString() const {
        std::ostringstream os;
        os << "CoordSys(" << getName() << ", " << getDate() << ")";
        return os.str();
    }

    boost::shared_ptr<coordConv::CoordSys> makeCoordSys(std::string const &name, double date) {
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
        } else {
            std::ostringstream os;
            os << "Unknown coordinate system name: " << name;
            throw std::runtime_error(os.str());
        }
    }

}
