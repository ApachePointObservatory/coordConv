#include <sstream>
#include <stdexcept>
#include "boost/make_shared.hpp"
#include "coordConv/coordSys.h"

namespace coordConv {
    
    Coord CoordSys::convertFrom(CoordSys const &fromCoordSys, Coord const &fromCoord, Site const &site) const {
        Coord icrsCoord = fromCoordSys.toFK5J2000(fromCoord, site);
        return fromFK5J2000(icrsCoord, site);
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
