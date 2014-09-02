#include <stdexcept>
#include <sstream>
#include "coordConv/time.h"
#include "coordConv/appGeoFromAppTopo.h"
#include "coordConv/appTopoFromAppGeo.h"
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
            _appGeoCoordSys.setCurrDate(_appGeoCoordSys.dateFromTAI(date));
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
        return appTopoFromAppGeo(appGeoCoord, site, _date);
    };

    Coord AppTopoCoordSys::toFK5J2000(Coord const &coord, Site const &site) const {
        Coord appGeoCoord = appGeoFromAppTopo(coord, site, _date);
        return _appGeoCoordSys.toFK5J2000(appGeoCoord, site);
    };

    std::string AppTopoCoordSys::__repr__() const {
        std::ostringstream os;
        os << "AppTopoCoordSys(" << getDate() << ")";
        return os.str();
    }

}
