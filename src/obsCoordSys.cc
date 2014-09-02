#include <stdexcept>
#include <sstream>
#include "coordConv/mathUtils.h"
#include "coordConv/appTopoFromObs.h"
#include "coordConv/obsFromAppTopo.h"
#include "coordConv/coordSys.h"

namespace coordConv {

    ObsCoordSys::ObsCoordSys(double date)
    :
        ApparentCoordSys("obs", date),
        _appTopoCoordSys()
    {
        setDate(date);
    };
    
    void ObsCoordSys::_setDate(double date) const {
        if (date > 0) {
            _appTopoCoordSys.setCurrDate(date);
        }
        this->_date = date;
    }
    
    CoordSys::Ptr ObsCoordSys::clone() const {
        return clone(getDate());
    }

    CoordSys::Ptr ObsCoordSys::clone(double date) const {
        return CoordSys::Ptr(new ObsCoordSys(date));
    };

    Coord ObsCoordSys::fromFK5J2000(Coord const &coord, Site const &site) const {
        Coord appTopoCoord = _appTopoCoordSys.fromFK5J2000(coord, site);
        return obsFromAppTopo(appTopoCoord, site);
    }

    Coord ObsCoordSys::toFK5J2000(Coord const &coord, Site const &site) const {
        Coord appTopoCoord = appTopoFromObs(coord, site);
        return _appTopoCoordSys.toFK5J2000(appTopoCoord, site);
    }

    std::string ObsCoordSys::__repr__() const {
        std::ostringstream os;
        os << "ObsCoordSys(" << getDate() << ")";
        return os.str();
    }

}

