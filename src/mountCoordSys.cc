#include <stdexcept>
#include <sstream>
#include "coordConv/mathUtils.h"
#include "coordConv/coordSys.h"

namespace coordConv {

    MountCoordSys::MountCoordSys(double date)
    :
        ApparentCoordSys("mount", date)
    {};
    
    CoordSys::Ptr MountCoordSys::clone() const {
        return clone(getDate());
    }

    CoordSys::Ptr MountCoordSys::clone(double date) const {
        return CoordSys::Ptr(new MountCoordSys(date));
    };

    Coord MountCoordSys::fromFK5J2000(Coord const &coord, Site const &site) const {
        return Coord(); // null Coord
    }

    Coord MountCoordSys::toFK5J2000(Coord const &coord, Site const &site) const {
        return Coord(); // null Coord
    }

}

