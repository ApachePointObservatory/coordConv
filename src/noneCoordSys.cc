#include <stdexcept>
#include <sstream>
#include "coordConv/mathUtils.h"
#include "coordConv/coordSys.h"

namespace coordConv {

    NoneCoordSys::NoneCoordSys(double date)
    :
        ApparentCoordSys("none", date)
    {};
    
    CoordSys::Ptr NoneCoordSys::clone() const {
        return clone(getDate());
    }

    CoordSys::Ptr NoneCoordSys::clone(double date) const {
        return CoordSys::Ptr(new NoneCoordSys(date));
    };

    Coord NoneCoordSys::fromFK5J2000(Coord const &coord, Site const &site) const {
        return Coord(); // null Coord
    }

    Coord NoneCoordSys::toFK5J2000(Coord const &coord, Site const &site) const {
        return Coord(); // null Coord
    }

    std::string NoneCoordSys::__repr__() const {
        return std::string("NoneCoordSys()");
    }

}

