#include <stdexcept>
#include <sstream>
#include "coordConv/mathUtils.h"
#include "coordConv/coordSys.h"

namespace coordConv {

    NoneCoordSys::NoneCoordSys(double date)
    :
        OtherCoordSys("none", date, false)
    {};
    
    CoordSys::Ptr NoneCoordSys::clone() const {
        return clone(getDate());
    }

    CoordSys::Ptr NoneCoordSys::clone(double date) const {
        return CoordSys::Ptr(new NoneCoordSys(date));
    };

    std::string NoneCoordSys::__repr__() const {
        return std::string("NoneCoordSys()");
    }

}

