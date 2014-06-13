#include <stdexcept>
#include <sstream>
#include "coordConv/mathUtils.h"
#include "coordConv/coordSys.h"

namespace coordConv {

    OtherCoordSys::OtherCoordSys(std::string const &name, double date, DateTypeEnum dateType, bool isMean)
    :
        CoordSys(name, date, dateType, isMean, false)
    {};
    
    CoordSys::Ptr OtherCoordSys::clone() const {
        return clone(getDate());
    }

    CoordSys::Ptr OtherCoordSys::clone(double date) const {
        return CoordSys::Ptr(new OtherCoordSys(_name, date, _dateType, _isMean));
    };

    double OtherCoordSys::dateFromTAI(double tai) const {
        return tai;
    }

    Coord OtherCoordSys::fromFK5J2000(Coord const &coord, Site const &site) const {
        return Coord(); // null Coord
    }

    Coord OtherCoordSys::toFK5J2000(Coord const &coord, Site const &site) const {
        return Coord(); // null Coord
    }

    std::string OtherCoordSys::__repr__() const {
        std::ostringstream os;
        os << "OtherCoordSys(name=" << getName() << ", isMean=" << isMean() << ", date=" << getDate() << ")";
        return os.str();
    }

}

