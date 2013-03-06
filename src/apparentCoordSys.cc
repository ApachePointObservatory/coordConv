#include "coordConv/coordSys.h"

namespace coordConv {

    ApparentCoordSys::ApparentCoordSys(std::string const &name, double date)
    :
        CoordSys(name, false, date)
    {}
    
    // works for apparent topocentric and observed; override for geocentric
    double ApparentCoordSys::dateFromTAI(double tai) const {
        return tai;
    }

    // warning: simply zeros the velocity component (uninteresting)
    Coord ApparentCoordSys::removePM(Coord const &coord, double tai) const {
        Coord zpmCoord = Coord(coord.getVecPos());
        return zpmCoord;
    }

}
