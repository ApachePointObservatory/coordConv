#include "coordConv/coordSys.h"

namespace coordConv {

    ApparentCoordSys::ApparentCoordSys(std::string const &name, double date)
    :
        CoordSys(name, date, false, true)
    {}
    
    // works for apparent topocentric and observed; override for geocentric
    double ApparentCoordSys::dateFromTAI(double tai) const {
        return tai;
    }

}
