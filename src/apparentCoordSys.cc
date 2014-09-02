#include "coordConv/coordSys.h"

namespace coordConv {

    ApparentCoordSys::ApparentCoordSys(std::string const &name, double date, DateTypeEnum dateType)
    :
        CoordSys(name, date, dateType, false, true)
    {}

}
