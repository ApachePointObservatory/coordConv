#include <stdexcept>
#include "coordConv/angSideAng.h"
#include "coordConv/offsetSph.h"

namespace coordConv {

    void offsetSph(
        double &destEquatAng,
        double &destPolarAng,
        double &destOrient,
        double srcEquatAng,
        double srcPolarAng,
        double srcOrient,
        double dist
    ) {
    
        double sideA = 90.0 - srcPolarAng;
        double angB = 90.0 - srcOrient;
        double sideC = dist;
        double angA, sideB, angC;
        if (angSideAng(angA, sideB, angC, sideA, angB, sideC)) {
            throw std::runtime_error("Cannot compute");
        }
        destEquatAng = srcEquatAng + angC;
        destPolarAng = 90.0 - sideB;
        destOrient = angA - 90.0;
    }

}
