#pragma once

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>
/*
Define inline math utilities
*/
namespace coordConv {
    
    inline void rot2D(double &rotX, double &rotY, double x, double y, double ang) {
        double sinAng = sind(ang);
        double cosAng = cosd(ang);

        rotX = cosAng * x - sinAng * y;
        rotY = sinAng * x + cosAng * y;
    }

    inline double wrapPos(double ang) {
        // put angle into range (-360, 360), then finish the job
        double wrappedAng = std::fmod(ang, 360);
        if (wrappedAng < 0) {
            wrappedAng += 360;
        }
        if (wrappedAng == 360) {
            // this can happen if wrappedAng is so small that wrappedAng + 360 rounds to 360
            wrappedAng = 0;
        };
        return wrappedAng;
    }

    inline double wrapCtr(double ang) {
        // put angle into range (-360, 360), then finish the job
        double wrappedAng = std::fmod(ang, 360);
        if (wrappedAng >= 180) {
            wrappedAng -= 360;
            if (wrappedAng < -180) {
                // handle roundoff error
                wrappedAng = -180;
            }
        } else if (wrappedAng < -180) {
            wrappedAng += 360;
            if (wrappedAng >= 180) {
                // handle roundoff error
                wrappedAng = -180;
            }
        }
        return wrappedAng;
    }

    inline double wrapNear(double ang, double refAng) {
        double wrappedAng = refAng + wrapCtr(ang - refAng);

        // roundoff error can cause slightly out-of-range values; the following fixes those
        // (and even seems to preserve ang - refAng < 180, though I'm not sure why)
        if (wrappedAng - refAng >= 180) {
            wrappedAng -= 360;
        }
        // avoid if-else in case wrappedAng -= 360 results in wrappedAng - refAng slightly less than -180;
        // maximum relative roundoff error for addition is 2 epsilon
        if (wrappedAng - refAng < -180) {
            wrappedAng -= wrappedAng * 2.0 * DoubleEpsilon;
        }
        return wrappedAng;
    }
}
