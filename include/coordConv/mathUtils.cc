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
        if (wrappedAng < -180) {
            wrappedAng += 360;
        } else if (wrappedAng > 180) {
            wrappedAng -= 360;
        }
        if (wrappedAng == 180) {
            wrappedAng = -180;
        };
        return wrappedAng;
    }

    inline double wrapNear(double ang, double nearAng) {
        double wrappedAng = nearAng + wrapCtr(ang - nearAng);
        double const delta = wrappedAng - nearAng;
        // roundoff error can cause out-of-range values; the following fixes those
        // and even seems to preserve ang - nearAng < 180, though I'm not sure why
        if (delta < -180) {
            wrappedAng += 360;
        } else if (delta >= 180) {
            wrappedAng -= 360;
        }
        return wrappedAng;
    }
}
