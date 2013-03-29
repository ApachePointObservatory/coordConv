#pragma once

#include <cmath>
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
        double wrappedAng = std::fmod(ang, 360.0);
        if (wrappedAng < 0.0) {
            wrappedAng += 360.0;
        }
        if (wrappedAng == 360.0) {
            // this can happen if wrappedAng is so small that wrappedAng + 360 rounds to 360
            wrappedAng = 0.0;
        };
        return wrappedAng;
    }

    inline double wrapCtr(double ang) {
        return wrapPos(ang + 180.0) - 180.0;
    }

    inline double wrapNear(double ang, double nearAng) {
        double wrappedAng = nearAng + wrapCtr(ang - nearAng);
        if ((wrappedAng - nearAng) < -180) {
            // I'm not sure exactly how this happens, but unit tests showed it can
            wrappedAng += 360;
        }
        return wrappedAng;
    }
}
