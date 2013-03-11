#include "coordConv/pvt.h"

namespace {
    const double DeltaT = 0.01;
}

namespace coordConv {

    bool polarFromXY(PVT &r, PVT &theta, PVT const &x, PVT const &y, double tai) {
        double rArr[2], thetaArr[2];
        bool atPole;
        for (int i = 0; i < 2; ++i) {
            double tempTAI = tai + (i * DeltaT);
            atPole |= polarFromXY(rArr[i], thetaArr[i], x.getPos(tai), y.getPos(tai));
        }
        r.setFromAnglePair(rArr, tai, DeltaT);
        theta.setFromAnglePair(thetaArr, tai, DeltaT);
        return atPole;
    }

    void xyFromPolar(PVT &x, PVT &y, PVT const &r, PVT const &theta, double tai) {
        double xArr[2], yArr[2];
        for (int i = 0; i < 2; ++i) {
            double tempTAI = tai + (i * DeltaT);
            xyFromPolar(xArr[i], yArr[i], r.getPos(tai), theta.getPos(tai));
        }
        x.setFromAnglePair(xArr, tai, DeltaT);
        y.setFromAnglePair(yArr, tai, DeltaT);
    }

}
