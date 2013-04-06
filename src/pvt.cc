#include "coordConv/pvt.h"

namespace {
    const double DeltaT = 0.01;
}

namespace coordConv {

    bool polarFromXY(PVT &r, PVT &theta, PVT const &x, PVT const &y, double tai) {
        double rArr[2], thetaArr[2];
        bool atPole = false;
        for (int i = 0; i < 2; ++i) {
            double tempTAI = tai + (i * DeltaT);
            atPole |= polarFromXY(rArr[i], thetaArr[i], x.getPos(tempTAI), y.getPos(tempTAI));
        }
        r.setFromPair(rArr, tai, DeltaT, false);
        theta.setFromPair(thetaArr, tai, DeltaT, true);
        return atPole;
    }

    void xyFromPolar(PVT &x, PVT &y, PVT const &r, PVT const &theta, double tai) {
        double xArr[2], yArr[2];
        for (int i = 0; i < 2; ++i) {
            double tempTAI = tai + (i * DeltaT);
            xyFromPolar(xArr[i], yArr[i], r.getPos(tempTAI), theta.getPos(tempTAI));
        }
        x.setFromPair(xArr, tai, DeltaT, false);
        y.setFromPair(yArr, tai, DeltaT, false);
    }

}
