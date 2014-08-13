#include <iomanip>
#include <sstream>
#include "coordConv/pvt.h"

namespace {
    const double DeltaT = 0.01;
}

namespace coordConv {

    std::string PVT::__repr__() const {
        std::ostringstream os;
        os << *this;
        return os.str();
    }

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

    void rot2D(PVT &rotX, PVT &rotY, PVT const &x, PVT const &y, double ang, double tai) {
        double rotXArr[2], rotYArr[2];
        for (int i = 0; i < 2; ++i) {
            double tempTAI = tai + (i * DeltaT);
            rot2D(rotXArr[i], rotYArr[i], x.getPos(tempTAI), y.getPos(tempTAI), ang);
        }
        rotX.setFromPair(rotXArr, tai, DeltaT, false);
        rotY.setFromPair(rotYArr, tai, DeltaT, false);
    }

    std::ostream &operator<<(std::ostream &os, PVT const &pvt) {
        std::ios_base::fmtflags oldFlags = os.flags();
        std::streamsize const oldPrecision = os.precision();
        os << std::fixed
            << "PVT(" << std::setprecision(7) << pvt.pos << ", " << pvt.vel << ", "
            << std::setprecision(6) << pvt.t
            << ")" << std::setprecision(oldPrecision);
        os.flags(oldFlags);
        return os;
    }

}
