#include <algorithm>
#include <cmath>
#include <limits>
#include "coordConv/mathUtils.h"

namespace coordConv {

    double hypot(double x, double y) {
        // assume that Eigen does it as well as anybody
        Eigen::Vector2d vec(x, y);
        return vec.norm();
    }

    bool polarFromXY (double &r, double &theta, double x, double y) {
        r = hypot(x, y);

        bool atOrig;
        if (r < std::numeric_limits<double>::epsilon()) {
           atOrig = true;
           theta = 0.0;
        } else {
           atOrig = false;
           theta = atan2d(y, x);
        }
        return atOrig;
    }

    void xyFromPolar(double &x, double &y, double r, double theta) {
        x = r * cosd(theta);
        y = r * sind(theta);
    }

    void computeRotationMatrix(Eigen::Matrix3d &rotMat, Eigen::Vector3d const &axis, double rotAngle) {
        // this code is a minor adaptation of LSST afw Coord::rotate
        double const c = cosd(rotAngle);
        double const mc = 1.0 - c;
        double const s = sind(rotAngle);
    
        double axisMag = axis.norm();
        double const ux = axis(0) / axisMag;
        double const uy = axis(1) / axisMag;
        double const uz = axis(2) / axisMag;
        
        rotMat <<
            (ux*ux + (1.0 - ux*ux)*c),  (ux*uy*mc - uz*s),          (ux*uz*mc + uy*s),
            (ux*uy*mc + uz*s),          (uy*uy + (1.0 - uy*uy)*c),  (uy*uz*mc - ux*s),
            (uz*ux*mc - uy*s),          (uy*uz*mc + ux*s),          (uz*uz + (1.0 - uz*uz)*c);
    }

}
