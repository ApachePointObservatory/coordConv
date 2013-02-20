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

}
