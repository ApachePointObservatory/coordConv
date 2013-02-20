#include "coordConv/mathUtils.h"

namespace coordConv {

    void rotXY(Eigen::Vector3d &toVec, Eigen::Vector3d const &fromVec, double xAng, double yAng) {

        double const sinX = sind(xAng);
        double const cosX = cosd(xAng);
        double const sinY = sind(yAng);
        double const cosY = cosd(yAng);

        Eigen::Matrix3d rotMat;
        rotMat <<
              cosY,  sinX * sinY,    cosX * sinY,
               0.0,  cosX,         - sinX,
            - sinY,  sinX * cosY,    cosX * cosY;

        toVec = rotMat * fromVec;
    }
    
}
