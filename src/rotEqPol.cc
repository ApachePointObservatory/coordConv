#include "coordConv/mathUtils.h"

namespace coordConv {

    void rotEqPol(Eigen::Vector3d &toVec, Eigen::Vector3d const &fromVec, double eqAng, double polarAng) {
        double const sinEq = sind(eqAng);
        double const cosEq = cosd(eqAng);
        double const sinPol = sind(polarAng);
        double const cosPol = cosd(polarAng);

        Eigen::Matrix3d rotMat;
        rotMat <<
             (sinEq * sinEq) + (cosEq * cosEq * cosPol),  - sinEq * cosEq * (1 - cosPol),              cosEq * sinPol,
            - sinEq * cosEq * (1 - cosPol),                (cosEq * cosEq) + (sinEq * sinEq * cosPol), sinEq * sinPol,
            - cosEq * sinPol,                             - sinEq * sinPol,                            cosPol;

        toVec = rotMat * fromVec;
    }
    
}
