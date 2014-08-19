#pragma once

#include "Eigen/Dense"

namespace coordConv {

    /**
    Rotate a 3-vector described by equatorial and polar angles, as follows:
    The plane of rotation contains the z axis and a line in the x-y plane
    at angle eqAng from the x axis towards y. The amount of rotation is
    angle polarAng from the z axis towards the line in the x-y plane.

    @param[out] toVec  rotated 3-vector
    @param[in] fromVec  input 3-vector
    @param[in] eqAng  angle of line in x-y plane (from x to y);
               the plane of rotation includes this line and z
    @param[in] polarAng  angle of rotation (from z axis to line in x-y plane)
    */
    void rotEqPol(
        Eigen::Vector3d &toVec,
        Eigen::Vector3d const &fromVec,
        double eqAng,
        double polarAng
    );
    
}
