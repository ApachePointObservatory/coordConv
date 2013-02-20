#pragma once

#include "Eigen/Dense"

namespace coordConv {

    /**
    Rotate a 3-vector, first about the x axis, then about the y axis.
    Warning: toVec cannot replace fromVec

    Inputs:
    @param[out] toVec: rotated 3-vector
    @param[in] fromVec: input 3-vector
    @param[in] xAng: angle about x axis (deg); positive rotation is from y to z.
    @param[in] yAng: angle about original y axis (deg); positive rotation is from z to x.
    */
    void rotXY(
        Eigen::Vector3d &toVec,
        Eigen::Vector3d const &fromVec,
        double xAng,
        double yAng
    );
    
}
