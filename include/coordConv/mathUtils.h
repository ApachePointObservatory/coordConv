#pragma once

#include <cmath>
#include <limits>
#include "Eigen/Dense"
#include "coordConv/physConst.h"
#include "coordConv/mathUtils.cc"
/*
Math utilities
*/
namespace coordConv {
    
    // useful constants (especially useful in Python, since there is no other easy way to get to them)
    const double DoubleEpsilon = std::numeric_limits<double>::epsilon();
    const double DoubleMax = std::numeric_limits<double>::max();
    const double DoubleMin = std::numeric_limits<double>::min();
    const double DoubleNaN = std::numeric_limits<double>::quiet_NaN();

    /**
    Hypotenuse function (also present in C++11, which I'm not using yet)
    
    @param[in] x, y: sides of right triangle
    @return hypotenuse of right triangle
    */
    double hypot(double x, double y);

    inline bool isFinite(double val) {
        return ((val > std::numeric_limits<double>::min()) && (val < std::numeric_limits<double>::max()));
    }

    /**
    Rotate a 2-dimensional vector by a given angle.

    @param[out] rotVec: rotated vector
    @param[in]  vec: vector to rotate
    @param[in]  ang: angle by which to rotate (rad)

    Changing coordinate systems:
    Given a point P whose position in coordinate system A is P_A_xy
    and another coordinate system B whose angle with respect to A is B_A_ang
    and whose position with respect to A is B_A_xy,
    then P_B_xy, the position of P in coordinate system B is:
        P_B_xy = (P_A_xy - B_A_xy) rotated by -B_A_ang
    */
//    inline void rot2D(Eigen::Vector2d &rotVec, Eigen::Vector2d const &vec, double ang);

    /**
    Compute angle wrapped into range: 0 <= wrapped ang < 360 deg
    
    @param[in] ang: angle to wrap (deg)
    @return wrapped angle (deg)
    */
    double wrapPos(double ang);

    /**
    Compute angle wrapped into range: -180 <= wrapped ang < 180 deg
    
    @param[in] ang: angle to wrap (deg)
    @return wrapped angle (deg)
    */
    double wrapCtr(double ang);
    
    /**
    Wrap one angle to be within 180 degrees of another: 180 <= wrapped ang - nearAng < 180
    
    @param[in] ang: angle to wrap (deg)
    @param[in] nearAng: result is wrapped to be near this angle (deg)
    @return wrapped angle (deg)
    */
    inline double wrapNear(double ang, double nearAng);

    /// sine of angle in degrees
    inline double sind(double ang) { return std::sin(ang * RadPerDeg); }

    /// cosine of angle in degrees
    inline double cosd(double ang) { return std::cos(ang * RadPerDeg); }

    /// tangent of angle in degrees
    inline double tand(double ang) { return std::tan(ang * RadPerDeg); }

    /// arcsine in degrees
    inline double asind(double x) { return std::asin(x) / RadPerDeg; }

    /// arccosine in degrees
    inline double acosd(double x) { return std::acos(x) / RadPerDeg; }

    /// arctangent in degrees
    inline double atand(double x) { return std::atan(x) / RadPerDeg; }

    /// arctangent2 in degrees
    inline double atan2d(double x, double y) { return std::atan2(x, y) / RadPerDeg; }

    /**
    Convert cartesian coordinates to polar coordinates.

    @param[out] r: magnitude of vector (same units as "x" and "y")
    @param[out] theta: angle of vector (degrees)
               0 along x, 90 along y and in the range (-180, 180)
    @param[in] x: x component of vector (arbitrary units)
    @param[in] y: y component of vector (same units as "x")
    @return true if |r| is so small that theta cannot be computed and sets theta to 0
    */
    bool polarFromXY (double &r, double &theta, double x, double y);

    /**
    Convert polar coordinates to cartesian coordinates.

    @param[out] x: x component of vector (same units as "r")
    @param[out] y: y component of vector (same units as "r")
    @param[in] r: magnitude of vector (arbitrary units)
    @param[in] theta: angle of vector from x axis (degrees)
    */
    void xyFromPolar(double &x, double &y, double r, double theta);
        
    /**
    Compute a rotation matrix given an axis and rotation angle
    
    @param[out] rotMat: rotation matrix
    @param[in] axis: axis of rotation (magnitude is ignored, but must be finite and nonzero)
    @param[in] rotAngle: rotation angle (deg)
    */
    void computeRotationMatrix(Eigen::Matrix3d &rotMat, Eigen::Vector3d const &axis, double rotAngle);
}
