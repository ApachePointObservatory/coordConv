#pragma once

#include <string>
#include <iostream>
#include "coordConv/mathUtils.h"

namespace coordConv {

    /**
    Position, velocity and time
    
    Requirements:
    * position is in degrees (though that is only required by setFromPair with isAngle=true).
    * velocity is in degrees/unit of time.
    However, as used within this package, time is always TAI (MJD, seconds).
    */
    class PVT {
    public:
        double pos, vel, t;

        /**
        Construct from specified position, velocity and time
        
        @param pos: position (deg)
        @param vel: velocity (deg/sec)
        @param t:   TAI date (MJD, sec)
        */
        explicit PVT(double pos, double vel, double t)
        :
            pos(pos),
            vel(vel),
            t(t)
        { };
        
        /**
        Construct a null PVT
        */
        explicit PVT()
        :
            pos(DoubleNaN),
            vel(DoubleNaN),
            t(DoubleNaN)
        { };
        
        virtual ~PVT() {};
        
        /// Return a copy
        PVT copy() const {
            return PVT(pos, vel, t);
        }
        
        /// Return a copy with a specified time
        PVT copy(double t) const {
            return PVT(getPos(t), vel, t);
        }
        
        /// Return the position at the specified time; return NaN if unknown
        double getPos(double t) const {
            return pos + (vel * (t - this->t));
        }
        
        /// Set PVT invalid at the specified time (which defaults to NaN)
        void invalidate(double t = DoubleNaN) {
            pos = DoubleNaN;
            vel = DoubleNaN;
            this->t = t;
        }
        
        /// Are all values finite? (Does it have finite pos, vel and t)?
        bool isfinite() const {
            return std::isfinite(pos) &&  std::isfinite(vel) && std::isfinite(t);
        }

        PVT & operator+=(PVT const &rhs) {
            pos += rhs.getPos(t);
            vel += rhs.vel;
            return *this;
        }

        PVT & operator-=(PVT const &rhs) {
            pos -= rhs.getPos(t);
            vel -= rhs.vel;
            return *this;
        }

        PVT & operator+=(double const &rhs) {
            pos += rhs;
            return *this;
        }

        PVT & operator-=(double const &rhs) {
            pos -= rhs;
            return *this;
        }

        PVT & operator*=(double const &rhs) {
            pos *= rhs;
            vel *= rhs;
            return *this;
        }

        PVT & operator/=(double const &rhs) {
            pos /= rhs;
            vel /= rhs;
            return *this;
        }


        PVT operator+(PVT const &rhs) const {
            PVT retVal(*this);
            retVal += rhs;
            return retVal;
        }

        PVT operator-(PVT const &rhs) const {
            PVT retVal(*this);
            retVal -= rhs;
            return retVal;
        }

        PVT operator+(double const &rhs) const {
            PVT retVal(*this);
            retVal += rhs;
            return retVal;
        }

        PVT operator-(double const &rhs) const {
            PVT retVal(*this);
            retVal -= rhs;
            return retVal;
        }

        PVT operator*(double const &rhs) const {
            PVT retVal(*this);
            retVal *= rhs;
            return retVal;
        }

        PVT operator/(double const &rhs) const {
            PVT retVal(*this);
            retVal /= rhs;
            return retVal;
        }
        
        PVT operator-() const {
            return PVT(-pos, -vel, t);
        }

        bool operator==(PVT const &rhs) {
            return (this->pos == rhs.pos) && (this->vel == rhs.vel) && (this->t == rhs.t);
        }

        bool operator!=(PVT const &rhs) {
            return !operator==(rhs);
        }
        
        /**
        Return a string representation
        */
        virtual std::string __repr__() const;
        
        /**
        Set from a pair of positions or angles computed at different times
        
        @param[in] posPair: pair of positions, where:
            posPair[0] is computed at time t
            posPair[1] is computed at time t + deltaT
        @param[in] t: time at which posPair[0] is computed
        @param[in] deltaT: time difference between the two positions
        @param[in] isAngle: if true then posPair values are treated as angles, in degrees,
            and the resulting velocity is computed using:
            posPair[1] - posPair[0] wrapped into the range [-180, 180)
        
        @warning if isAngle true then posPair must be in degrees
        */
        void setFromPair(double const posPair[2], double t, double deltaT, bool isAngle) {
            pos = posPair[0];
            if (isAngle) {
                vel = coordConv::wrapCtr(posPair[1] - posPair[0]) / deltaT;
            } else {
                vel = (posPair[1] - posPair[0]) / deltaT;
            }
            this->t = t;
        }
    };

    /**
    Convert cartesian coordinates to polar coordinates.

    @param[out] r: magnitude of vector (same units as "x" and "y")
    @param[out] theta: angle of vector (degrees)
               0 along x, 90 along y and in the range (-180, 180)
    @param[in] x: x component of vector (arbitrary units)
    @param[in] y: y component of vector (same units as "x")
    @param[in] tai: TAI date (MJD, sec)
    @return true if |r| is so small that theta cannot be computed and sets theta to 0
    */
    bool polarFromXY(PVT &r, PVT &theta, PVT const &x, PVT const &y, double tai);

    /**
    Convert polar coordinates to cartesian coordinates.

    @param[out] x: x component of vector (same units as "r")
    @param[out] y: y component of vector (same units as "r")
    @param[in] r: magnitude of vector (arbitrary units)
    @param[in] theta: angle of vector from x axis (degrees)
    @param[in] tai: TAI date (MJD, sec)
    */
    void xyFromPolar(PVT &x, PVT &y, PVT const &r, PVT const &theta, double tai);

    /**
    Rotate a 2-dimensional PVT vector by a given angle.

    @param[out] rotX: rotated x value
    @param[out] rotY: rotated y value
    @param[in]  x: unrotated x value
    @param[in]  y: unrotated y value
    @param[in]  ang: angle by which to rotate (deg)
    @param[in] tai: TAI date (MJD, sec)

    Changing coordinate systems:
    Given a point P whose position in coordinate system A is P_A_xy
    and another coordinate system B whose angle with respect to A is B_A_ang
    and whose position with respect to A is B_A_xy,
    then P_B_xy, the position of P in coordinate system B is:
        P_B_xy = (P_A_xy - B_A_xy) rotated by -B_A_ang
    */
    void rot2D(PVT &rotX, PVT &rotY, PVT const &x, PVT const &y, double ang, double tai);

    /**
    Compute PVT angle wrapped into range [0, 360) deg; only the pos differs
    */
    inline PVT wrapPos(
        PVT const &pvt  ///< input PVT angle (pos in deg)
    ) {
        PVT ret = pvt;
        ret.pos = coordConv::wrapPos(ret.pos);
        return ret;
    }

    /**
    Compute PVT angle wrapped into range [-180, 180) deg; only the pos differs
    */
    inline PVT wrapCtr(
        PVT const &pvt  ///< input PVT angle (pos in deg)
    ) {
        PVT ret = pvt;
        ret.pos = coordConv::wrapCtr(ret.pos);
        return ret;
    }

    std::ostream &operator<<(std::ostream &os, PVT const &pvt);

}
