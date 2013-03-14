#pragma once

#include "coordConv/mathUtils.h"

namespace coordConv {

    /**
    Position, velocity and time
    
    Time is TAI (MJD seconds)
    
    These are designed to hold angles (in deg), though this is only assumed in a few cases.
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
        
        /// Return the position at the specified time; return NaN if unknown
        double getPos(double tai) const {
            return pos + (vel * (tai - t));
        }
        
        /// Is this PVT valid? (Does it have finite pos, vel and t)?
        bool isValid() const {
            return std::isfinite(t) && std::isfinite(vel) && std::isfinite(pos);
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
        
        /**
        Set from a pair of angles (in deg) computed at different times
        
        It matters that these are angles because velocity is computed using
        anglePair[1] - anglePair[0] wrapped into the range [-180, 180)
    
        @param[out] pvt: PVT to set
        @param[in] anglePair: pair of positions (angles in degrees), where:
            anglePair[0] is computed at time tai
            anglePair[1] is computed at time tai + deltaTAI
        @param[in] tai: TAI date (MJD sec)
        @param[in] deltaTAI: time difference (sec)
        */
        void setFromAnglePair(double anglePair[2], double tai, double deltaTAI) {
            pos = anglePair[0];
            vel = coordConv::wrapCtr(anglePair[1] - anglePair[0]) / deltaTAI;
            t = tai;
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

}
