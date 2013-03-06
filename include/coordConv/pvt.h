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
        
        /**
        Set from a pair of angles (in deg) computed at different times
        
        The only reason it matters that these are angles is that velocity is computed
        using anglePair[1] - anglePair[0] wrapped into the range [-180, 180)
    
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
