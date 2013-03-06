#pragma once

#include "Eigen/Dense"
#include "coordConv/pvt.h"
#include "coordConv/coord.h"

namespace coordConv {
    
    /**
    PVT coordinates represent a coordinate at a particular time with a known velocity
    */
    class PVTCoord {
    public:
        /**
        Construct a PVTCoord from a pair of coords
        
        @param[in] coord0: coordinate at time tai
        @param[in] coord1: coordinate at time tai + deltaT; proper motion and radial velocity are ignored
        @param[in] tai: TAI date of PVTCoord (MJD seconds)
        @param[in] deltaT: TAI of coord1 - TAI of coord0
        
        @warning: if deltaT is large the computed velocity will be inaccurate.
        */
        explicit PVTCoord(Coord const &coord0, Coord const &coord1, double tai, double deltaT);

        /**
        Construct a PVTCoord with zero velocity from a single Coord
        
        @param[in] coord: coordinate at time tai
        @param[in] tai: time of PVTCoord (TAI, MJD seconds)
        */
        explicit PVTCoord(Coord const &coord, double tai);
        
        virtual ~PVTCoord() {};
        
        /**
        Get the TAI date of this coordinate
        */
        virtual double getTAI() const { return _tai; };
        
        /**
        Retrieve the coord at specified time
        
        @warning: if tai is not near the date of this coordinate (see getTAI)
        the result will be inaccurate. The point of allowing different times is to 
        compute velocity at the date of this coordinate, not to extrapolate position
        over long periods of time.
        */
        virtual Coord getCoord(double tai) const { return Coord(_pos + (_vel * (tai - _tai)), _pm); };
        
        /**
        Retrieve spherical position
        
        @param[out] equatPVT: equatorial PVT (deg, deg/sec, TAI date)
        @param[out] polarPVT: polar PVT (deg, deg/sec, TAI date)
        @return atPole: true if so near the pole that equatorial angle could not be computed.
        */
        virtual bool getSphPVT(PVT &equatPVT, PVT &polarPVT) const;
    
    private:
        double _tai; // TAI date (MJD, seconds)
        Eigen::Vector3d _pos;   // position (AU)
        Eigen::Vector3d _vel;   // d_pos/dt (AU/sec)
        Eigen::Vector3d _pm;    // proper motion and radial velocity (AU/year)
    };

}
