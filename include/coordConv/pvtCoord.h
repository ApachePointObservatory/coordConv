#pragma once

#include "coordConv/pvt.h"
#include "coordConv/coord.h"

namespace coordConv {
    
    /**
    PVT coordinates represent a coordinate moving at constant speed along a great circle
    */
    class PVTCoord {
    public:
        /**
        Construct a PVTCoord from a coord, orientation, velocity and time
        
        @param[in] coord: coordinate at time tai
        @param[in] orient: orientation of arc of motion at coord (deg)
            see Coord.offset for a diagram explaining orientation
        @param[in] vel: speed of motion along arc of great circle (deg/sec)
        @param[in] tai: TAI date of coord (MJD seconds)
        */
        explicit PVTCoord(Coord const &coord, double orient, double vel, double tai);

        /**
        Construct a PVTCoord from a pair of coords
        
        @param[in] coord0: coordinate at time tai
        @param[in] coord1: coordinate at time tai + deltaT; proper motion and radial velocity are ignored
        @param[in] tai: initial TAI date of PVTCoord (MJD seconds)
        @param[in] deltaT: TAI of coord1 - TAI of coord0
        */
        explicit PVTCoord(Coord const &coord0, Coord const &coord1, double tai, double deltaT);
        
        /**
        Construct a PVTCoord with all NaN data
        */
        explicit PVTCoord();

        virtual ~PVTCoord() {};
        
        /**
        Get orientation at initial TAI date
        */
        virtual double getOrient() const { return _orient; };
        
        /**
        Get velocity
        */
        virtual double getVel() const { return _vel; };
        
        /**
        Get initial TAI date
        */
        virtual double getTAI() const { return _tai; };
        
        /**
        Return Coord at initial TAI date
        */
        virtual Coord getCoord() const { return _coord; };
        
        /**
        Compute the coord at a specified TAI date
        
        @param[in] tai: TAI date at which to compute coord (MJD, sec)
        */
        virtual Coord getCoord(double tai) const;
        
        /**
        Retrieve spherical position
        
        @param[out] equatPVT: equatorial PVT (deg, deg/sec, TAI date)
        @param[out] polarPVT: polar PVT (deg, deg/sec, TAI date)
        @return atPole: true if so near the pole that equatorial angle could not be computed.
        */
        virtual bool getSphPVT(PVT &equatPVT, PVT &polarPVT) const;
        
        /**
        Offset a PVTCoord; see Coord::offset for a full explanation.

        @param[out] toOrient: orientation of offset arc at offset position (deg)
        @param[in] fromOrient: orientation of offset arc at this position (deg)
        @param[in] dist: offset distance as the length of the arc of a great circle (deg)
        @return offset coord
        
        @note The result is PVTCoord at tai, offset by dist at tai in direction fromOrient at tai.

        @raise runtime_error if this coord is too near a pole
        */
        virtual PVTCoord offset(PVT &toOrient, PVT const &fromOrient, PVT const &dist, double tai) const;
    
    private:
        Coord _coord;   // coordinate at initial time
        double _orient; // orientation of great circle arc at initial time (deg)
        double _vel;    // velocity along the great circle (deg/sec)
        double _tai;    // initial TAI date (MJD, seconds)
    };

}
