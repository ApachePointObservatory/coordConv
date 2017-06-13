#pragma once

#include "coordConv/pvt.h"
#include "coordConv/coord.h"

namespace coordConv {
    
    /**
    A coordinate moving at a constant cartesian velocity

    Primarily intended for computing instantaneous velocity by comparing two Coord at nearby time.
    */
    class PVTCoord {
    public:
        /**
        Construct a PVTCoord from a coord, vector velocity and time
        
        @param[in] coord  coordinate at time tai
        @param[in] vel  cartesian velocity (deg/sec)
        @param[in] tai  TAI date of coord (MJD seconds)

        @throw std::runtime_error if coord is at pole at time tai and vel nonzero.
        */
        explicit PVTCoord(Coord const &coord, Eigen::Vector3d const &vel, double tai);

        /**
        Construct a PVTCoord from a pair of coords
        
        @param[in] coord0  coordinate at time tai
        @param[in] coord1  coordinate at time tai + deltaT; proper motion and radial velocity are ignored
        @param[in] tai  initial TAI date of PVTCoord (MJD seconds)
        @param[in] deltaT  TAI of coord1 - TAI of coord0; must be nonzero
        
        @throw std::runtime_error if deltaT = 0
        */
        explicit PVTCoord(Coord const &coord0, Coord const &coord1, double tai, double deltaT);

        /**
        Construct a PVTCoord from spherical PVTs and distance at TAI date equatPVT.t

        @param[in] equatPVT  equatorial angle (e.g. RA, Long, Az) (degrees)
        @param[in] polarPVT  polar angle (e.g. Dec, Latitude, Alt) (degrees)
        @param[in] distPVT  distance (AU); if invalid (!distPVT.isfinite()) then infinity is assumed

        @throw std::runtime_error if equatPVT.t != polarPVT.t
        */
        explicit PVTCoord(PVT const &equatPVT, PVT const &polarPVT, PVT const &distPVT=PVT());

        /**
        Construct a PVTCoord from spherical PVTs, parallax, proper motion and radial velocity
        at TAI date equatPVT.t
        
        @param[in] equatPVT  equatorial angle (e.g. RA, Long, Az) (degrees)
        @param[in] polarPVT  polar angle (e.g. Dec, Latitude, Alt) (degrees)
        @param[in] distPVT  distance (AU)
        @param[in] equatPM  equatorial proper motion (arcsec/century);
            this is dEquatAng/dt, so it gets large near the pole
        @param[in] polarPM  polar proper motion (arcsec/century)
        @param[in] radVel  radial velocity (km/sec, positive receding)

        @throw std::runtime_error if equatPVT.t != polarPVT.t
        */
        explicit PVTCoord(PVT const &equatPVT, PVT const &polarPVT, PVT const &distPVT, double equatPM, double polarPM, double radVel);
        
        /**
        Construct a PVTCoord with all NaN data
        */
        explicit PVTCoord();
        
        /**
        Copy PVTCoord at the specified date

        @param[in] tai  TAI date at which to retrieve position (MJD, sec)
        */
        PVTCoord copy(double tai) const;

        /**
        Get orientation at initial TAI date
        */
//        double getOrient() const { return _orient; };
        
        /**
        Get velocity
        */
        Eigen::Vector3d getVel() const { return _vel; };
        
        /**
        Get initial TAI date
        */
        double getTAI() const { return _tai; };
        
        /**
        Return Coord at initial TAI date
        */
        Coord getCoord() const { return _coord; };
        
        /**
        Compute the coord at a specified TAI date
        
        @param[in] tai  TAI date at which to compute coord (MJD, sec)
        */
        Coord getCoord(double tai) const;

        /**
        Return True if all values are finite
        */
        bool isfinite() const;
       
        /**
        Get spherical position
        
        The returned velocities are the dEquat/dt and dPolar/dt at the pvtCoord's TAI date
        
        @param[out] equatPVT  equatorial PVT (deg, deg/sec, TAI date)
        @param[out] polarPVT  polar PVT (deg, deg/sec, TAI date)
        @return atPole: true if so near the pole that equatorial angle could not be computed.
        */
        bool getSphPVT(PVT &equatPVT, PVT &polarPVT) const;

        /**
        Get distance in AU
        
        @return distance, in AU; if getCoord().atInfinity() then the value is not inf,
        but will often be approximately AUPerParsec/MinParallax.
        */
        PVT getDistance() const;

        /**
        Compute the angular separation from another PVTCoord

        The separation is computed at the date of this pvtCoord

        @param[in] pvtCoord  PVT coord to which to measure angular separation
        
        @return angular separation (deg)
        */
        PVT angularSeparation(PVTCoord const &pvtCoord) const;

        /**
        Compute the orientation of a great circle offset to another PVTCoord

        The orientation is computed at the date of this pvtCoord

        @param[in] pvtCoord  PVT coord to which to measure orientation
        
        In detail: computes the orientation at this point of a great circle connecting this PVT coord
        to another PVT coord. The orientation is 0 if the great circle lies along the direction of
        increasing equatorial angle, 90 if it lies along the direction increasing polar angle.
        
        @return orientation (deg), or NaN if the angular separation is too near 0 or 180 at tai
        */
        PVT orientationTo(PVTCoord const &pvtCoord) const;
        
        /**
        Offset a PVTCoord by a specified distance in a specified direction; see Coord::offset for a full explanation.

        @param[out] toOrient  orientation of offset arc at offset position (deg)
        @param[in] fromOrient  orientation of offset arc at this position (deg)
        @param[in] dist  offset distance as the length of the arc of a great circle (deg)
        @return offset PVT coord

        @throw std::runtime_error if this PVTCoord is too near a pole to be offset
        */
        PVTCoord offset(PVT &toOrient, PVT const &fromOrient, PVT const &dist) const;

        bool operator==(PVTCoord const &rhs) {
            return (this->getCoord()  == rhs.getCoord())
                && (this->getVel()    == rhs.getVel())
                && (this->getTAI()    == rhs.getTAI());
        }

        bool operator!=(PVTCoord const &rhs) {
            return !operator==(rhs);
        }

        /**
        Return a string representation
        */
        std::string __repr__() const;
    
    private:
        Coord _coord;   // coordinate at initial time
        Eigen::Vector3d _vel;    // vector velocity (AU/sec)
        double _tai;    // initial TAI date (MJD, seconds)

        /**
        Set fields based on a pair of coords
        
        @param[in] coord0  coordinate at time tai
        @param[in] coord1  coordinate at time tai + deltaT; proper motion and radial velocity are ignored
        @param[in] tai  initial TAI date of PVTCoord (MJD seconds)
        @param[in] deltaT  TAI of coord1 - TAI of coord0; must be nonzero
        
        @throw std::runtime_error if:
        - deltaT = 0

        @note this method can go away once we switch to C++11 and can have one constructor call another.
        */
        void _setFromCoordPair(Coord const &coord0, Coord const &coord1, double tai, double deltaT);

    };


    std::ostream &operator<<(std::ostream &os, PVTCoord const &pvtCoord);

}
