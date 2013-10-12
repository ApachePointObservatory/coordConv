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
            see Coord.offset for an explanation of orientation
        @param[in] vel: speed of motion along arc of great circle (deg/sec)
        @param[in] tai: TAI date of coord (MJD seconds)
        */
        explicit PVTCoord(Coord const &coord, double orient, double vel, double tai);

        /**
        Construct a PVTCoord from a pair of coords
        
        @param[in] coord0: coordinate at time tai
        @param[in] coord1: coordinate at time tai + deltaT; proper motion and radial velocity are ignored
        @param[in] tai: initial TAI date of PVTCoord (MJD seconds)
        @param[in] deltaT: TAI of coord1 - TAI of coord0; must be nonzero
        @param[in] defOrient: default orientation (deg); if the velocity is so low that orientation
            cannot be computed then orientation=defOrient and vel=0
            See Coord.offset for an explanation of orientation
        
        @raise std::runtime_error if deltaT = 0
        */
        explicit PVTCoord(Coord const &coord0, Coord const &coord1, double tai, double deltaT, double defOrient=0);

        /**
        Construct a PVTCoord from spherical PVTs

        @param[in] equatAng: equatorial angle (e.g. RA, Long, Az) (degrees)
        @param[in] polarAng: polar angle (e.g. Dec, Latitude, Alt) (degrees)
        @param[in] tai: date at which to evaluate equatAng and polarAng, and initial TAI date of PVTCoord (MJD seconds)
        @param[in] parallax: parallax (arcsec)
        @param[in] defOrient: default orientation (deg); if the velocity is so low that orientation
            cannot be computed then orientation=defOrient and vel=0
            See Coord.offset for an explanation of orientation

        The provided velocities are treated as the instantaneous dEquat/dt and dPolar/dt along the great circle arc, so:
        - orientation = atan2d(equat space vel, polar vel)
        - velocity = hypot(equat space vel, polar vel)
        where equat space vel = equat vel * cosd(polar pos at tai)

        @raise std::runtime_error if too near pole and equat or polar abs(vel) > DoubleEpsilon
        */
        explicit PVTCoord(PVT const &equatPVT, PVT const &polarPVT, double tai, double parallax=0, double defOrient=0);

        /**
        Construct a PVTCoord from spherical PVTs and proper motion
        
        @param[in] equatAng: equatorial angle (e.g. RA, Long, Az) (degrees)
        @param[in] polarAng: polar angle (e.g. Dec, Latitude, Alt) (degrees)
        @param[in] tai: date at which to evaluate equatAng and polarAng, and initial TAI date of PVTCoord (MJD seconds)
        @param[in] parallax: parallax (arcsec)
        @param[in] equatPM: equatorial proper motion (arcsec/century);
            this is dEquatAng/dt, so it gets large near the pole
        @param[in] polarPM: polar proper motion (arcsec/century)
        @param[in] radVel: radial velocity (km/sec, positive receding)
        @param[in] defOrient: default orientation (deg); if the velocity is so low that orientation
            cannot be computed then orientation=defOrient and vel=0
            See Coord.offset for an explanation of orientation

        The PVT velocities are treated as the instantaneous dEquat/dt and dPolar/dt along the great circle arc:
        - orientation = atan2d(equat space vel, polar vel)
        - velocity = hypot(equat space vel, polar vel)
        where equat space vel = equat vel * cosd(polar pos at tai)

        @raise std::runtime_error if too near pole and equat or polar abs(vel) > DoubleEpsilon
        */
        explicit PVTCoord(PVT const &equatPVT, PVT const &polarPVT, double tai, double parallax, double equatPM, double polarPM, double radVel, double defOrient=0);
        
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
        Return True if all values are finite
        */
        bool isfinite() const;
       
        /**
        Retrieve spherical position at a specified TAI date
        
        The returned velocities are the instantaneous dEquat/dt and dPolar/dt along the great circle arc
        
        @param[out] equatPVT: equatorial PVT (deg, deg/sec, TAI date)
        @param[out] polarPVT: polar PVT (deg, deg/sec, TAI date)
        @param[in] tai: TAI date at which to retrieve position (MJD, sec)
        @return atPole: true if so near the pole that equatorial angle could not be computed.
        */
        virtual bool getSphPVT(PVT &equatPVT, PVT &polarPVT, double tai) const;

        /**
        Compute the angular separation from another PVTCoord

        @param[in] coord: PVT coord to which to measure angular separation
        @param[in] tai: TAI date at which to retrieve position (MJD, sec)
        
        @return angular separation (deg)
        */
        PVT angularSeparation(PVTCoord const &pvtCoord, double tai) const;

        /**
        Compute the orientation of a great circle offset to another PVTCoord

        @param[in] pvtCoord: PVT coord to which to measure orientation
        @param[in] tai: TAI date at which to retrieve position (MJD, sec)
        
        In detail: computes the orientation at this point of a great circle connecting this PVT coord
        to another PVT coord. The orientation is 0 if the great circle lies along the direction of
        increasing equatorial angle, 90 if it lies along the direction increasing polar angle.
        
        @return orientation (deg), or NaN if the two PVTCoord are too close together at tai
        */
        PVT orientationTo(PVTCoord const &pvtCoord, double tai) const;
        
        /**
        Offset a PVTCoord; see Coord::offset for a full explanation.

        @param[out] toOrient: orientation of offset arc at offset position (deg)
        @param[in] fromOrient: orientation of offset arc at this position (deg)
        @param[in] dist: offset distance as the length of the arc of a great circle (deg)
        @param[in] tai: TAI date at which to compute offset (MJD, sec)
        @return offset PVT coord
        
        @note The result is PVTCoord at tai, offset by dist at tai in direction fromOrient at tai.

        @raise std::runtime_error if this PVT coord is too near a pole at tai
        */
        virtual PVTCoord offset(PVT &toOrient, PVT const &fromOrient, PVT const &dist, double tai) const;

        bool operator==(PVTCoord const &rhs) {
            return (this->getCoord()  == rhs.getCoord())
                && (this->getOrient() == rhs.getOrient())
                && (this->getVel()    == rhs.getVel())
                && (this->getTAI()    == rhs.getTAI());
        }

        bool operator!=(PVTCoord const &rhs) {
            return !operator==(rhs);
        }

        /**
        Return a string representation
        */
        virtual std::string __repr__() const;
    
    private:
        Coord _coord;   // coordinate at initial time
        double _orient; // orientation of great circle arc at initial time (deg)
        double _vel;    // velocity along the great circle (deg/sec)
        double _tai;    // initial TAI date (MJD, seconds)

        /**
        Set _orient and _vel from spherical info
        
        @param[in] polarPos: polar angle
        @param[in] equatVel: equatorial velocity (dEquatPos/dt)
        @param[in] polarVel: polar velocity
        @param[in] defOrient: default orientation (deg); if _coord.atPole()
            then orientation=defOrient and vel=0

        The provided velocities are treated as the instantaneous dEquat/dt and dPolar/dt along the great circle arc, so:
        - orientation = atan2d(equatSpaceVel, polarVel)
        - velocity = hypot(equatSpaceVel, polarVel)
        where equatSpaceVel = equatVel * cosd(polarPos)

        @raise std::runtime_error if _coord.atPole() and |equatVel| or |polarVel| > DoubleEpsilon.
        */
        void _setOrientVelFromSph(double polarPos, double equatVel, double polarVel, double defOrient);
    };


    std::ostream &operator<<(std::ostream &os, PVTCoord const &pvtCoord);

}
