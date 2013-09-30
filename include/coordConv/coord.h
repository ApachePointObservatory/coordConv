#pragma once

#include <limits>
#include <map>
#include <string>
#include "Eigen/Dense"

namespace coordConv {
    
    const double MinParallax = 1e-7; // minimum parallax (arcsec)
    
    /**
    Coordinates represent target position and proper motion.
    
    The coordinate system is always right handed, so azimuth is 0 south, 90 east
    (and hour angle is negated, though it is not visible in very much of the API).
    
    Access is available as spherical coordinates and cartesian vectors.
    
    If parallax < MinParallax / 0.9 then isInfinity() returns true and parallax is reported as 0.
    Having a lower limit prevents vector operations from overflowing.
    */
    class Coord {
    public:
        /**
        Construct a Coord from spherical position
        
        @param[in] equatAng: equatorial angle (e.g. RA, Long, Az) (degrees)
        @param[in] polarAng: polar angle (e.g. Dec, Latitude, Alt) (degrees)
        @param[in] parallax: parallax (arcsec)
        */
        explicit Coord(double equatAng, double polarAng, double parallax=0);

        /**
        Construct a Coord from spherical position and proper motion
        
        @param[in] equatAng: equatorial angle (e.g. RA, Long, Az) (degrees)
        @param[in] polarAng: polar angle (e.g. Dec, Latitude, Alt) (degrees)
        @param[in] parallax: parallax (arcsec)
        @param[in] equatPM: equatorial proper motion (arcsec/century);
            this is dEquatAng/dt, so it gets large near the pole
        @param[in] polarPM: polar proper motion (arcsec/century)
        @param[in] radVel: radial velocity (km/sec, positive receding)
        */
        explicit Coord(double equatAng, double polarAng, double parallax, double equatPM, double polarPM, double radVel);
        
        /**
        Construct a Coord from cartesian position
        
        @param[in] pos: cartesian position (AU)
        @warning: distance is not constrained
        */
        explicit Coord(Eigen::Vector3d const &pos);
        
        /**
        Construct a Coord from cartesian position and velocity
        
        @param[in] pos: cartesian position (AU)
        @param[in] vel: cartesian velocity (AU/year)
        @warning: distance is not constrained
        */
        explicit Coord(Eigen::Vector3d const &pos, Eigen::Vector3d const &vel);
        
        /**
        Construct a Coord with unknown position and proper motion
        */
        explicit Coord();
        
        ~Coord() {};
    
        /**
        Return true if object is so far away that it is considered to be at infinity

        @return true if distance is within 90% of MinParallax (includes some slop to handle coordinate conversion
        that reduce distance slightly).
        */
        bool atInfinity() const { return _atInfinity; }

        /**
        Return true if so near the pole that equatorial angles cannot be computed
        */
        bool atPole() const { return _atPole; }
        
        /**
        Get distance in AU
        
        @return distance, in AU; if atInfinity() then the value is not inf, but will often be approximately
        AUPerParsec/MinParallax.
        */
        double getDist() const { return _dist; }
        
        /**
        Get parallax in arcsec
        
        @return parallax (arcsec), or 0 if atInfinity()
        */
        double getParallax() const;
    
        /**
        Get spherical position
        
        @param[out] equatAng: equatorial angle (e.g. RA, Long, Az) (degrees)
        @param[out] polarAng: polar angle (e.g. Dec, Latitude, Alt) (degrees)
        @return atPole; if true then equatAng is arbitrarily set to 0
        */
        bool getSphPos(double &equatAng, double &polarAng) const;
        
        /**
        Get proper motion
        
        @param[out] equatPM: equatorial component of proper motion (e.g. dRA/dt) (arcsec/century)
        @param[out] polarPM: polar component of proper motion (e.g. dDec/dt) (arcsec/century)
        @return atPole; if true then equatPM and polarPM are arbitrarily set to 0
        */
        bool getPM(double &equatPM, double &polarPM) const;
        
        /**
        Get radial velocity
        
        @return radVel: radial velocity (km/sec, positive receding)
        @warning if atInfinity then radVel may be surprisingly changed by a coordinate transformation;
            you may want to report toRadVel = fromRadVel in that situation
        */
        double getRadVel() const;

        /**
        Get the cartesian position and velocity
        
        @return cartesian position (AU)
        */
        Eigen::Vector3d const getVecPos() const { return _pos; }

        /**
        Get the cartesian proper motion and radial velocity
        
        @return cartesian proper motion and radial velocity (AU/year)
        */
        Eigen::Vector3d const getVecPM() const { return _pm; }
        
        /**
        Test if all values are finite
        */
        bool isfinite() const;

        /**
        Compute the angular separation from another coord
        
        @return angular separation (deg)
        */
        double angularSeparation(Coord const &coord) const;

        /**
        Compute the orientation of a great circle offset to another coord
        
        In detail: computes the orientation at this point of a great circle connecting this coord
        to another coord. The orientation is 0 if the great circle lies along the direction of
        increasing equatorial angle, 90 if it lies along the direction increasing polar angle.
        
        @return orientation (deg), or NaN if the two coords are too close together
        */
        double orientationTo(Coord const &coord) const;
        
        /**
        Compute a new coord offset from this coord along the arc of a great circle
        
        @param[out] toOrient: orientation of offset arc at offset position (deg)
        @param[in] fromOrient: orientation of offset arc at this position (deg)
        @param[in] dist: offset distance as the length of the arc of a great circle (deg)
        @return offset coord

        @raise runtime_error if this coord is too near a pole
    
        This diagram may help:

                                                                .
                                                          .      ⎞ toOrient
                                                   *--------------> dir. of increasing equatorial angle at offset coord
                                           .    offset coord
                                  .
                        .
               .         ⎞ fromOrient
         *----------------> dir. of increasing equatorial angle at this coord
        this coord
        */
        Coord offset(double &toOrient, double fromOrient, double dist) const;
        
        /**
        Return a string representation
        */
        virtual std::string __repr__() const;

    private:
        Eigen::Vector3d _pos;   // vector position (AU)
        Eigen::Vector3d _pm;    // vector proper motion and radial velocity (AU/year)
        double _dist;           // distance (AU); a cache of _pos.norm()
        bool _atInfinity;       // true if distance far enough; a cached value
        bool _atPole;           // true if very near the pole; a cached value
        
        /**
        Set _pos from spherical position; used by several constructors
        
        @param[in] equatAng: equatorial angle (e.g. RA, Long, Az) (degrees)
        @param[in] polarAng: polar angle (e.g. Dec, Latitude, Alt) (degrees)
        @param[in] parallax: parallax (arcsec)

        @warning Does not set _pm
        */
        void _setPosFromSph(double equatAng, double polarAng, double parallax);
        
        /**
        Set various cached information
        
        _pos must be set before you call this; _pm need not be set.
        */
        void _setCache();
    };

    inline bool operator==(Coord const &lhs, Coord const &rhs) {
        return (lhs.getVecPos() == rhs.getVecPos()) && (lhs.getVecPM() == rhs.getVecPM());
    }

    inline bool operator!=(Coord const &lhs, Coord const &rhs) {
        return !operator==(lhs, rhs);
    }
    
    std::ostream &operator<<(std::ostream &out, Coord const &coord);

}
