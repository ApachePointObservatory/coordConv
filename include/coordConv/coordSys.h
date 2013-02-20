#pragma once

#include <limits>
#include <map>
#include <string>
#include "boost/shared_ptr.hpp"
#include "coordConv/site.h"
#include "coordConv/coord.h"

namespace coordConv {
    
    /**
    Abstract base class for coordinate systems
    */
    class CoordSys {
    public:
        ///< construct a CoordSys
        explicit CoordSys(std::string const &name, double date) : _name(name), _date() { setDate(date); };
        
        ///< destructor
        virtual ~CoordSys() { };
        
        ///< return a copy with the same date
        virtual boost::shared_ptr<CoordSys> clone() const = 0;
        
        ///< return a copy with a specified date
        virtual boost::shared_ptr<CoordSys> clone(double date) const = 0;

        ///< get the name of this coordinate system
        std::string getName() const { return _name; };
        
        ///< get the date of this coordinate system
        virtual double getDate() const { return _date; };
        
        ///< set the date of this coordinate system
        virtual void setDate(double date) { _date = date; };
        
        ///< convert a coordinate to ICRS at date of observation 2000.0 from this system at this date
        virtual Coord toICRS(Coord const &coord, Site const &site) const = 0;
        
        ///< convert a coordinate from ICRS at date of observation J2000 to this system at this date
        virtual Coord fromICRS(Coord const &coord, Site const &site) const = 0;
        
        ///< convert a coordinate from another coordinate system to this system
        virtual Coord convertFrom(CoordSys const &fromCoordSys, Coord const &fromCoord, Site const &site) const;
        
        ///< get string representation
        virtual std::string asString() const;
    
    protected:
        std::string _name;  /// name of coordinate system
        double _date;       /// date of coordinate system (units depend on coordinate system)
    };

    /**
    ICRS RA, Dec; date is Julian years
    */
    class ICRSCoordSys: public CoordSys {
    public:
        /**
        Construct an ICRSCoordSys
        
        @param[in] date: date of observation in Julian years
        */
        explicit ICRSCoordSys(double date=2000.0);
        virtual ~ICRSCoordSys() {};
        virtual boost::shared_ptr<CoordSys> clone() const;
        virtual boost::shared_ptr<CoordSys> clone(double date) const;
        virtual Coord fromICRS(Coord const &coord, Site const &site) const;
        virtual Coord toICRS(Coord const &coord, Site const &site) const;
    };
    
    /**
    FK5 RA, Dec; date is Julian years, and is both the date of observation and the date of equinox
    */
    class FK5CoordSys: public CoordSys {
    public:
        /**
        Construct an FK5CoordSys
        
        @param[in] date: date of equinox and date of observation in Julian years
        */
        explicit FK5CoordSys(double date=2000.0);
        virtual ~FK5CoordSys() {};
        virtual boost::shared_ptr<CoordSys> clone() const;
        virtual boost::shared_ptr<CoordSys> clone(double date) const;
        virtual void setDate(double date);
        virtual Coord fromICRS(Coord const &coord, Site const &site) const;
        virtual Coord toICRS(Coord const &coord, Site const &site) const;

    private:
        Eigen::Matrix3d _to2000PrecMat; /// precession matrix from date to J2000.0
    };
    
    /**
    FK4 RA, Dec; date is Besselian years, and is both the date of observation and the date of equinox
    
    @warning: the FK4 system has significant fictitious proper motion. Coords will be treated as fixed
    (the fictitious proper motion removed) if you specify proper motion AND radial velocity as zero.
    If any component of propoer motion or radial velocity is nonzero, then all components are treated
    as correct. Thus it is usually safest not to specify radial velocity for FK4 targets.
    */
    class FK4CoordSys: public CoordSys {
    public:
        /**
        Construct an FK4CoordSys
        
        @param[in] date: date of equinox and date of observation in Besselian years
        */
        explicit FK4CoordSys(double date=1950.0);
        virtual ~FK4CoordSys() {};
        virtual boost::shared_ptr<CoordSys> clone() const;
        virtual boost::shared_ptr<CoordSys> clone(double date) const;
        virtual void setDate(double date);
        virtual Coord fromICRS(Coord const &coord, Site const &site) const;
        virtual Coord toICRS(Coord const &coord, Site const &site) const;

    private:
        Eigen::Vector3d _eTerms;
        Eigen::Matrix3d _From1950PrecMat, _To1950PrecMat;
    };

    /**
    IAU 1958 Galactic longitude, latitude; date is Julian years

    References:
      P.T. Wallace, slaEqGal, a SLALIB subroutine; Starlink, RGO
      Blaauw et al, Mon.Not.R.Astron.Soc.,121,123 (1960)
    */
    class GalCoordSys: public CoordSys {
    public:
        /**
        Construct a GalCoordSys
        
        @param[in] date: date of observation in Julian years
        */
        explicit GalCoordSys(double date=2000.0);
        virtual ~GalCoordSys() {};
        virtual boost::shared_ptr<CoordSys> clone() const;
        virtual boost::shared_ptr<CoordSys> clone(double date) const;
        virtual Coord fromICRS(Coord const &coord, Site const &site) const;
        virtual Coord toICRS(Coord const &coord, Site const &site) const;
    };

    /**
    Apparent Geocentric RA, Dec coordinates; date is TDB date in Julian years
    
    Technically the date should be TDB, but TT will always do. Note that TT has a fixed offset from TAI
    (TT_TAI in physConst.h, in seconds).
    
    @warning:
    - Not fully accurate for solar system objects.
    - This transformation requires iteration, so it will be somewhat slow.

    The following approximations have been used:
    - The annual aberration correction is not accurate for solar system objects.
    - No correction is applied for the bending of light by sun's gravity.
    This introduces errors on the order of 0.02" at a distance of 20 degrees from the sun (Wallace, 1986)
    */    
    class AppGeoCoordSys: public CoordSys {
    public:
        /**
        Construct an AppGeoCoordSys
        
        @param[in] date: TDB date in Julian years (but TT will always do)
        @param[in] maxAge: maximum delta date (years) before setDate will update an internal cache
        */
        explicit AppGeoCoordSys(double date=std::numeric_limits<double>::quiet_NaN(), double maxAge=1.5e-5);
        virtual ~AppGeoCoordSys() {};
        virtual boost::shared_ptr<CoordSys> clone() const;
        virtual boost::shared_ptr<CoordSys> clone(double date) const;
        virtual void setDate(double date);
        
        /// return maximum delta date (years) before setDate will update the internal cache
        virtual double getMaxAge() const { return _maxAge; };
        virtual Coord fromICRS(Coord const &coord, Site const &site) const;
        virtual Coord toICRS(Coord const &coord, Site const &site) const;

    private:
        double _maxAge;             ///< maximum date differential to reuse cache (years)
	    double _cachedDate;         ///< date at which data computed; 0 if never computed
	    double _pmSpan;             ///< time over which to correct for proper motion (Julian years)
	    Eigen::Vector3d _bcPos;     ///< barycentric position of Earth (au)
	    Eigen::Vector3d _hcDir;     ///< heliocentric position of Earth (unit vector)
	    double _gravRad;            ///< gravitational radius of sun * 2 / sun-earth distance
	    Eigen::Vector3d _bcBeta;    ///< barycentric velocity of the Earth (c)
	    double _gammaI;             ///< sqrt(1 - bcBeta^2)
	    Eigen::Matrix3d _pnMat;     ///< precession/nutation matrix
    };

    /**
    Apparent Topocentric Az,Alt; date is TAI (MJD, seconds)
    */
    class AppTopoCoordSys: public CoordSys {
    public:
        /**
        Construct an AppTopoCoordSys
        
        @param[in] date: date as TAI (MJD, seconds)
        */
        explicit AppTopoCoordSys(double date=std::numeric_limits<double>::quiet_NaN());
        virtual ~AppTopoCoordSys() {};
        virtual boost::shared_ptr<CoordSys> clone() const;
        virtual boost::shared_ptr<CoordSys> clone(double date) const;
        virtual void setDate(double date) { setDate(date, false); };
        /**
        Set the date and optionally freeze cached apparent geocentric data
        
        @param[in] date: date as TAI (MJD, seconds)
        @param[in] freezeCache: if true, do not update cached apparent geocentric data
        
        The reason for freezeCache is to support computing velocity by computing position at two nearby times.
        Set freezeCache true for the second computation to avoid an unexpected cache update causing
        a mis-computation of velocity.
        
        @raise std::runtime_error if freezeCache true and delta date > AppGeoCoordSys's default MaxAge * 2;
        this is only intended to catch gross errors
        */
        virtual void setDate(double date, bool freezeCache);
        virtual Coord fromICRS(Coord const &coord, Site const &site) const;
        virtual Coord toICRS(Coord const &coord, Site const &site) const;
        virtual Coord fromAppGeo(Coord const &coord, Site const &site) const;
        virtual Coord toAppGeo(Coord const &coord, Site const &site) const;

    private:
        AppGeoCoordSys _appGeoCoordSys;
    };

    /**
    Observed Az,Alt (Apparent Topocentric with atmospheric refraction); date is TAI (MJD, seconds)
    */
    class ObsCoordSys: public CoordSys {
    public:
        /**
        Construct an ObsCoordSys
        
        @param[in] date: date as TAI (MJD, seconds)
        */
        explicit ObsCoordSys(double date=std::numeric_limits<double>::quiet_NaN());
        virtual ~ObsCoordSys() {};
        virtual boost::shared_ptr<CoordSys> clone() const;
        virtual boost::shared_ptr<CoordSys> clone(double date) const;
        /**
        Set the date and optionally freeze cached apparent geocentric data
        
        @param[in] date: date as TAI (MJD, seconds)
        @param[in] freezeCache: if true, do not update cached apparent geocentric data
        
        The reason for freezeCache is to support computing velocity by computing position at two nearby times.
        Set freezeCache true for the second computation to avoid an unexpected cache update causing
        a mis-computation of velocity.
        
        @raise std::runtime_error if freezeCache true and delta date > AppGeoCoordSys's default MaxAge * 2;
        this is only intended to catch gross errors
        */
        virtual void setDate(double date) { setDate(date, false); };
        virtual void setDate(double date, bool freezeCache);
        virtual Coord fromICRS(Coord const &coord, Site const &site) const;
        virtual Coord toICRS(Coord const &coord, Site const &site) const;
        virtual Coord fromAppTopo(Coord const &coord, Site const &site) const;
        virtual Coord toAppTopo(Coord const &coord, Site const &site) const;

    private:
        AppTopoCoordSys _appTopoCoordSys;
    };
    
    boost::shared_ptr<CoordSys> makeCoordSys(std::string const &name, double date);
}
