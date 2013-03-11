#pragma once

#include <limits>
#include <map>
#include <string>
#include "boost/shared_ptr.hpp"
#include "coordConv/site.h"
#include "coordConv/coord.h"
#include "coordConv/pvtCoord.h"

namespace coordConv {
    
    /**
    Abstract base class for coordinate systems
    
    Subclasses must define fromFK5J2000 and toFK5J2000, and should override setDate if information cached based on date.
    */
    class CoordSys {
    public:
        /**
        Construct a CoordSys given a name and date
        */
        explicit CoordSys(std::string const &name, bool isMean, double date) :
            _name(name), _isMean(isMean), _date() { setDate(date); };
        
        ///< Destructor
        virtual ~CoordSys() { };
        
        /**
        Return a copy with the same date
        */
        virtual boost::shared_ptr<CoordSys> clone() const = 0;
        
        /**
        Return a copy with a specified date
        */
        virtual boost::shared_ptr<CoordSys> clone(double date) const = 0;

        /**
        Get the name of this coordinate system (all lowercase)
        */
        std::string getName() const { return _name; };
        
        /**
        Return true for a mean coordinate system, false for an apparent coordinate system
        */
        bool isMean() const { return _isMean; };
        
        /**
        Get the date of this coordinate system
        
        The units depend on the specific coordinate system
        */
        double getDate() const { return _date; };
        
        /**
        Set the date of this coordinate system
        
        The units depend on the specific coordinate system
        */
        virtual void setDate(double date) { _date = date; };
        
        /**
        Convert a coordinate to FK5 at date of observation J2000 from this coordinate system at this date
        
        @param[in] coord: position in this coordinate system at this date
        @param[in] site: site information
        @return position in ICRS coordinates at date of observation J2000
        */
        virtual Coord toFK5J2000(Coord const &coord, Site const &site) const = 0;
        
        /**
        Convert a coordinate from FK5 at date of observation J2000 to this system at this date
        
        @param[in] coord: position in ICRS coordinates at date of observation J2000
        @param[in] site: site information
        @return position in this coordinate system at this date
        */
        virtual Coord fromFK5J2000(Coord const &coord, Site const &site) const = 0;
        
        /**
        Convert a coordinate from another coordinate system to this system
        
        @param[in] fromCoordSys: initial coordinate system and date
        @param[in] fromCoord: initial position
        @param[in] site: site information
        @return position in this coordinate system at this date
        */
        virtual Coord convertFrom(CoordSys const &fromCoordSys, Coord const &fromCoord, Site const &site) const;
        
        /**
        Convert a PVTCoord from another coordinate system to this system
        
        @param[in] fromCoordSys: initial coordinate system and date
        @param[in] fromPVTCoord: initial PVTCoord
        @param[in] site: site information
        @return PVTCoord in this coordinate system at this date
        */
        virtual PVTCoord convertFrom(CoordSys const &fromCoordSys, PVTCoord const &fromPVTCoord, Site const &site) const;
        
        /**
        Convert a coordinate from another coordinate system to this system, including orientation
        
        @param[out] toDir: orientation in this coordinate system (deg; 0 along increasing equatorial angle, 90 along increasing polar angle)
        @param[out] scaleChange: change in scale: output delta sky/input delta sky, measured along the specified direction
        @param[in] fromCoordSys: initial coordinate system
        @param[in] fromCoord: initial position
        @param[in] fromDir: initial orientation (deg; 0 along increasing equatorial angle, 90 along increasing polar angle)
        @param[in] site: site information
        @return position in this coordinate system
        */
        virtual Coord convertFrom(double &toDir, double &scaleChange, CoordSys const &fromCoordSys, Coord const &fromCoord, double fromDir, Site const &site) const;
        
        /**
        Convert a PVTCoord from another coordinate system to this system, including orientation
        
        @param[out] toDir: orientation in this coordinate system (deg; 0 along increasing equatorial angle, 90 along increasing polar angle)
        @param[out] scaleChange: change in scale: output delta sky/input delta sky, measured along the specified direction
        @param[in] fromCoordSys: initial coordinate system
        @param[in] fromCoord: initial position
        @param[in] fromDir: initial orientation (deg; 0 along increasing equatorial angle, 90 along increasing polar angle)
        @param[in] site: site information
        @return position in this coordinate system
        */
        virtual PVTCoord convertFrom(PVT &toDir, double &scaleChange, CoordSys const &fromCoordSys, PVTCoord const &fromPVTCoord, PVT const &fromDir, Site const &site) const;
        
        /**
        Convert TAI (MJD, seconds) to a suitable date for this coordinate system
        
        @param[in] tai: TAI date (MJD, seconds)
        @return date in appropriate units for this coordinate system
        */
        virtual double dateFromTAI(double tai) const = 0;
        
        ///< get string representation
        virtual std::string asString() const;
    
    protected:
        std::string _name;  /// name of coordinate system
        bool _isMean;       /// true for mean coordinate systems
        double _date;       /// date of coordinate system (units depend on coordinate system)
    };
    
    class MeanCoordSys: public CoordSys {
    public:
        explicit MeanCoordSys(std::string const &name, double date);
        virtual ~MeanCoordSys() {};
        virtual double dateFromTAI(double tai) const;
        
        /**
        Remove the effects of proper motion and radial velocity to the specified TAI date
        
        @param[in] tai: TAI date (MJD, seconds)
        @return coord with proper motion and radial velocity (but not parallax) removed
        */
        virtual Coord removePM(Coord const &coord, double tai) const;
    };
    
    class ApparentCoordSys: public CoordSys {
    public:
        explicit ApparentCoordSys(std::string const &name, double date);
        virtual ~ApparentCoordSys() {};
        virtual double dateFromTAI(double tai) const;
    };

    /**
    ICRS RA, Dec; date is Julian years
    */
    class ICRSCoordSys: public MeanCoordSys {
    public:
        /**
        Construct an ICRSCoordSys
        
        @param[in] date: date of observation in Julian years
        */
        explicit ICRSCoordSys(double date=2000.0);
        virtual ~ICRSCoordSys() {};
        virtual boost::shared_ptr<CoordSys> clone() const;
        virtual boost::shared_ptr<CoordSys> clone(double date) const;
        virtual Coord fromFK5J2000(Coord const &coord, Site const &site) const;
        virtual Coord toFK5J2000(Coord const &coord, Site const &site) const;
    };
    
    /**
    FK5 RA, Dec; date is Julian years, and is both the date of observation and the date of equinox
    */
    class FK5CoordSys: public MeanCoordSys {
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
        virtual Coord fromFK5J2000(Coord const &coord, Site const &site) const;
        virtual Coord toFK5J2000(Coord const &coord, Site const &site) const;

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
    class FK4CoordSys: public MeanCoordSys {
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
        virtual Coord fromFK5J2000(Coord const &coord, Site const &site) const;
        virtual Coord toFK5J2000(Coord const &coord, Site const &site) const;
        virtual double dateFromTAI(double tai) const;

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
    class GalCoordSys: public MeanCoordSys {
    public:
        /**
        Construct a GalCoordSys
        
        @param[in] date: date of observation in Julian years
        */
        explicit GalCoordSys(double date=2000.0);
        virtual ~GalCoordSys() {};
        virtual boost::shared_ptr<CoordSys> clone() const;
        virtual boost::shared_ptr<CoordSys> clone(double date) const;
        virtual Coord fromFK5J2000(Coord const &coord, Site const &site) const;
        virtual Coord toFK5J2000(Coord const &coord, Site const &site) const;
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
    class AppGeoCoordSys: public ApparentCoordSys {
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
        virtual Coord fromFK5J2000(Coord const &coord, Site const &site) const;
        virtual Coord toFK5J2000(Coord const &coord, Site const &site) const;
        virtual double dateFromTAI(double tai) const;
        
        /// return maximum delta date (years) before setDate will update the internal cache
        virtual double getMaxAge() const { return _maxAge; };

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
    class AppTopoCoordSys: public ApparentCoordSys {
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
        virtual Coord fromFK5J2000(Coord const &coord, Site const &site) const;
        virtual Coord toFK5J2000(Coord const &coord, Site const &site) const;

        /**
        Set the date and optionally freeze cached apparent geocentric data
        
        @param[in] date: date as TAI (MJD, seconds)
        @param[in] freezeCache: if true, do not update cached apparent geocentric data
        
        The reason for freezeCache is to support computing velocity by computing position at two nearby times.
        Set freezeCache true for the second computation to avoid an unexpected cache update causing
        a mis-computation of velocity.
        
        @note The standard setDate, with no freezeCache, argument does not freeze the cache.
        
        @raise std::runtime_error if freezeCache true and delta date > AppGeoCoordSys's default MaxAge * 2;
        this is only intended to catch gross errors
        */
        virtual void setDate(double date, bool freezeCache);

        /**
        Convert from apparent geocentric coordinates at this date

        @param[in] coord: initial position
        @param[in] site: site information
        @return position in this coordinate system at this date
        */
        virtual Coord fromAppGeo(Coord const &coord, Site const &site) const;

        /**
        Convert to apparent geocentric coordinates at this date

        @param[in] coord: initial position
        @param[in] site: site information
        @return position in apparent geocentric coordinates at this date
        */
        virtual Coord toAppGeo(Coord const &coord, Site const &site) const;

    private:
        AppGeoCoordSys _appGeoCoordSys;
    };

    /**
    Observed Az,Alt (Apparent Topocentric with atmospheric refraction); date is TAI (MJD, seconds)
    */
    class ObsCoordSys: public ApparentCoordSys {
    public:
        /**
        Construct an ObsCoordSys
        
        @param[in] date: date as TAI (MJD, seconds)
        */
        explicit ObsCoordSys(double date=std::numeric_limits<double>::quiet_NaN());
        virtual ~ObsCoordSys() {};
        virtual boost::shared_ptr<CoordSys> clone() const;
        virtual boost::shared_ptr<CoordSys> clone(double date) const;
        virtual void setDate(double date) { setDate(date, false); };
        virtual Coord fromFK5J2000(Coord const &coord, Site const &site) const;
        virtual Coord toFK5J2000(Coord const &coord, Site const &site) const;
        virtual Coord fromAppTopo(Coord const &coord, Site const &site) const;
        virtual Coord toAppTopo(Coord const &coord, Site const &site) const;

        /**
        Set the date and optionally freeze cached apparent geocentric data
        
        @param[in] date: date as TAI (MJD, seconds)
        @param[in] freezeCache: if true, do not update cached apparent geocentric data
        
        The reason for freezeCache is to support computing velocity by computing position at two nearby times.
        Set freezeCache true for the second computation to avoid an unexpected cache update causing
        a mis-computation of velocity.
        
        @note The standard setDate, with no freezeCache, argument does not freeze the cache.
        
        @raise std::runtime_error if freezeCache true and delta date > AppGeoCoordSys's default MaxAge * 2;
        this is only intended to catch gross errors
        */
        virtual void setDate(double date, bool freezeCache);

    private:
        AppTopoCoordSys _appTopoCoordSys;
    };
    
    boost::shared_ptr<CoordSys> makeCoordSys(std::string const &name, double date);
}
