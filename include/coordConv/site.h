#pragma once

#include "Eigen/Dense"

namespace coordConv {

    /**
    Site information
    */
    class Site {
    public:
    
        double meanLong;    ///< site longitude, ignoring pole wander (deg, positive eastward)
        double meanLat;     ///< site latitude, ignoring pole wander (deg)
        double elev;        ///< geodetic elevation (meters above reference spheroid)
        double poleX;       ///< pole wander, X (deg)
        double poleY;       ///< pole wander, Y (deg)
        double ut1_tai;     ///< UT1-TAI (seconds) at date of coordinate conversion
        double utc_tai;     ///< UTC-TAI (seconds) at date of coordinate conversion
            ///< not used by coordConv, but often useful and easily computed at the same time as ut1_tai
        double corrLong;    ///< longitude corrected for pole wander (deg)
        double corrLat;     ///< latitude corrected for pole wander (deg)
        double wavelen;     ///< wavelength for which to compute refraction coefficients (Angstroms)
        double refCoA, refCoB;  ///< A, B refraction coefficients (deg), where:
            ///< zdSpace = refCoA tan zdEarth + refCoB tan^3 zdEarth
            ///< zdEarth is the zenith distance of an object as observed through the atmosphere (radians)
            ///< zdSpace is the zenith distance of an object if there was no refraction (radians)
        double azCorr;      ///< azimuth correction (terrestrial-celestial, deg)
        double diurAbMag;   ///< magnitude of diurnal aberration vector:
            /// speed of rotation of observatory / speed of light (radians/au)
        Eigen::Vector3d pos;    ///< cartesian position of observatory (au)
        
        /**
        Construct a new Site.
        
        Initialize with no refraction correction and no polar wander.
        
        @param[in] meanLong;    ///< site longitude, ignoring pole wander (deg, positive eastward)
        @param[in] meanLat;     ///< site latitude, ignoring pole wander (deg)
        @param[in] elev;        ///< site elevation (meters above sea level)
        
        @throw std::range_error if meanLat is not in range [-90, 90]
        */
        Site(double meanLong, double meanLat, double elev);
        
        /**
        Set current pole wander, based on the USNO earth orientation bulletin
        
        @param[in] x: x correction (deg)
        @param[in] y: y correction (deg)
        
        Updates the following fields:
        - poleX
        - poleY
        - corrLong
        - corrLat
        - azCorr
        - diurAbMag
        - pos
        */
        void setPoleWander(double x, double y);
        
        /**
        Print a string representation
        */
        std::string __repr__() const;
    };

    std::ostream &operator<<(std::ostream &out, Site const &site);

}
