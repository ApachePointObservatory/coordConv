#include <stdexcept>
#include <sstream>
#include "coordConv/mathUtils.h"
#include "coordConv/physConst.h"
#include "coordConv/coordSys.h"

namespace coordConv {

    Coord::Coord(double equatAng, double polarAng, double parallax) {
        _setPosFromSph(equatAng, polarAng, parallax);
        _setCache();
        _vel.setZero();
    }

    Coord::Coord(double equatAng, double polarAng, double parallax, double equatPM, double polarPM, double radVel) {
        _setPosFromSph(equatAng, polarAng, parallax);
        _setCache();

        const double RadPerYear_per_ArcsecPerCentury = RadPerDeg / (ArcsecPerDeg * 100.0);
        const double AUPerYear_per_KmPerSec = SecPerDay * DaysPerYear / KmPerAU;

        const double sinEquat = sind(equatAng);
        const double cosEquat = cosd(equatAng);
        const double sinPolar = sind(polarAng);
        const double cosPolar = cosd(polarAng);;

        // change units of proper motion from arcsec/century to au/year
        // (multiply by distance and fix the units)
        const double pmAUPerYear1 = equatPM * _dist * RadPerYear_per_ArcsecPerCentury;
        const double pmAUPerYear2 = polarPM * _dist * RadPerYear_per_ArcsecPerCentury;

        // change units of radial velocity from km/sec to au/year
        const double radVelAUPerYear = radVel * AUPerYear_per_KmPerSec;
//        const double radVelAUPerYear = _atInfinity ? 0 : radVel * AUPerYear_per_KmPerSec;

        // compute velocity vector in au/year
        _vel <<
            - (pmAUPerYear2 * sinPolar * cosEquat) - (pmAUPerYear1 * cosPolar * sinEquat) + (radVelAUPerYear * cosPolar * cosEquat),
            - (pmAUPerYear2 * sinPolar * sinEquat) + (pmAUPerYear1 * cosPolar * cosEquat) + (radVelAUPerYear * cosPolar * sinEquat),
            + (pmAUPerYear2 * cosPolar)                                                   + (radVelAUPerYear * sinPolar);
    }

    Coord::Coord(Eigen::Vector3d &pos, Eigen::Vector3d &vel)
    :
        _pos(pos),
        _vel(vel)
    {
        _setCache();
    }

    Coord::Coord(Eigen::Vector3d &pos)
    :
        _pos(pos),
        _vel()
    {
        _setCache();
        _vel.setZero();
    }
    
    double Coord::angularSeparation(Coord const &coord) const {
        double crossMag = _pos.cross(coord.getVecPos()).norm();
        double dotProd = _pos.dot(coord.getVecPos());
        return atan2d(crossMag, dotProd);
    }

    double Coord::getParallax() const {
        return _atInfinity ? 0 : AUPerParsec / _dist;
    }

    bool Coord::getSphPos(double &equatAng, double &polarAng) const {
        double x = _pos(0);
        double y = _pos(1);
        double z = _pos(2);

        if (_atPole) {
            equatAng = 0.0;
            polarAng = (z > 0.0) ? 90.0 : -90.0;
        } else {
            double posXYMag = hypot(x, y);
            equatAng = wrapPos(atan2d(y, x));
            polarAng = atan2d(z, posXYMag);
        }
        return _atPole;
    }
    
    bool Coord::getPM(double &equatPM, double &polarPM) const {
        double const ArcsecPerCentury_per_RadPerYear = 100.0 * ArcsecPerDeg / RadPerDeg;

        double x  = _pos(0);
        double y  = _pos(1);
        double z  = _pos(2);
        double vX = _vel(0);
        double vY = _vel(1);
        double vZ = _vel(2);
    
        // now that radial velocity has been computed
        // handle the "at pole" case
        if (_atPole) {
           equatPM = 0.0;
           polarPM = 0.0;
           return _atPole;
        }

        // useful quantities
        double magPxy = hypot(x, y);
        double magPxySq  = magPxy * magPxy;
        double magPSq = _dist * _dist;

        // compute proper motion in rad per year,
        // then convert to arcsec per century;
        // the divisions are safe because:
        // - magPxySq must have some reasonable minimum value,
        //   else scFromCC would have set atPole true,
        //   and that case has already been handled above
        // - magPSq must have some reasonable minimum value,
        //   else scFromCC would have set isOK false
        equatPM = (((x * vY) - (y * vX)) / magPxySq) * ArcsecPerCentury_per_RadPerYear;
        polarPM = (((vZ * magPxy) - ((z / magPxy) * ((x * vX) + (y * vY)))) / magPSq) * ArcsecPerCentury_per_RadPerYear;
        return _atPole;
    }
    
    double Coord::getRadVel() const {
        // compute radial velocity in (au/year) and convert to (km/s)
        double const KMPerSec_per_AUPerYear = KmPerAU / (DaysPerYear * SecPerDay);
        return (_pos / _dist).dot(_vel) * KMPerSec_per_AUPerYear;
//        return _atInfinity ? 0 : _pos.dot(_vel) * KMPerSec_per_AUPerYear / _dist;
    }

    void Coord::_setPosFromSph(double equatAng, double polarAng, double parallax) {
        if ((polarAng < -90.0) || (polarAng > 90.0)) {
            std::ostringstream os;
            os << "polarAng = " << polarAng << " not in range [-90, 90]";
            throw std::runtime_error(os.str());
        }
        
        double dist = AUPerParsec / std::max(parallax, MinParallax);

        _pos <<
            dist * cosd(polarAng) * cosd(equatAng),
            dist * cosd(polarAng) * sind(equatAng),
            dist * sind(polarAng);
    }
    
    void Coord::_setCache() {
        _dist = _pos.norm();

        // make sure |_pos| is large enough to compute with
        if (_dist * _dist < std::numeric_limits<double>::min()) {
            std::ostringstream os;
            os << "Magnitude of _pos = (" << _pos(0) << ", " << _pos(1) << ", " << _pos(2) << ") too small";
            throw std::runtime_error(os.str());
        }
        _atInfinity = _dist > 0.9 * AUPerParsec / MinParallax;
        // this test for atPole is based on reliably round-tripping equatPM to 6 digits
        double xyFracMag = hypot(_pos(0), _pos(1)) / _dist;
        _atPole = (xyFracMag * xyFracMag < std::numeric_limits<double>::epsilon());
    }

}
