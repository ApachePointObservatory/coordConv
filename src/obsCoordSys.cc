#include <stdexcept>
#include <sstream>
#include "coordConv/mathUtils.h"
#include "coordConv/coordSys.h"

namespace coordConv {

    ObsCoordSys::ObsCoordSys(double date)
    :
        ApparentCoordSys("obs", date),
        _appTopoCoordSys()
    {
        setDate(date);
    };
    
    void ObsCoordSys::setDate(double date, bool freezeCache) {
        this->_date = date;
        _appTopoCoordSys.setDate(date, freezeCache);
    }
    
    CoordSys::Ptr ObsCoordSys::clone() const {
        return clone(getDate());
    }

    CoordSys::Ptr ObsCoordSys::clone(double date) const {
        return CoordSys::Ptr(new ObsCoordSys(date));
    };

    Coord ObsCoordSys::fromFK5J2000(Coord const &coord, Site const &site) const {
        Coord appTopoCoord = _appTopoCoordSys.fromFK5J2000(coord, site);
        return fromAppTopo(appTopoCoord, site);
    }

    Coord ObsCoordSys::toFK5J2000(Coord const &coord, Site const &site) const {
        Coord appTopoCoord = toAppTopo(coord, site);
        return _appTopoCoordSys.toFK5J2000(appTopoCoord, site);
    }

    Coord ObsCoordSys::fromAppTopo(Coord const &coord, Site const &site) const {
        Eigen::Vector3d appTopoPos = coord.getVecPos();

        // For zdu > ZDu_Max the correction is computed at ZDu_Max.
        // This is unphysical, but allows working with arbitrary positions.
        // The model used (at this writing) is not much good beyond 83 degrees
        // and going beyond ~87 requires more iterations to give reversibility
        double const ZDu_Max = 85.0;

        // convert inputs to easy-to-read variables
        double const xu = appTopoPos(0);
        double const yu = appTopoPos(1);
        double const zu = appTopoPos(2);

        // useful quantities
        double const rxymag = hypot(xu, yu);
        double const rxysq = rxymag * rxymag;

        Eigen::Vector3d obsPos;
        if (rxysq * std::numeric_limits<double>::epsilon() <= std::numeric_limits<double>::min()) {
            if ((rxysq + (zu * zu)) * std::numeric_limits<double>::epsilon() <= std::numeric_limits<double>::min()) {
                // |R| is too small to use -- probably a bug in the calling software
                throw std::runtime_error("appTopoPos too short");
            } else {
                // at zenith; set output = input
                obsPos = appTopoPos;
            }
        } else {
            // unrefracted zenith distance
            double zdu = atan2d(rxymag, zu);

            // Compute the refraction correction using an iterative approximation (see details).
            // Compute it at the unrefracted zenith distance, unless that ZD is too large,
            // in which case compute the correction at the max unrefracted zenith distance.
            double zdr_u = 0.0;
            double zdu_iter = zdu;
            if (zdu_iter > ZDu_Max) {
               zdu_iter = ZDu_Max;
            }
            for (int iter = 0; iter < 2; ++iter) {
               double zdr_iter = zdu_iter + zdr_u;
               double cosZD = cosd(zdr_iter);
               double tanZD = tand(zdr_iter);
               zdr_u = zdr_u - ((zdr_u + (site.refCoA * tanZD) + (site.refCoB * tanZD * tanZD * tanZD)) /
                    (1.0 + (RadPerDeg * (site.refCoA + (3.0 * site.refCoB * tanZD * tanZD)) / (cosZD * cosZD))));
            }

            // compute refracted position as a cartesian vector
            double zdr = zdu + zdr_u;
            obsPos <<
                xu,
                yu,
                rxymag * tand(90.0 - zdr);
        }
        return Coord(obsPos);
    }
    
    Coord ObsCoordSys::toAppTopo(Coord const &coord, Site const &site) const {

        // For zdu > ZDu_Max the correction is computed at ZDu_Max.
        // This is unphysical, but allows working with arbitrary positions.
        // The model used (at this writing) is not much good beyond 83 degrees
        // and going beyond ~87 requires more iterations to give reversibility
        const double ZDu_Max = 85.0;
        
        Eigen::Vector3d obsPos = coord.getVecPos();

        // convert inputs to easy-to-read variables
        double const xr = obsPos(0);
        double const yr = obsPos(1);
        double const zr = obsPos(2);

        // useful quantities
        double const rxymag = hypot(xr, yr);
        double const rxysq = rxymag * rxymag;

        // test input vector
        if (rxysq * std::numeric_limits<double>::epsilon() <= std::numeric_limits<double>::min()) {
           if ((rxysq + (zr * zr)) * std::numeric_limits<double>::epsilon() <= std::numeric_limits<double>::min()) {
                // |R| is too small to use -- probably a bug in the calling software
                throw std::runtime_error("obsPos too short");
            } else {
                // at zenith; return obsPos
                return Coord(obsPos);
            }
        }

        double zdr = atan2d(rxymag, zr); // refracted zenith distance
        double zdu; // unrefracted zenith distance

        // Compute the refraction correction. Compute it at the refracted zenith distance,
        // unless that ZD is too large, in which case compute the correction at the
        // maximum UNrefracted ZD (this provides reversibility with refract).
        bool tooLow = false;
        if (zdr > ZDu_Max) {
           // zdr < zdu, so we're certainly past the limit
           // don't even bother to try computing the standard correction
           tooLow = true;
        } else {
            double tanZD = tand(zdr);
            zdu = zdr + (site.refCoA * tanZD) + (site.refCoB * tanZD * tanZD * tanZD);
            if (zdu > ZDu_Max) {
                tooLow = true;
            }
        }

        if (tooLow) {
            // compute correction at zdu = ZDu_Max and use that instead
            // (iteration is required because we want the correction at a known zdu, not at a known zdr)
            double ZDr_u = 0.0;
            double ZDu_iter = ZDu_Max;
            for (int iter = 0; iter < 2; ++iter) {
                double ZDr_iter = ZDu_iter + ZDr_u;
                double cosZD = cosd(ZDr_iter);
                double tanZD = tand(ZDr_iter);
                ZDr_u = ZDr_u - ((ZDr_u + (site.refCoA * tanZD) + (site.refCoB * tanZD * tanZD * tanZD)) /
                    (1.0 + (RadPerDeg * (site.refCoA + (3.0 * site.refCoB * tanZD * tanZD)) / (cosZD * cosZD))));
            }
            zdu = zdr - ZDr_u;
        }

        // compute unrefracted position as a cartesian vector
        Eigen::Vector3d appTopoPos;
        appTopoPos <<
            xr,
            yr,
            rxymag * tand(90.0 - zdu);
        return Coord(appTopoPos);
    }

}

