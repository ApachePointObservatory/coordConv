#include "coordConv/time.h"
#include "coordConv/haDecFromAzAlt.h"
#include "coordConv/mathUtils.h"
#include "coordConv/appGeoFromAppTopo.h"

namespace coordConv {

    Coord appGeoFromAppTopo(Coord const &coord, Site const &site, double tai) {

        double const last = lastFromTAI(tai, site);
        double const sinLAST = sind(last);
        double const cosLAST = cosd(last);
        
        Eigen::Vector3d appTopoPos = coord.getVecPos();

        // rotate position from alt/az to -HA/Dec; use latitude corrected for pole wander
        Eigen::Vector3d posC;
        haDecFromAzAlt(posC, appTopoPos, site.corrLat);

        // remove correction for diurnal aberration
        // following Pat Wallace's slaOAPQK, the same equation is used
        // as applying the correction, but the sign of diurAbMag is reversed
        Eigen::Vector3d posB;
        double bMag = posC.norm();
        double diurAbScaleCorr = 1.0 + (site.diurAbMag * (posC(1) / bMag));
        posB(0) =  posC(0) * diurAbScaleCorr;
        posB(1) = (posC(1) - (site.diurAbMag * bMag)) * diurAbScaleCorr;
        posB(2) =  posC(2) * diurAbScaleCorr;

        // correct position for diurnal parallax (needed for planets, not stars)
        Eigen::Vector3d posA = posB + site.pos;

        // rotate position from -HA/Dec to RA/Dec
        Eigen::Vector3d appGeoPos;
        appGeoPos(0) = + (cosLAST * posA(0)) - (sinLAST * posA(1));
        appGeoPos(1) = + (sinLAST * posA(0)) + (cosLAST * posA(1));
        appGeoPos(2) =    posA(2);
        return Coord(appGeoPos);
    }

}
