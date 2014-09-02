#include "coordConv/time.h"
#include "coordConv/azAltFromHADec.h"
#include "coordConv/mathUtils.h"
#include "coordConv/appTopoFromAppGeo.h"

namespace coordConv {

    Coord appTopoFromAppGeo(Coord const &coord, Site const &site, double tai) {
        Eigen::Vector3d appGeoPos = coord.getVecPos();

        double const last = lastFromTAI(tai, site);
        double const sinLAST = sind(last);
        double const cosLAST = cosd(last);

        // rotate position and offset from RA/Dec to (-HA)/Dec
        Eigen::Vector3d posA;
        posA <<
            + (cosLAST * appGeoPos(0)) + (sinLAST * appGeoPos(1)),
            - (sinLAST * appGeoPos(0)) + (cosLAST * appGeoPos(1)),
               appGeoPos(2);

        // correct position for diurnal parallax
        Eigen::Vector3d posB = posA - site.pos;

        // correct position for diurnal aberration; follows Pat Wallace's slaAOPQK
        Eigen::Vector3d posC;
        double bMag = posB.norm();
        double diurAbScaleCorr = 1.0 - (site.diurAbMag * (posB(1) / bMag));
        posC(0) =  posB(0) * diurAbScaleCorr;
        posC(1) = (posB(1) + (site.diurAbMag * bMag)) * diurAbScaleCorr;
        posC(2) =  posB(2) * diurAbScaleCorr;

        // rotate position from -HA/Dec to alt/az; use latitude corrected for pole wander
        Eigen::Vector3d appTopoPos;
        azAltFromHADec(appTopoPos, posC, site.corrLat);
        
        return Coord(appTopoPos);
    }

}
