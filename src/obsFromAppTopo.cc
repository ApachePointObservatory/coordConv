#include <stdexcept>
#include "coordConv/mathUtils.h"
#include "coordConv/obsFromAppTopo.h"

namespace coordConv {

    Coord obsFromAppTopo(Coord const &appTopoCoord, Site const &site) {
        Eigen::Vector3d appTopoPos = appTopoCoord.getVecPos();

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

            // Compute the refraction correction using an iterative approximation;
            // based on tests 2 iterations is plenty, but do one more for paranoia's sake.
            // Compute it at the unrefracted zenith distance, unless that ZD is too large,
            // in which case compute the correction at the max unrefracted zenith distance.
            double zdr_u = 0.0;
            double zdu_iter = zdu;
            if (zdu_iter > ZDu_Max) {
               zdu_iter = ZDu_Max;
            }
            for (int iter = 0; iter < 3; ++iter) {
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

}

