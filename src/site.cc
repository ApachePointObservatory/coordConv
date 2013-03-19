#include "slalib.h"
#include "coordConv/physConst.h"
#include "coordConv/site.h"

namespace coordConv {

    void Site::setPoleWander(double x, double y) {
        double corrLongRad, corrLatRad, azCorrRad;
        slaPolmo(
            meanLong * RadPerDeg,
            meanLat * RadPerDeg,
            x * RadPerDeg,
            y * RadPerDeg,
            &corrLongRad,
            &corrLatRad,
            &azCorrRad
        );
        corrLong = corrLongRad / RadPerDeg;
        corrLat = corrLatRad / RadPerDeg;
        azCorr = azCorrRad / RadPerDeg;
        double polarDist;   ///< distance of observatory from Earth's axis (au)
        double zDist;   ///< distance of observatory from plane of Earth's equator (au)
        slaGeoc(corrLongRad, elev, &polarDist, &zDist);
        
        pos << polarDist, 0.0, zDist;

        const double SidRate = 2 * Pi * SiderealPerSolar / SecPerDay;   ///< sidereal rate (rad/sec)
        diurAbMag = polarDist * SidRate / VLight;
    }

}
