#include <sstream>
#include "slalib.h"
#include "coordConv/physConst.h"
#include "coordConv/site.h"

namespace coordConv {

    void Site::setPoleWander(double x, double y) {
        poleX = x;
        poleY = y;
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

    std::string Site::__repr__() const {
        std::ostringstream os;
        os << "site.meanLong=" << meanLong << "\n";
        os << "site.meanLat=" << meanLat << "\n";
        os << "site.elev=" << elev << "\n";
        os << "site.poleX=" << poleX << "\n";
        os << "site.poleY=" << poleY << "\n";
        os << "site.ut1_tai=" << ut1_tai << "\n";
        os << "site.utc_tai=" << utc_tai << "\n";
        os << "site.corrLong=" << corrLong << "\n";
        os << "site.corrLat=" << corrLat << "\n";
        os << "site.wavelen=" << wavelen << "\n";
        os << "site.refCoA=" << refCoA << "\n";
        os << "site.refCoB=" << refCoB << "\n";
        os << "site.azCorr=" << azCorr << "\n";
        os << "site.diurAbMag=" << diurAbMag << "\n";
        os << "site.pos=" << pos(0) << ", " << pos(1) << ", " << pos(2);
        return os.str();
    }

    std::ostream &operator<<(std::ostream &os, Site const &site) {
        // use overloaded __repr__
        os << site.__repr__();
        return os;
    }

}
