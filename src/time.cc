#include "slalib.h"
#include "coordConv/mathUtils.h"
#include "coordConv/physConst.h"
#include "coordConv/time.h"

namespace {
    inline double ttDaysFromTAI(double tai) {
        return (tai + coordConv::TT_TAI) / coordConv::SecPerDay;
    }

    inline double taiFromTTDays(double ttDays) {
        return (ttDays * coordConv::SecPerDay) - coordConv::TT_TAI;
    }
}

namespace coordConv {

    double lastFromTAI(double tai, Site const &site) {
        // compute Greenwich mean sidereal time, in degrees
        double ut1Days = (tai + site.ut1_tai) / SecPerDay;
        double gmst = slaGmst(ut1Days) / RadPerDeg;

        // compute apparent - mean sidereal time, in degrees
        double ttDays = (tai + TT_TAI) / SecPerDay;
        double appMinusMean = slaEqeqx(ttDays) / RadPerDeg;

        // compute local apparent sideral time, in degrees, in range [0, 360)
        return wrapPos(gmst + site.corrLong + appMinusMean);
    }
    
    double julianEpochFromTAI(double tai) {
        return 2000.0 + ((ttDaysFromTAI(tai) - MJDJ2000) / DaysPerYear);
    }
    
    double taiFromJulianEpoch(double julianEpoch) {
        return taiFromTTDays(MJDJ2000 + ((julianEpoch - 2000.0) * DaysPerYear));
    }

    double besselianEpochFromTAI(double tai) {
        return 1900.0 + (ttDaysFromTAI(tai) - 15019.81352 ) / 365.242198781;
    }

    double taiFromBesselianEpoch(double date) {
        return taiFromTTDays(15019.81352 + (date - 1900.0) * 365.242198781);
    }

}
