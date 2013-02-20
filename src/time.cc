#import "slalib.h"
#import "coordConv/mathUtils.h"
#import "coordConv/physConst.h"
#import "coordConv/time.h"

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
    
    double julianEpochFromMJDSec(double mjdSec) {
        double mjdDays = mjdSec / SecPerDay;
        return 2000.0 + ((mjdDays - 51544.5) / DaysPerYear);
    }
    
    double mjdSecFromJulianEpoch(double julianEpoch) {
        double mjdDays = 51544.5 + ((julianEpoch - 2000) * DaysPerYear);
        return mjdDays * SecPerDay;
    }
        
}
