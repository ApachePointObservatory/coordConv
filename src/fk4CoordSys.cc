#include <cmath>
#include <stdexcept>
#include <sstream>
#include "slalib.h"
#include "coordConv/mathUtils.h"
#include "coordConv/physConst.h"
#include "coordConv/coordSys.h"

namespace {
    // WARNING: this data is all transposed from how it will be in the Eigen matrix
    static double const _fromFK5J2000PPCArr[9] = {
        +9.999256794999100e-01, -1.118148278880500e-02, -4.859004008828000e-03,
        +1.118148284078200e-02, +9.999374848980310e-01, -2.715574495700000e-05,
        +4.859003889183000e-03, -2.717714350100000e-05, +9.999881946018790e-01
    };

    static double const _fromFK5J2000PVCArr[9] = {
        -4.999649348699710e+01, +5.590898123577490e-01, +2.429267666920290e-01,
        -5.590898123577490e-01, -4.999708365186070e+01, +1.358274375617750e-03,
        -2.429267666920290e-01, +1.358171243214630e-03, -4.999961940950930e+01
    };

    static double const _fromFK5J2000VPCArr[9] = {
        -2.625948691330010e-11, -1.153677408253260e-08, +2.114845700512120e-08,
        +1.153432497186020e-08, -1.289946907269570e-10, -4.139139814877260e-10,
        -2.114281972092300e-08, +5.943248704376379e-10, +1.027351973197720e-10
    };

    static double const _fromFK5J2000VVCArr[9] = {
        +9.999043220431060e-01, -1.118145160106900e-02, -4.858519608686000e-03,
        +1.118145160896800e-02, +9.999161253401070e-01, -2.716261435500000e-05,
        +4.858519590501000e-03, -2.716586669100000e-05, +9.999668381314190e-01
    };

    static double const _toFK5J2000PPCArr[9] = {
        +9.999256781869020e-01, +1.118205957176600e-02, +4.857946721186000e-03,
        -1.118205964224700e-02, +9.999374784481320e-01, -2.714742649800000e-05,
        -4.857946558960000e-03, -2.717644118500000e-05, +9.999881997387700e-01
    };

    static double const _toFK5J2000PVCArr[9] = {
        +4.999756134052550e+01, +5.591143166167311e-01, +2.429089660392500e-01,
        -5.591143166167311e-01, +4.999815140225670e+01, -1.357552448795890e-03,
        -2.429089454127690e-01, -1.358748784672120e-03, +5.000068746930250e+01
    };

    static double const _toFK5J2000VPCArr[9] = {
        -2.626004779032070e-11, +1.153457133383040e-08, -2.114327131099750e-08,
        -1.153702049680800e-08, -1.289974459280040e-10, +5.943375646390270e-10,
        +2.114890871560100e-08, -4.139228222879730e-10, +1.027373916437010e-10
    };

    static double const _toFK5J2000VVCArr[9] = {
        +9.999470351546140e-01, +1.118250600724200e-02, +4.857669948650000e-03,
        -1.118250612180500e-02, +9.999588338188330e-01, -2.713730953900000e-05,
        -4.857669684959000e-03, -2.718447137100000e-05, +1.000009560363560e+00
    };

    static const Eigen::Map<const Eigen::Matrix3d> fromFK5J2000PP(_fromFK5J2000PPCArr);
    static const Eigen::Map<const Eigen::Matrix3d> fromFK5J2000PV(_fromFK5J2000PVCArr);
    static const Eigen::Map<const Eigen::Matrix3d> fromFK5J2000VP(_fromFK5J2000VPCArr);
    static const Eigen::Map<const Eigen::Matrix3d> fromFK5J2000VV(_fromFK5J2000VVCArr);

    static const Eigen::Map<const Eigen::Matrix3d> toFK5J2000PP(_toFK5J2000PPCArr);
    static const Eigen::Map<const Eigen::Matrix3d> toFK5J2000PV(_toFK5J2000PVCArr);
    static const Eigen::Map<const Eigen::Matrix3d> toFK5J2000VP(_toFK5J2000VPCArr);
    static const Eigen::Map<const Eigen::Matrix3d> toFK5J2000VV(_toFK5J2000VVCArr);
}

namespace coordConv {

    FK4CoordSys::FK4CoordSys(double date)
    :
        MeanCoordSys("fk4", date),
        _eTerms(),
        _From1950PrecMat(),
        _To1950PrecMat()
    {
        setDate(date);
    };

    CoordSys::Ptr FK4CoordSys::clone() const {
        return clone(getDate());
    }

    CoordSys::Ptr FK4CoordSys::clone(double date) const {
        return CoordSys::Ptr(new FK4CoordSys(date));
    };
    
    void FK4CoordSys::setDate(double date) {
        this->_date = date;
        if (std::isfinite(date)) {
            // note: slaEtrms and slaPrebn both want Besselian date
            double eTermsCArr[3], from1950PrecCArr[3][3], to1950PrecCArr[3][3];
            slaEtrms(date, eTermsCArr);
            slaPrebn(1950.0, date, from1950PrecCArr);
            slaPrebn(date, 1950.0, to1950PrecCArr);
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    _From1950PrecMat(i,j) = from1950PrecCArr[i][j];
                      _To1950PrecMat(i,j) =   to1950PrecCArr[i][j];
                }
                _eTerms(i) = eTermsCArr[i];
            }
        }
    }

    Coord FK4CoordSys::fromFK5J2000(Coord const &coord, Site const &site) const {
        Eigen::Vector3d fk5J2000Pos = coord.getVecPos();
        Eigen::Vector3d fk5J2000PM = coord.getVecPM();

        // convert position and velocity from J2000.0 to B1950
        Eigen::Vector3d b1950Pos = (fromFK5J2000PP * fk5J2000Pos) + (fromFK5J2000PV * fk5J2000PM);
        Eigen::Vector3d b1950Vel = (fromFK5J2000VP * fk5J2000Pos) + (fromFK5J2000VV * fk5J2000PM);

        // correct position for velocity (PM and rad. vel.) from 1950 to date
        Eigen::Vector3d tempPos = b1950Pos + ((this->_date - 1950.0) * b1950Vel);

        // precess position and velocity from 1950 to date
        Eigen::Vector3d meanToPos = _From1950PrecMat * tempPos;
        Eigen::Vector3d fk4PM     = _From1950PrecMat * b1950Vel;

        // add e-terms to mean position, iterating thrice (should be plenty!)
        // to get mean catalog place. As a minor approximation,
        // we don't bother to add variation in e-terms to the velocity.
        Eigen::Vector3d fk4Pos = meanToPos;
        for (int iter = 0; iter < 3; ++iter) {
            double magPos = fk4Pos.norm();
            fk4Pos = meanToPos + (magPos * _eTerms);
        }
        
        return Coord(fk4Pos, fk4PM);
    };
     
    Coord FK4CoordSys::toFK5J2000(Coord const &coord, Site const &site) const {
        // use the excellent approximation that ICRS = FK4 J2000
        Eigen::Vector3d fk4Pos = coord.getVecPos();
        Eigen::Vector3d fk4PM = coord.getVecPM();

        // subtract e-terms from position
        double magPos = fk4Pos.norm();
        Eigen::Vector3d meanFK4Pos = fk4Pos - (magPos * _eTerms);
        
        if ((fk4PM.array() == 0.0).all()) {
            // object is fixed on the sky; handle FK4 fictitious proper motion
            
            // precess position and velocity to B1950
            Eigen::Vector3d fk41950Pos = _To1950PrecMat * meanFK4Pos;

            // convert position to J2000.0 and compute fictitious velocity
            Eigen::Vector3d tempPos = toFK5J2000PP * fk41950Pos;
            Eigen::Vector3d ficVel  = toFK5J2000VP * fk41950Pos;
        
            // subtract fictitious velocity over the period this->_date to J2000
            double period = 2000.0 - slaEpj(slaEpb2d(this->_date));
            Eigen::Vector3d fk5J2000Pos = tempPos - (ficVel * period);
            
            return Coord(fk5J2000Pos);
        
        } else {
            // proper motion specified

            // correct position for velocity (proper motion and radial velocity) to B1950
            Eigen::Vector3d corrPos = meanFK4Pos + (1950.0 - this->_date) * fk4PM;

            // precess position and velocity to B1950
            Eigen::Vector3d fk41950Pos = _To1950PrecMat * corrPos;
            Eigen::Vector3d fk41950Vel = _To1950PrecMat * fk4PM;

            // convert position and velocity to J2000.0
            Eigen::Vector3d fk5J2000Pos = (toFK5J2000PP * fk41950Pos) + (toFK5J2000PV * fk41950Vel);
            Eigen::Vector3d fk5J2000PM = (toFK5J2000VP * fk41950Pos) + (toFK5J2000VV * fk41950Vel);

            return Coord(fk5J2000Pos, fk5J2000PM);
        }
    };
   
    double FK4CoordSys::dateFromTAI(double tai) const {
        return slaEpb((tai + TT_TAI) / SecPerDay); // slaEpb converts MJD days to Besselian epoch
    }

    std::string FK4CoordSys::__repr__() const {
        std::ostringstream os;
        os << "FK4CoordSys(" << getDate() << ")";
        return os.str();
    }

}
