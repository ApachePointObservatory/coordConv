#include <stdexcept>
#include <sstream>
#include "slalib.h"
#include "coordConv/mathUtils.h"
#include "coordConv/coordSys.h"

namespace {
    // WARNING: this data is all transposed from how it will be in the Eigen matrix
    static double const _FromICRSPPCArr[9] = {
        +9.999256794999100e-01, -1.118148278880500e-02, -4.859004008828000e-03,
        +1.118148284078200e-02, +9.999374848980310e-01, -2.715574495700000e-05,
        +4.859003889183000e-03, -2.717714350100000e-05, +9.999881946018790e-01
    };

    static double const _FromICRSPVCArr[9] = {
        -4.999649348699710e+01, +5.590898123577490e-01, +2.429267666920290e-01,
        -5.590898123577490e-01, -4.999708365186070e+01, +1.358274375617750e-03,
        -2.429267666920290e-01, +1.358171243214630e-03, -4.999961940950930e+01
    };

    static double const _FromICRSVPCArr[9] = {
        -2.625948691330010e-11, -1.153677408253260e-08, +2.114845700512120e-08,
        +1.153432497186020e-08, -1.289946907269570e-10, -4.139139814877260e-10,
        -2.114281972092300e-08, +5.943248704376379e-10, +1.027351973197720e-10
    };

    static double const _FromICRSVVCArr[9] = {
        +9.999043220431060e-01, -1.118145160106900e-02, -4.858519608686000e-03,
        +1.118145160896800e-02, +9.999161253401070e-01, -2.716261435500000e-05,
        +4.858519590501000e-03, -2.716586669100000e-05, +9.999668381314190e-01
    };

    static double const _ToICRSPPCArr[9] = {
        +9.999256781869020e-01, +1.118205957176600e-02, +4.857946721186000e-03,
        -1.118205964224700e-02, +9.999374784481320e-01, -2.714742649800000e-05,
        -4.857946558960000e-03, -2.717644118500000e-05, +9.999881997387700e-01
    };

    static double const _ToICRSPVCArr[9] = {
        +4.999756134052550e+01, +5.591143166167311e-01, +2.429089660392500e-01,
        -5.591143166167311e-01, +4.999815140225670e+01, -1.357552448795890e-03,
        -2.429089454127690e-01, -1.358748784672120e-03, +5.000068746930250e+01
    };

    static double const _ToICRSVPCArr[9] = {
        -2.626004779032070e-11, +1.153457133383040e-08, -2.114327131099750e-08,
        -1.153702049680800e-08, -1.289974459280040e-10, +5.943375646390270e-10,
        +2.114890871560100e-08, -4.139228222879730e-10, +1.027373916437010e-10
    };

    static double const _ToICRSVVCArr[9] = {
        +9.999470351546140e-01, +1.118250600724200e-02, +4.857669948650000e-03,
        -1.118250612180500e-02, +9.999588338188330e-01, -2.713730953900000e-05,
        -4.857669684959000e-03, -2.718447137100000e-05, +1.000009560363560e+00
    };

    static const Eigen::Map<const Eigen::Matrix3d> FromICRSPP(_FromICRSPPCArr);
    static const Eigen::Map<const Eigen::Matrix3d> FromICRSPV(_FromICRSPVCArr);
    static const Eigen::Map<const Eigen::Matrix3d> FromICRSVP(_FromICRSVPCArr);
    static const Eigen::Map<const Eigen::Matrix3d> FromICRSVV(_FromICRSVVCArr);

    static const Eigen::Map<const Eigen::Matrix3d> ToICRSPP(_ToICRSPPCArr);
    static const Eigen::Map<const Eigen::Matrix3d> ToICRSPV(_ToICRSPVCArr);
    static const Eigen::Map<const Eigen::Matrix3d> ToICRSVP(_ToICRSVPCArr);
    static const Eigen::Map<const Eigen::Matrix3d> ToICRSVV(_ToICRSVVCArr);
}

namespace coordConv {

    FK4CoordSys::FK4CoordSys(double date)
    :
        CoordSys("fk4", date),
        _eTerms(),
        _From1950PrecMat(),
        _To1950PrecMat()
    {
        setDate(date);
    };

    boost::shared_ptr<CoordSys> FK4CoordSys::clone() const {
        return clone(getDate());
    }

    boost::shared_ptr<CoordSys> FK4CoordSys::clone(double date) const {
        return boost::shared_ptr<CoordSys>(new FK4CoordSys(date));
    };
    
    void FK4CoordSys::setDate(double date) {
        this->_date = date;
        if (isFinite(date)) {
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

    Coord FK4CoordSys::fromICRS(Coord const &coord, Site const &site) const {
        Eigen::Vector3d icrsPos = coord.getVecPos();
        Eigen::Vector3d icrsVel = coord.getVecVel();

        // convert position and velocity from J2000.0 to B1950
        Eigen::Vector3d b1950Pos = (FromICRSPP * icrsPos) + (FromICRSPV * icrsVel);
        Eigen::Vector3d b1950Vel = (FromICRSVP * icrsPos) + (FromICRSVV * icrsVel);

        // correct position for velocity (PM and rad. vel.) from 1950 to date
        Eigen::Vector3d tempPos = b1950Pos + ((this->_date - 1950.0) * b1950Vel);

        // precess position and velocity from 1950 to date
        Eigen::Vector3d meanToPos = _From1950PrecMat * tempPos;
        Eigen::Vector3d fk4Vel    = _From1950PrecMat * b1950Vel;

        // add e-terms to mean position, iterating thrice (should be plenty!)
        // to get mean catalog place. As a minor approximation,
        // we don't bother to add variation in e-terms to the velocity.
        Eigen::Vector3d fk4Pos = meanToPos;
        for (int iter = 0; iter < 3; ++iter) {
            double magPos = fk4Pos.norm();
            fk4Pos = meanToPos + (magPos * _eTerms);
        }
        
        return Coord(fk4Pos, fk4Vel);
    };
     
    Coord FK4CoordSys::toICRS(Coord const &coord, Site const &site) const {
        // use the excellent approximation that ICRS = FK4 J2000
        Eigen::Vector3d fk4Pos = coord.getVecPos();
        Eigen::Vector3d fk4Vel = coord.getVecVel();

        // subtract e-terms from position
        double magPos = fk4Pos.norm();
        Eigen::Vector3d meanFK4Pos = fk4Pos - (magPos * _eTerms);
        
        if ((fk4Vel.array() == 0.0).all()) {
            // object is fixed on the sky; handle FK4 fictitious proper motion
            
            // precess position and velocity to B1950
            Eigen::Vector3d fk41950Pos = _To1950PrecMat * meanFK4Pos;

            // convert position to J2000.0 and compute fictitious velocity
            Eigen::Vector3d tempPos = ToICRSPP * fk41950Pos;
            Eigen::Vector3d ficVel  = ToICRSVP * fk41950Pos;
        
            // subtract fictitious velocity over the period this->_date to J2000
            double period = 2000.0 - slaEpj(slaEpb2d(this->_date));
            Eigen::Vector3d icrsPos = tempPos - (ficVel * period);
            
            return Coord(icrsPos);
        
        } else {
            // proper motion specified

            // correct position for velocity (proper motion and radial velocity) to B1950
            Eigen::Vector3d corrPos = meanFK4Pos + (1950.0 - this->_date) * fk4Vel;

            // precess position and velocity to B1950
            Eigen::Vector3d fk41950Pos = _To1950PrecMat * corrPos;
            Eigen::Vector3d fk41950Vel = _To1950PrecMat * fk4Vel;

            // convert position and velocity to J2000.0
            Eigen::Vector3d icrsPos = (ToICRSPP * fk41950Pos) + (ToICRSPV * fk41950Vel);
            Eigen::Vector3d icrsVel = (ToICRSVP * fk41950Pos) + (ToICRSVV * fk41950Vel);

            return Coord(icrsPos, icrsVel);
        }
    };
   
}
