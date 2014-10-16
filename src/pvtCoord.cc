#include <cmath>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <vector>
#include <boost/tr1/array.hpp> // so array works with old and new compilers

#include "coordConv/mathUtils.h"
#include "coordConv/pvtCoord.h"

static const double DeltaT = 0.01;

namespace coordConv {

    PVTCoord::PVTCoord(Coord const &coord, Eigen::Vector3d const &vel, double tai) :
        _coord(coord),
        _vel(vel),
        _tai(tai)
    { };
    
    PVTCoord::PVTCoord(Coord const &coord0, Coord const &coord1, double tai, double deltaT) {
        _setFromCoordPair(coord0, coord1, tai, deltaT);
    }

    PVTCoord::PVTCoord(PVT const &equatPVT, PVT const &polarPVT, PVT const &distPVT) {
        if (equatPVT.t != polarPVT.t) {
            std::ostringstream os;
            os << "equatPVT.t = " << equatPVT.t << " != " << polarPVT.t << " = polarPVT.t";
            throw std::runtime_error(os.str());
        }
        if (distPVT.isfinite() and (distPVT.t != equatPVT.t)) {
            std::ostringstream os;
            os << "distPVT is finite and distPVT.t = " << distPVT.t << " != " << equatPVT.t << " = equatPVT.t";
            throw std::runtime_error(os.str());
        }
        double const tai = equatPVT.t;
        std::vector<Coord> coordArr;
        for (int i = 0; i < 2; ++i) {
            double evalTAI = tai + (i * DeltaT);
            double parallax = 0;
            if (distPVT.isfinite()) {
                parallax = parallaxFromDistance(distPVT.getPos(evalTAI));
            }
            coordArr.push_back(Coord(equatPVT.getPos(evalTAI), polarPVT.getPos(evalTAI), parallax));
        }
        _setFromCoordPair(coordArr[0], coordArr[1], tai, DeltaT);
    }

    PVTCoord::PVTCoord(PVT const &equatPVT, PVT const &polarPVT, PVT const &distPVT, double equatPM, double polarPM, double radVel) {
        if (equatPVT.t != polarPVT.t) {
            std::ostringstream os;
            os << "equatPVT.t = " << equatPVT.t << " != " << polarPVT.t << " = polarPVT.t";
            throw std::runtime_error(os.str());
        }
        if (distPVT.isfinite() and (distPVT.t != equatPVT.t)) {
            std::ostringstream os;
            os << "distPVT is finite and distPVT.t = " << distPVT.t << " != " << equatPVT.t << " = equatPVT.t";
            throw std::runtime_error(os.str());
        }
        std::vector<Coord> coordArr;
        double const tai = equatPVT.t;
        for (int i = 0; i < 2; ++i) {
            double evalTAI = tai + (i * DeltaT);
            double parallax = 0;
            if (distPVT.isfinite()) {
                parallax = parallaxFromDistance(distPVT.getPos(evalTAI));
            }
            coordArr.push_back(Coord(equatPVT.getPos(evalTAI), polarPVT.getPos(evalTAI), parallax, equatPM, polarPM, radVel));
        }
        _setFromCoordPair(coordArr[0], coordArr[1], tai, DeltaT);
    }

    PVTCoord::PVTCoord() :
        _coord(Coord()),
        _vel(),
        _tai(DoubleNaN)
    { };

    PVTCoord PVTCoord::copy(double tai) const {
        return PVTCoord(getCoord(tai), _vel, tai);
    }
    
    Coord PVTCoord::getCoord(double tai) const {
        if (tai == _tai) {
            return _coord;
        }
        Eigen::Vector3d pos = _coord.getVecPos() + (_vel * (tai - _tai));
        return Coord(pos, _coord.getVecPM());
    }    
    
    bool PVTCoord::getSphPVT(PVT &equatPVT, PVT &polarPVT) const {
        double equatPos[2], polarPos[2];
        bool atPole = false;
        for (int i = 0; i < 2; ++i) {
            double evalDate = _tai + (i * DeltaT);
            Coord coord = getCoord(evalDate);
            atPole |= coord.getSphPos(equatPos[i], polarPos[i]);
        }
        equatPVT.setFromPair(equatPos, _tai, DeltaT, true);
        polarPVT.setFromPair(polarPos, _tai, DeltaT, false);
        return atPole;
    }

    PVT PVTCoord::getDistance() const {
        double dist[2];
        for (int i = 0; i < 2; ++i) {
            double evalTAI = _tai + (i * DeltaT);
            dist[i] = getCoord(evalTAI).getDistance();
        }
        PVT distPVT;
        distPVT.setFromPair(dist, _tai, DeltaT, false);
        return distPVT;
    }

    bool PVTCoord::isfinite() const {
        return _coord.isfinite()
            && std::isfinite(_vel(0)) && std::isfinite(_vel(1)) && std::isfinite(_vel(2))
            && std::isfinite(_tai);
    }

    PVT PVTCoord::angularSeparation(PVTCoord const &pvtCoord) const {
        double posArr[2];
        for (int i = 0; i < 2; ++i) {
            double evalTAI = _tai + (i * DeltaT);
            Coord thisCoord = getCoord(evalTAI);
            Coord otherCoord = pvtCoord.getCoord(evalTAI);
            posArr[i] = thisCoord.angularSeparation(otherCoord);
        }
        PVT res = PVT();
        res.setFromPair(posArr, _tai, DeltaT, false);
        return res;
    }

    PVT PVTCoord::orientationTo(PVTCoord const &pvtCoord) const {
        double posArr[2];
        for (int i = 0; i < 2; ++i) {
            double evalTAI = _tai + (i * DeltaT);
            Coord thisCoord = getCoord(evalTAI);
            Coord otherCoord = pvtCoord.getCoord(evalTAI);
            posArr[i] = thisCoord.orientationTo(otherCoord);
        }
        PVT res = PVT();
        // if the orientation is only finite at one of the two times
        // then the distance is 0 at the other time and the orientation is fixed
        // (since the relative velocity vector goes through this PVT)
        if (std::isfinite(posArr[0]) && std::isfinite(posArr[1])) {
            res.setFromPair(posArr, _tai, DeltaT, true);
        } else if (std::isfinite(posArr[0])) {
            res.pos = posArr[0];
            res.vel = 0;
            res.t = _tai;
        } else if (std::isfinite(posArr[1])) {
            res.pos = posArr[1];
            res.vel = 0;
            res.t = _tai;
        } else {
            res.pos = DoubleNaN;
            res.vel = DoubleNaN;
            res.t = _tai;
        }
        return res;
    }

    PVTCoord PVTCoord::offset(PVT &toOrient, PVT const &fromOrient, PVT const &dist) const {
        std::vector<Coord> coordArr;
        double toOrientArr[2];
        for (int i = 0; i < 2; ++i) {
            double evalTAI = _tai + (i * DeltaT);
            Coord unoffCoord = getCoord(evalTAI);
            coordArr.push_back(unoffCoord.offset(toOrientArr[i], fromOrient.getPos(evalTAI), dist.getPos(evalTAI)));
        }
        toOrient.setFromPair(toOrientArr, _tai, DeltaT, true);
        return PVTCoord(coordArr[0], coordArr[1], _tai, DeltaT);
    }

    std::string PVTCoord::__repr__() const {
        std::ostringstream os;
        os << *this;
        return os.str();
    }

    void PVTCoord::_setFromCoordPair(Coord const &coord0, Coord const &coord1, double tai, double deltaT) {
        if (deltaT == 0) {
            throw std::runtime_error("deltaT must be nonzero");
        }
        _coord = coord0;
        _tai = tai;
        _vel = (coord1.getVecPos() - coord0.getVecPos()) / deltaT;
    }

    std::ostream &operator<<(std::ostream &os, PVTCoord const &pvtCoord) {
        Coord coord = pvtCoord.getCoord();
        Eigen::Vector3d vel = pvtCoord.getVel();
        double tai = pvtCoord.getTAI();
        std::ios_base::fmtflags oldFlags = os.flags();
        std::streamsize const oldPrecision = os.precision();
        os << std::fixed
            << "PVTCoord(" << coord
            << ", (" << vel(0) << ", " << vel(1) << ", " << vel(2) << "), "
            << std::setprecision(7) << tai
            << ")" << std::setprecision(oldPrecision);
        os.flags(oldFlags);
        return os;
    }

}
