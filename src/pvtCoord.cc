#include <cmath>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <vector>
#include <tr1/array>

#include "coordConv/mathUtils.h"
#include "coordConv/pvtCoord.h"

static const double DeltaT = 0.01;

namespace coordConv {

    PVTCoord::PVTCoord(Coord const &coord, double orient, double vel, double tai) :
        _coord(coord),
        _orient(orient),
        _vel(vel),
        _tai(tai)
    {
        _checkCoord();
    }
    
    PVTCoord::PVTCoord(Coord const &coord0, Coord const &coord1, double tai, double deltaT, double defOrient) :
        _coord(coord0),
        _orient(defOrient),    // assume default for now
        _vel(0),            // assume default for now
        _tai(tai)
    {
        if (deltaT == 0) {
            throw std::runtime_error("deltaT must be nonzero");
        }
        double dist = coord0.angularSeparation(coord1);
        double orient = coord0.orientationTo(coord1);
        if (std::isfinite(orient)) {
            _vel = dist / deltaT;
            _orient = orient;
        }
        _checkCoord();
    }

    PVTCoord::PVTCoord(PVT const &equatPVT, PVT const &polarPVT, double tai, double parallax, double defOrient) {
        double polarPos = polarPVT.getPos(tai);
        _coord = Coord(equatPVT.getPos(tai), polarPos, parallax);
        _tai = tai;
        _setOrientVelFromSph(polarPos, equatPVT.vel, polarPVT.vel, defOrient);
        _checkCoord();
    }

    PVTCoord::PVTCoord(PVT const &equatPVT, PVT const &polarPVT, double tai, double parallax, double equatPM, double polarPM, double radVel, double defOrient) {
        double polarPos = polarPVT.getPos(tai);
        _coord = Coord(equatPVT.getPos(tai), polarPos, parallax, equatPM, polarPM, radVel);
        _tai = tai;
        _setOrientVelFromSph(polarPos, equatPVT.vel, polarPVT.vel, defOrient);
        _checkCoord();
    }

    PVTCoord::PVTCoord() :
        _coord(Coord()),
        _orient(DoubleNaN),
        _vel(DoubleNaN),
        _tai(DoubleNaN)
    { };
    
    Coord PVTCoord::getCoord(double tai) const {
        if (_vel == 0.0) {
            return _coord;
        } else {
            double dist = _vel * (tai - _tai);
            double newOrient;
            return _coord.offset(newOrient, _orient, dist);
        }
    }    
    
    bool PVTCoord::getSphPVT(PVT &equatPVT, PVT &polarPVT, double tai) const {
        double equatPos[2], polarPos[2];
        bool atPole = false;
        for (int i = 0; i < 2; ++i) {
            double evalDate = tai + (i * DeltaT);
            Coord coord = getCoord(evalDate);
            atPole |= coord.getSphPos(equatPos[i], polarPos[i]);
        }
        equatPVT.setFromPair(equatPos, tai, DeltaT, true);
        polarPVT.setFromPair(polarPos, tai, DeltaT, false);
        return atPole;
    }
    
    bool PVTCoord::isfinite() const {
        return _coord.isfinite() && std::isfinite(_orient) && std::isfinite(_vel) && std::isfinite(_tai);
    }

    PVT PVTCoord::angularSeparation(PVTCoord const &pvtCoord, double tai) const {
        double posArr[2];
        for (int i = 0; i < 2; ++i) {
            double tempTAI = tai + (i * DeltaT);
            Coord thisCoord = getCoord(tempTAI);
            Coord otherCoord = pvtCoord.getCoord(tempTAI);
            posArr[i] = thisCoord.angularSeparation(otherCoord);
        }
        PVT res = PVT();
        res.setFromPair(posArr, tai, DeltaT, false);
        return res;
    }

    PVT PVTCoord::orientationTo(PVTCoord const &pvtCoord, double tai) const {
        double posArr[2];
        for (int i = 0; i < 2; ++i) {
            double tempTAI = tai + (i * DeltaT);
            Coord thisCoord = getCoord(tempTAI);
            Coord otherCoord = pvtCoord.getCoord(tempTAI);
            posArr[i] = thisCoord.orientationTo(otherCoord);
        }
        PVT res = PVT();
        // if the orientation is only finite at one of the two times
        // then the distance is 0 at the other time and the orientation is fixed
        // (since the relative velocity vector goes through this PVT)
        if (std::isfinite(posArr[0]) && std::isfinite(posArr[1])) {
            res.setFromPair(posArr, tai, DeltaT, true);
        } else if (std::isfinite(posArr[0])) {
            res.pos = posArr[0];
            res.vel = 0;
            res.t = tai;
        } else if (std::isfinite(posArr[1])) {
            res.pos = posArr[1];
            res.vel = 0;
            res.t = tai;
        } else {
            res.pos = DoubleNaN;
            res.vel = DoubleNaN;
            res.t = tai;
        }
        return res;
    }

    PVTCoord PVTCoord::offset(PVT &toOrient, PVT const &fromOrient, PVT const &dist, double tai) const {
        std::vector<Coord> coordArr;
        double toOrientArr[2];
        for (int i = 0; i < 2; ++i) {
            double tempTAI = tai + (i * DeltaT);
            Coord unoffCoord = getCoord(tempTAI);
            coordArr.push_back(unoffCoord.offset(toOrientArr[i], fromOrient.getPos(tempTAI), dist.getPos(tempTAI)));
        }
        toOrient.setFromPair(toOrientArr, tai, DeltaT, true);
        return PVTCoord(coordArr[0], coordArr[1], tai, DeltaT, fromOrient.getPos(tai));
    }

    std::string PVTCoord::__repr__() const {
        std::ostringstream os;
        os << *this;
        return os.str();
    }

    void PVTCoord::_setOrientVelFromSph(double polarPos, double equatVel, double polarVel, double defOrient) {
        if (_coord.atPole()) {
            if ((std::abs(polarVel) > DoubleEpsilon) || (std::abs(equatVel) > DoubleEpsilon)) {
                throw std::runtime_error("Too near pole; cannot specify nonzero velocity");
            }
            _orient = defOrient;
            _vel = 0;
        } else {
            double equatSpaceVel = equatVel * cosd(polarPos);
            _orient = atan2d(polarVel, equatSpaceVel);
            _vel = hypot(polarVel, equatSpaceVel);
        }
    }

    void PVTCoord::_checkCoord() const {
        if (_coord.atPole() and _vel != 0) {
            throw std::runtime_error("Invalid PVTCoord: coord at pole but velocity nonzero");
        }
    }

    std::ostream &operator<<(std::ostream &os, PVTCoord const &pvtCoord) {
        Coord coord = pvtCoord.getCoord();
        double orient = pvtCoord.getOrient();
        double vel = pvtCoord.getVel();
        double tai = pvtCoord.getTAI();
        std::ios_base::fmtflags oldFlags = os.flags();
        std::streamsize const oldPrecision = os.precision();
        os << std::fixed
            << "PVTCoord(" << coord << ", " 
            << std::setprecision(5) << orient << ", " << vel << ", "
            << std::setprecision(7) << tai
            << ")" << std::setprecision(oldPrecision);
        os.flags(oldFlags);
        return os;
    }

}
