#include <cmath>
#include <stdexcept>
#include <vector>
#include <tr1/array>
#include "coordConv/pvtCoord.h"
#include "coordConv/mathUtils.h"

static const double DeltaT = 0.01;

namespace coordConv {

    PVTCoord::PVTCoord(Coord const &coord, double orient, double vel, double tai) :
        _coord(coord),
        _orient(orient),
        _vel(vel),
        _tai(tai)
    {}
    
    PVTCoord::PVTCoord(Coord const &coord0, Coord const &coord1, double tai, double deltaT, double defOrient) :
        _coord(coord0),
        _orient(defOrient), // assume insufficient velocity to compute orientation
        _vel(0),            // assume insufficient velocity to compute orientation
        _tai(tai)
    {
        if (deltaT == 0) {
            throw std::runtime_error("deltaT must be nonzero");
        }
        double dist = coord0.angularSeparation(coord1);
        double vel = dist / deltaT;
        double orient = coord0.orientationTo(coord1);
        if (std::isfinite(orient)) {
            _orient = orient;
            _vel = vel;
        }
    }

    PVTCoord::PVTCoord(PVT const &polarPVT, PVT const &equatPVT, double tai, double parallax, double defOrient) {
        _coord = Coord(polarPVT.getPos(tai), equatPVT.getPos(tai), parallax);
        _vel = hypot(equatPVT.vel, polarPVT.vel);
        if (std::abs(_vel) > DoubleEpsilon) {
            _orient = atan2d(equatPVT.vel, polarPVT.vel);
        } else {
            _orient = defOrient;
            _vel = 0;
        }
        _tai = tai;
    }

    PVTCoord::PVTCoord(PVT const &polarPVT, PVT const &equatPVT, double tai, double parallax, double equatPM, double polarPM, double radVel, double defOrient) {
        _coord = Coord(polarPVT.getPos(tai), equatPVT.getPos(tai), parallax, equatPM, polarPM, radVel);
        _vel = hypot(equatPVT.vel, polarPVT.vel);
        if (std::abs(_vel) > DoubleEpsilon) {
            _orient = atan2d(equatPVT.vel, polarPVT.vel);
        } else {
            _orient = defOrient;
            _vel = 0;
        }
        _tai = tai;
    }

    PVTCoord::PVTCoord() :
        _coord(Coord()),
        _orient(DoubleNaN),
        _vel(DoubleNaN),
        _tai(DoubleNaN)
    { };
    
    Coord PVTCoord::getCoord(double tai) const {
        double dist = _vel * (tai - _tai);
        double newOrient;
        return _coord.offset(newOrient, _orient, dist);
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
        res.setFromPair(posArr, tai, DeltaT, true);
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

    std::ostream &operator<<(std::ostream &os, PVTCoord const &pvtCoord) {
        Coord coord = pvtCoord.getCoord();
        double orient = pvtCoord.getOrient();
        double vel = pvtCoord.getVel();
        double tai = pvtCoord.getTAI();
        os << "PVTCoord(" << coord << ", " << orient << ", " << vel << ", " << tai << ")";
        return os;
    }

}
