#include <cmath>
#include <stdexcept>
#include <vector>
#include "coordConv/pvtCoord.h"

static const double DeltaT = 0.01;

namespace coordConv {

    PVTCoord::PVTCoord(Coord const &coord, double orient, double vel, double tai) :
        _coord(coord),
        _orient(orient),
        _vel(vel),
        _tai(tai)
    {}
    
    PVTCoord::PVTCoord(Coord const &coord0, Coord const &coord1, double tai, double deltaT) :
        _coord(coord0),
        _orient(),
        _vel(),
        _tai(tai)
    {
        if (deltaT == 0) {
            throw std::runtime_error("deltaT must be nonzero");
        }
        double dist = coord0.angularSeparation(coord1);
        _vel = dist / deltaT;
        _orient = coord0.orientationTo(coord1);
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

    bool PVTCoord::getSphPVT(PVT &equatPVT, PVT &polarPVT) const {
        double equatPos[2], polarPos[2];
        Coord laterCoord = getCoord(_tai + DeltaT);
        bool atPole = _coord.getSphPos(equatPos[0], polarPos[0]);
        atPole |= laterCoord.getSphPos(equatPos[1], polarPos[1]);
        
        equatPVT.setFromPair(equatPos, _tai, DeltaT, true);
        polarPVT.setFromPair(polarPos, _tai, DeltaT, false);
        return atPole;
    }
    
    bool PVTCoord::isfinite() const {
        return _coord.isfinite() && std::isfinite(_orient) && std::isfinite(_vel) && std::isfinite(_tai);
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
        return PVTCoord(coordArr[0], coordArr[1], tai, DeltaT);
    }

}
