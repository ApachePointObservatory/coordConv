#include <stdexcept>
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
    
    Coord PVTCoord::getCoord(double tai) const {
        double dist = _vel * (tai - _tai);
        double newOrient;
        return _coord.offset(newOrient, _orient, dist);
    }    

//     PVTCoord PVTCoord::getPVTCoord(double tai) const {
//         double dist = _vel * (tai - _tai);
//         double newOrient;
//         Coord newCoord = _coord.offset(newOrient, _orient, dist);
//         return PVTCoord(newCoord, newOrient, _vel, tai);
//     }    

    bool PVTCoord::getSphPVT(PVT &equatPVT, PVT &polarPVT) const {
        double equatPos[2], polarPos[2];
        Coord laterCoord = getCoord(_tai + DeltaT);
        bool atPole = _coord.getSphPos(equatPos[0], polarPos[0]);
        atPole |= laterCoord.getSphPos(equatPos[1], polarPos[1]);
        
        equatPVT.setFromAnglePair(equatPos, _tai, DeltaT);
        polarPVT.setFromAnglePair(polarPos, _tai, DeltaT);
        return atPole;
    }

}
