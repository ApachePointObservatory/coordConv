#include <stdexcept>
#include "coordConv/pvtCoord.h"

static const double DeltaT = 0.01;

namespace coordConv {
    
    PVTCoord::PVTCoord(Coord const &coord0, Coord const &coord1, double tai, double deltaT) :
        _tai(tai),
        _pos(coord0.getVecPos()),
        _vel((coord1.getVecPos() - _pos) / deltaT),
        _pm(coord0.getVecPM())
    {
        if (deltaT == 0) {
            throw std::runtime_error("deltaT must be nonzero");
        }
    }

    PVTCoord::PVTCoord(Coord const &coord, double tai) :
        _tai(tai),
        _pos(coord.getVecPos()),
        _vel(),
        _pm(coord.getVecPM())
    {
        _vel.setZero();
    }
    
    bool PVTCoord::getSphPVT(PVT &equatPVT, PVT &polarPVT) const {
        Coord currCoord(_pos, _pm);
        Coord laterCoord = getCoord(_tai + DeltaT);

        double equatPos[2], polarPos[2];
        bool atPole = currCoord.getSphPos(equatPos[0], polarPos[0]);
        atPole |= laterCoord.getSphPos(equatPos[1], polarPos[1]);
        
        equatPVT.setFromAnglePair(equatPos, _tai, DeltaT);
        polarPVT.setFromAnglePair(polarPos, _tai, DeltaT);
        return atPole;
    }

}
