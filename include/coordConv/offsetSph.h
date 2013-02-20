#pragma once

#include "Eigen/Dense"

namespace coordConv {

    /**
    Offset a spherical coordinate along a great cicle

    @param[out] destEquatAng: resulting equatorial angle (e.g. RA) (deg)
    @param[out] destPolarAng: resulting polar angle (e.g. Dec) (deg)
    @param[out] destOrient: angle of great circle at destination (deg)
        (angle from increasing toSphPos(0) to great circle)
    @param[in] srcEquatAng: initial equatorial angle (e.g. RA) (deg)
    @param[in] srcPolarAng: initial polar angle (e.g. Dec) (deg)
    @param[in] srcOrient: angle of great circle at source (deg)
        (angle from direction of increasing fromSphPos(0) to great circle)
    @param[in] dist: offset length (great circle arg length) (deg)

    @raise runtime_error if either srcPolarAng or destPolarAng too near a pole
    
    This diagram may help:

                                                            .
                                                      .      ⎞ toOrient
                                               *--------------> dir of increasing destEquatAng
                                       .      dest
                              .
                    .
           .         ⎞ fromOrient
     *----------------> dir of increasing srcEquatAng
    src
    */
    void offsetSph(
        double &destEquatAng,
        double &destPolarAng,
        double &destOrient,
        double srcEquatAng,
        double srcPolarAng,
        double srcOrient,
        double dist
    );

}
