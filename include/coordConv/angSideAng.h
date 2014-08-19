#pragma once

#include <cmath>

namespace coordConv {

    /**
    Solve for two sides and the included angle of a spherical triangle

    @param[out] angA  interior angle opposite side A (deg)
    @param[out] sideB  length of side B (deg)
    @param[out] angC  interior angle opposite side C (deg)
    @param[in] sideA  length of side A (deg)
    @param[in] angB  interior angle opposite side B (deg)
    @param[in] sideC  length of side C (deg)

    @return unknownAng: if true: sideB is so near 0 or 180 (see Special Cases below) that angA and angC cannot
    not be computed; angA and angC are set to 90 (hence their sum is 180, which is essentially correct)
    and sideB will be 0 or 180 (again, essentially correct).

    @throw runtime_error if the inputs are too small to allow computation

    Special Cases (in the order they are handled):
     sideA    angB    sideC       angA        sideB        angC
    ----------------------------------------------------------------
      ~0      any      ~0     unknown(90)       0        unknown(90)
      ~0      any     ~180    unknown(90)      180       unknown(90)
      ~0      any     !pole        0          sideC       180-angB

     ~180     any     ~0      unknown(90)      180       unknown(90)
     ~180     any    ~180     unknown(90)       0        unknown(90)
     ~180     any    !pole        180       180-sideC       angB

     !pole    any     ~0       180-angB       sideA          0
     !pole    any    ~180        angB       180-sideA       180

      any     ~0   ~=sideA    unknown(90)       0        unknown(90)
      any     ~0    <sideA        180      sideA-sideC       0
      any     ~0    >sideA         0       sideC-sideA      180

    where:
    - !pole means not nearly 0 and not nearly 180 (modulo 360)
    - unknown(90) means unknownAng = true is returned
    - All relations are modulo 360. For example ~0 means approximately zero, 360, etc.

    @note: allowing angles in the 3rd and 4th quadrants is unusual.

    References:
    Selby, Standard Math Tables, crc, 15th ed, 1967, p161 (Spherical Trig.)
    */
    bool angSideAng(double &angA, double &sideB, double &angC, double sideA, double angB, double sideC);

}
