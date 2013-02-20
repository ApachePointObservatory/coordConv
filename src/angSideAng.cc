#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include "coordConv/mathUtils.h"

namespace coordConv {

    bool angSideAng(double &angA, double &sideB, double &angC, double sideA, double angB, double sideC) {
        double const Epsilon = std::numeric_limits<double>::epsilon();

        bool unknownAng = false; // assume we can compute angA and angC

        double sinHalfAngB  = sind(angB  * 0.5);
        double cosHalfAngB  = cosd(angB  * 0.5);
        double sinHalfSideA = sind(sideA * 0.5);
        double cosHalfSideA = cosd(sideA * 0.5);
        double sinHalfSideC = sind(sideC * 0.5);
        double cosHalfSideC = cosd(sideC * 0.5);

        if (std::abs(sinHalfSideA) < Epsilon) {
            // sideA is nearly zero (modulo 360)
            if (std::abs(sinHalfSideC) < Epsilon) {
                // sideC is nearly 0 (modulo 360)
                angA = 90.0;
                sideB = 0.0;
                angC = 90.0;
                unknownAng = true;
            } else if (std::abs(cosHalfSideC) < Epsilon) {
                // sideC is nearly 180 (modulo 360)
                angA = 90.0;
                sideB = 180.0;
                angC = 90.0;
                unknownAng = true;
            } else {
                // sideC is not nearly 0 or 180
                angA = 0.0;
                sideB = sideC;
                angC = 180.0 - angB;
            }
        } else if (std::abs(cosHalfSideA) < Epsilon) {
            // sideA is nearly 180 (modulo 360)
            if (std::abs(cosHalfSideC) < Epsilon) {
                // sideC is nearly 180 (modulo 360)
                angA = 90.0;
                sideB = 0.0;
                angC = 90.0;
                unknownAng = true;
            } else if (std::abs(sinHalfSideC) < Epsilon) {
                // sideC is nearly 0 (modulo 360)
                angA = 90.0;
                sideB = 180.0;
                angC = 90.0;
                unknownAng = true;
            } else {
                // sideC is not nearly 0 or 180 (modulo 360)
                angA = 180.0;
                sideB = 180.0 - sideC;
                angC = angB;
            }
        } else if (std::abs(sinHalfSideC) < Epsilon) {
            // sideC is nearly zero (modulo 360) and sideA is not
            angA = 180.0 - angB;
            sideB = sideA;
            angC = 0.0;
        } else if (std::abs(cosHalfSideC) < Epsilon) {
            // sideC is nearly 180 (modulo 360) and sideA is not
            angA = angB;
            sideB = 180.0 - sideA;
            angC = 180.0;
        } else if (std::abs(sinHalfAngB) < Epsilon) {
            // B is nearly 0 (modulo 360)
            if (std::abs(sideA - sideC) < Epsilon) {
                // angB ~= 0 (modulo 360) and sideA ~= sideC (modulo 360); cannot compute angA or angC:
                angA = 90.0;
                sideB = 0.0;
                angC = 90.0;
                unknownAng = true;
            } else if (sideC < sideA) {
                angA = 180.0;
                sideB = sideA - sideC;
                angC = 0.0;
            } else {
                angA = 0.0;
                sideB = sideC - sideA;
                angC = 180.0;
            }
        } else {
            // +
            // compute angA and angC using Napier's analogies
            // -
            // compute numerator and denominator, where tan((a +/- c) / 2) = num/den
            double num1 = cosHalfAngB * ((cosHalfSideA * cosHalfSideC) + (sinHalfSideA * sinHalfSideC));
            double den1 = sinHalfAngB * ((cosHalfSideA * cosHalfSideC) - (sinHalfSideA * sinHalfSideC));
            double num2 = cosHalfAngB * ((sinHalfSideA * cosHalfSideC) - (cosHalfSideA * sinHalfSideC));
            double den2 = sinHalfAngB * ((sinHalfSideA * cosHalfSideC) + (cosHalfSideA * sinHalfSideC));

            // if numerator and denominator are too small
            // to accurately determine angle = atan2(num, den), give up
            if (   ((std::abs(num1) <= Epsilon) && (std::abs(den1) <= Epsilon))
                || ((std::abs(num2) <= Epsilon) && (std::abs(den2) <= Epsilon))) {
                std::ostringstream os;
                os << "Bug: can't compute angA and angC with sideA=" << sideA << ", angB=" << angB << ", sideC=" << sideC;
                throw std::runtime_error(os.str());
            }

            // compute (a +/- c) / 2, and use to compute angles a and c
            double h_sum_AC  = atan2d (num1, den1);
            double h_diff_AC = atan2d (num2, den2);

            angA = h_sum_AC + h_diff_AC;
            angC = h_sum_AC - h_diff_AC;

            // +
            // compute sideB using one of two Napier's analogies
            // (one is for sideB - sideA, one for sideB + sideA)
            // -
            double sinHalfAngA = sind(angA * 0.5);
            double cosHalfAngA = cosd(angA * 0.5);

            // numerator and denominator for Napier's analogy for sideB - sideA
            double num3 = sinHalfSideC * ((sinHalfAngB * cosHalfAngA) - (cosHalfAngB * sinHalfAngA));
            double den3 = cosHalfSideC * ((sinHalfAngB * cosHalfAngA) + (cosHalfAngB * sinHalfAngA));

            // numerator and denominator for Napier's analogy for sideB + sideA
            double num4 = sinHalfSideC * ((cosHalfAngB * cosHalfAngA) + (sinHalfAngB * sinHalfAngA));
            double den4 = cosHalfSideC * ((cosHalfAngB * cosHalfAngA) - (sinHalfAngB * sinHalfAngA));

            // pick the more suitable Napier's analogy to compute sideB
            if (std::abs(num3) + std::abs(den3) > std::abs(num4) + std::abs(den4)) {
                // use Napier's analogy for sideB - sideA
                sideB = 2.0 * atan2d(num3, den3) + sideA;
            } else {
                sideB = 2.0 * atan2d(num4, den4) - sideA;
            }
        }
        angA = wrapPos(angA);
        sideB = wrapPos(sideB);
        angC = wrapPos(angC);
        return unknownAng;
    }

}
