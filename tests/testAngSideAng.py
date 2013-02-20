#!/usr/bin/env python
import unittest
import numpy
from coordConv import angSideAng, sind, cosd, tand, atan2d

Eps = 2e-16
EpsTest = Eps * 1.001

class TestAngSideAng(unittest.TestCase):
    def testSmallSideA(self):
        """
        a ~ 0, B = various, c various:
        if c nearly 0 (modulo 360): expect C = 90, b = 0, A = 90, unknownAng
        if c nearly 180 (modulo 360): expect C = 90, b = 180, A = 90, unknownAng
        else: expect A = 0, b = a - c, C = 180 - B
        """
        for side_aa in (-Eps, 0.0, Eps):
            for ang_B in (0.0, Eps, 32.0, 97.0, 179.0, 180.0 - Eps, 180.0, 180.0 + Eps, 210.0, 360.0 - Eps, 360.0):
                for side_cc in (180.0, 180.0 - Eps, 179.0, 47.0, Eps, 0.0):
                    if abs(side_cc % 360.0) < EpsTest:
                        expRes = (True, 90.0, 0.0, 90.0)
                    elif abs((side_cc - 180) % 360.0) < EpsTest:
                        expRes = (True, 90.0, 180.0, 90.0)
                    else:
                        expRes = (0.0, side_cc - side_aa, 180.0 - ang_B)
                    self.checkOne((side_aa, ang_B, side_cc), expRes)

    def testBigSideA(self):
        """
        a ~ 180, B = various, c various:
        if c nearly 180 (modulo 360): expect C = 90, b = 0, A = 90, unknownAng
        if c nearly 0 (modulo 360): expect C = 90, b = 180, A = 90, unknownAng
        else: expect A = 180, b = 180 - c, C = B
        """
        for side_aa in (180.0 - Eps, 180.0, 180.0 + Eps):
            for ang_B in (0.0, Eps, 32.0, 97.0, 179.0, 180.0 - Eps, 180.0, 180.0 + Eps, 210.0, 360.0 - Eps, 360.0):
                for side_cc in (180.0, 180.0 - Eps, 179.0, 47.0, Eps, 0.0):
                    if abs((180.0 - side_cc) % 360.0) < EpsTest:
                        expRes = (True, 90.0, 0.0, 90.0)
                    elif abs(side_cc % 360.0) < EpsTest:
                        expRes = (True, 90.0, 180.0, 90.0)
                    else:
                        expRes = (180.0, 180.0 - side_cc, ang_B)
                    self.checkOne((side_aa, ang_B, side_cc), expRes)

    def testSmallSideC(self):
        """
        c ~ 0, B = various, a various:
        if a nearly 0: expect C = 90, b = 0, A = 90, unknownAng
        if a nearly 180: expect C = 90, b = 180, A = 90, unknownAng
        else: expect A = 180 - B, b = a, C = 0
        """
        for side_cc in (0.0, Eps):
            for ang_B in (0.0, Eps, 32.0, 97.0, 179.0, 180.0 - Eps, 180.0, 180.0 + Eps, 210.0, 360.0 - Eps, 360.0):
                for side_aa in (180.0, 180.0 - Eps, 179.0, 47.0, Eps, 0.0):
                    if abs(side_aa % 360.0) < EpsTest:
                        expRes = (True, 90.0, 0.0, 90.0)
                    elif abs((180.0 - side_aa) % 360.0) < EpsTest:
                        expRes = (True, 90.0, 180.0, 90.0)
                    else:
                        expRes = (180.0 - ang_B, side_aa, 0.0)
                    self.checkOne((side_aa, ang_B, side_cc), expRes)

    def testLargeSideC(self):
        """
        c ~ 180, B = various, a various:
        if a nearly 0 (modulo 360): expect C = 90, b = 180, A = 90, unknownAng
        if a nearly 180 (modulo 360): expect C = 90, b = 0, A = 90, unknownAng
        else: expect A = 180, b = 180 - c, C = B
        """
        for side_cc in (180.0 - Eps, 180.0):
            for ang_B in (0.0, Eps, 32.0, 97.0, 179.0, 180.0 - Eps, 180.0, 180.0 + Eps, 210.0, 360.0 - Eps, 360.0):
                for side_aa in (180.0, 180.0 - Eps, 179.0, 47.0, Eps, 0.0):
                    if side_aa < EpsTest:
                        expRes = (True, 90.0, 180.0, 90.0)
                    elif 180.0 - side_aa < EpsTest:
                        expRes = (True, 90.0, 0.0, 90.0)
                    else:
                        expRes = (ang_B, 180.0 - side_aa, 180.0)
                    self.checkOne((side_aa, ang_B, side_cc), expRes)

    def testSideA90(self):
        """
        a = 90, B varies but not nearly 0 or 360, c fairly small but >> Eps
        expect: A = 180 - B, b = a + c cos(B), C ~= 0
        """
        side_aa = 90.0
        for side_cc in (1.0e-12, 1.0e-10):
            for ang_B in (23, 90, 180 - Eps, 180, 180 + Eps, 256, 359):
                expRes = (180.0 - ang_B, side_aa + (side_cc * cosd(ang_B)), 0.0)
                self.checkOne((side_aa, ang_B, side_cc), expRes)

    def testSideC90(self):
        """
        a fairly small but >> Eps, B varies, c = 90
        expect: C = 180 - B, b = c + a cos(B), A ~= 0
        """
        side_cc = 90.0
        for side_aa in (1.0e-12, 1.0e-10):
            for ang_B in (23, 90, 180 - Eps, 180, 180 + Eps, 256, 359):
                expRes = (0.0, side_cc + (side_aa * cosd(ang_B)), 180.0 - ang_B)
                self.checkOne((side_aa, ang_B, side_cc), expRes)

    def testAngBSmall(self):
        """
        B small, a = any not small, c = any not small:
        if c != a: expect A = 90, b = 0, C = 90, unknown
        if c << a: expect A = 180,   b = c - a, C = 0
        if c >> a: expect A = 0, b = a - c, C = 180
        """
        for side_aa in (179.9, -27.0, 27.0, 0.1):
            for side_cc in (side_aa - 45.0, side_aa - Eps, side_aa, side_aa + Eps, side_aa + 45.0):
                if abs(side_cc - side_aa) < EpsTest:
                    expRes = (True, 90.0, 0.0, 90.0)
                elif side_cc < side_aa:
                    expRes = (180.0, side_aa - side_cc, 0.0)
                else:
                    expRes = (0.0, side_cc - side_aa, 180.0)
                for ang_B in (-Eps, 0.0, Eps):
                    self.checkOne((side_aa, ang_B, side_cc), expRes)

    def testRightTriangle(self):
        """
        right triangle: B = 90, a and c vary but avoid poles
        tan C = tan c / sin a
        tan c = (tan a / sinA * sinb)
        with some tweaks to handle the other quadrants
        """
        ang_B = 90.0
        for side_aa in (1.0, 20.0, 45.0, 90, 110.0, 179.0):
            for side_cc in (1.0, 20.0, 45.0, 90.0, 110.0, 179.0):
                ang_A = atan2d(tand(side_aa), sind(side_cc))
                ang_C = atan2d(tand(side_cc), sind(side_aa))
                side_bb = atan2d(tand(side_aa), sind(ang_A) * cosd(side_cc))
                # these tweaks handle other quadrants; they're based on what works, so are somewhat suspect
                if side_bb < 0:
                    side_bb = - side_bb
                if ang_A < 0:
                    ang_A = 180.0 + ang_A
                if ang_C < 0:
                    ang_C = 180.0 + ang_C
                self.checkOne((side_aa, ang_B, side_cc), (ang_A, side_bb, ang_C))

    def testSpecialCases(self):
        """A few hard-coded special cases
        """
        for testInput, expectedOutput in [
            # 90/90/90 triangle
            ((90, 90, 90), (90, 90, 90)),
    
            # inputs that might cause side_bb < 0, (but should not)
            ((45, 1, 45), (89.6464421219342, 0.707102293688337, 89.6464421219342)),
            ((45, -1, 45), (270.353557878066, 0.707102293688337, 270.353557878066)),
            ((135, 1, 135), (90.3535578780658, 0.707102293688337, 90.3535578780658)),
            ((135, -1, 135), (269.646442121934, 0.707102293688308, 269.646442121934)),
        ]:
            self.checkOne(testInput, expectedOutput)

    def checkOne(self, testInput, expectedOutput):
        """Check one case
        
        Inputs:
        - testInput: a vector of inputs (ang, side, ang)
        - expectedOutput: a vector of outputs;
            if length(3) then (side, ang, side) and unknownAng assumed False
            if length(r) then (unknownAng, side, ang, side)
        """
        if len(expectedOutput) == 3:
            expectedOutput = (False,) + tuple(expectedOutput)
        elif len(expectedOutput) != 4:
            raise RuntimeError("len(expectedOutput) = %d; must be 3 or 4" % (len(expectedOutput),))

        actualOutput = angSideAng(*testInput)

        # to handle angles comparing things like 359.999... to 0, compare sin and cos of ang_A and ang_C:
        procExpected = processOutput(expectedOutput)
        procActual = processOutput(actualOutput)
        if not numpy.allclose(procExpected, procActual, rtol=1.0e-10, atol=1.0e-10):
            self.fail("failed on input: %s; expected output = %s; actual output = %s" % \
                (testInput, expectedOutput, actualOutput))
        if actualOutput[1] < 0.0 or actualOutput[1] >= 360.0 \
            or actualOutput[2] < 0.0 or actualOutput[2] >= 360.0 \
            or actualOutput[3] < 0.0 or actualOutput[3] >= 360.0:
            self.fail("failed on input %s; one or more output angles out of range: %s" % \
                (testInput, actualOutput))
    
def processOutput(outputVec):
    return (
        outputVec[0],
        sind(outputVec[1]), cosd(outputVec[2]),
        outputVec[2],
        sind(outputVec[3]), cosd(outputVec[3]),
    )

if __name__ == '__main__':
    unittest.main()
