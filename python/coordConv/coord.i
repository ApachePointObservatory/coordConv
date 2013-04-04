%include "coordConv/coord.h"

%extend coordConv::Coord {
    %pythoncode {
        def __repr__(self):
            atPole, equatAng, polarAng = self.getSphPos()
            atPole, equatPM, polarPM = self.getPM()
            radVel = self.getRadVel()
            parallax = self.getParallax()
            return "Coord(equatAng=%s, polarAng=%s, parallax=%s, equatPM=%s, polarPM=%s, radVel=%s)" % \
                (equatAng, polarAng, parallax, equatPM, polarPM, radVel)

        def __str__(self):
            atPole, equatAng, polarAng = self.getSphPos()
            return "Coord(%s, %s)" % (equatAng, polarAng)
    }
}
