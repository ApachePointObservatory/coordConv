%include "coordConv/pvt.h"

%extend coordConv::PVT {
    %pythoncode {
        def __str__(self):
            return "PVT(%s, %s, %s)" % (self.pos, self.vel, self.t)
    }
}
