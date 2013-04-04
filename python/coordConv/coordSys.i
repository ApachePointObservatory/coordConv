%shared_ptr(coordConv::CoordSys);
%shared_ptr(coordConv::MeanCoordSys);
%shared_ptr(coordConv::ApparentCoordSys);
%shared_ptr(coordConv::ICRSCoordSys);
%shared_ptr(coordConv::FK5CoordSys);
%shared_ptr(coordConv::FK4CoordSys);
%shared_ptr(coordConv::GalCoordSys);
%shared_ptr(coordConv::AppGeoCoordSys);
%shared_ptr(coordConv::AppTopoCoordSys);
%shared_ptr(coordConv::ObsCoordSys);
%shared_ptr(coordConv::NoneCoordSys);
%shared_ptr(coordConv::MountCoordSys);

%include "coordConv/coordSys.h"

// note: I tried a %template but got complaints like the following:
// Warning 503: Can't wrap 'coordConv::FK5CoordSys' unless renamed to a valid identifier.
%extend coordConv::ICRSCoordSys {
    %pythoncode {
        def __repr__(self):
            return "%s(%s)" % (self.__class__.__name__, self.getDate())
    }
}
%extend coordConv::FK5CoordSys {
    %pythoncode {
        def __repr__(self):
            return "%s(%s)" % (self.__class__.__name__, self.getDate())
    }
}
%extend coordConv::FK4CoordSys {
    %pythoncode {
        def __repr__(self):
            return "%s(%s)" % (self.__class__.__name__, self.getDate())
    }
}
%extend coordConv::GalCoordSys {
    %pythoncode {
        def __repr__(self):
            return "%s(%s)" % (self.__class__.__name__, self.getDate())
    }
}
%extend coordConv::AppGeoCoordSys {
    %pythoncode {
        def __repr__(self):
            return "%s(%s)" % (self.__class__.__name__, self.getDate())
    }
}
%extend coordConv::AppTopoCoordSys {
    %pythoncode {
        def __repr__(self):
            return "%s(%s)" % (self.__class__.__name__, self.getDate())
    }
}
%extend coordConv::ObsCoordSys {
    %pythoncode {
        def __repr__(self):
            return "%s(%s)" % (self.__class__.__name__, self.getDate())
    }
}
%extend coordConv::NoneCoordSys {
    %pythoncode {
        def __repr__(self):
            return "%s()" % (self.__class__.__name__,)
    }
}
%extend coordConv::MountCoordSys {
    %pythoncode {
        def __repr__(self):
            return "%s(%s)" % (self.__class__.__name__, self.getDate())
    }
}
