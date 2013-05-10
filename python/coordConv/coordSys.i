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

%copyctor coordConv::ICRSCoordSys;
%copyctor coordConv::FK5CoordSys;
%copyctor coordConv::FK4CoordSys;
%copyctor coordConv::GalCoordSys;
%copyctor coordConv::AppGeoCoordSys;
%copyctor coordConv::AppTopoCoordSys;
%copyctor coordConv::ObsCoordSys;
%copyctor coordConv::NoneCoordSys;
%copyctor coordConv::MountCoordSys;

%include "coordConv/coordSys.h"
