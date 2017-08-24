function dist = pos2dist(lat1,lon1,lat2,lon2,method)
% function dist = pos2dist(lat1,lon1,lat2,lon2,method)
% calculate distance between two points on earth's surface
% given by their latitude-longitude pair.
% Input lat1,lon1,lat2,lon2 are in degrees, without 'NSWE' indicators.
% Input method is 1 or 2. Default is 1.
% Method 1 uses plane approximation,
% only for points within several tens of kilometers (angles in rads):
% d =
% sqrt(R_equator^2*(lat1-lat2)^2 + R_polar^2*(lon1-lon2)^2*cos((lat1+lat2)/2)^2)
% Method 2 calculates sphereic geodesic distance for points farther apart,
% but ignores flattening of the earth:
% d =
% R_aver * acos(cos(lat1)cos(lat2)cos(lon1-lon2)+sin(lat1)sin(lat2))
% Output dist is in km.
% Returns -99999 if input argument(s) is/are incorrect.
% Flora Sun, University of Toronto, Jun 12, 2004.
if nargin < 4
    dist = -99999;
    disp('Number of input arguments error! distance = -99999');
    return;
end
if abs(lat1)>90 | abs(lat2)>90 | abs(lon1)>360 | abs(lon2)>360
    dist = -99999;
    disp('Degree(s) illegal! distance = -99999');
    return;
end
if lon1 < 0
    lon1 = lon1 + 360;
end
if lon2 < 0
    lon2 = lon2 + 360;
end
% Default method is 1.
if nargin == 4
    method = 1;
end
if method == 1
    km_per_deg_la = 111.3237;
    km_per_deg_lo = 111.1350;
    km_la = km_per_deg_la * (lat1-lat2);
    % Always calculate the shorter arc.
    if abs(lon1-lon2) > 180
        dif_lo = abs(lon1-lon2)-180;
    else
        dif_lo = abs(lon1-lon2);
    end
    km_lo = km_per_deg_lo * dif_lo * cos((lat1+lat2)*pi/360);
    dist = sqrt(km_la^2 + km_lo^2);
else
    R_aver = 6374;
    deg2rad = pi/180;
    lat1 = lat1 * deg2rad;
    lon1 = lon1 * deg2rad;
    lat2 = lat2 * deg2rad;
    lon2 = lon2 * deg2rad;
    dist = R_aver * acos(cos(lat1)*cos(lat2)*cos(lon1-lon2) + sin(lat1)*sin(lat2));
end
