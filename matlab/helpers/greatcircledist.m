function rng = greatcircledist(lat1, lon1, lat2, lon2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we use lat in [0,pi], while this function expects [-pi/2,pi/2]
lat1 = lat1-pi/2;
lat2 = lat2-pi/2;

% Calculate great circle distance between points on a sphere using the
% Haversine Formula.  LAT1, LON1, LAT2, and LON2 are in radians.  RNG is a
% length and has the same units as the radius of the sphere, R.  (If R is
% 1, then RNG is effectively arc length in radians.)


a = sin((lat2-lat1)/2).^2 + cos(lat1) .* cos(lat2) .* sin((lon2-lon1)/2).^2;
rng = 2 * atan2(sqrt(a),sqrt(1 - a));
