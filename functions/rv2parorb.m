function [a,E,i,OM,om,teta] = rv2parorb(rr,vv, mu)

%transformation from Cartesian state to orbital elements
%
%[a, e, i, OM, om, theta] = rv2parorb (r, v, mu)
%
%input arguments: 

% rr  [3x1] posizion vector
% vv  [3x1] velocity vector
% mu  [3x1] gravitational parameter

%output arguments: 
% a     [1x1]  semi-major axis
% e     [1x1]  eccentricity
% i     [1x1]  inclination
% OM    [1x1]  RAAN
% om    [1x1]  pericenter anomaly
% theta [1x1]  true anomaly

II = [1 0 0]';
JJ = [0 1 0]';
KK = [0 0 1]';

R = norm(rr);
V = norm(vv);

a = 1/(2/R - V^2/mu);

hh = (cross(rr, vv));
H = norm(hh);

ee = ((1/mu) * (cross(vv,hh)) - rr/R);
E = norm(ee);

i = acos (hh(3)/H);

nn = (cross(KK,hh))/(norm((cross(KK,hh))));

if (nn(2) >= 0)
    OM = acos(nn(1));
else
    OM = 2*pi - acos(nn(1));
end

if (ee(3) >= 0)
    om = acos((dot(ee,nn))/E);
else
    om = 2*pi - acos((dot(ee,nn))/E);
end

if (dot(vv,rr) >= 0)
    teta = acos((dot(rr,ee))/(R*E));
else
    teta = 2 *pi - acos((dot(rr,ee))/(R*E));

end