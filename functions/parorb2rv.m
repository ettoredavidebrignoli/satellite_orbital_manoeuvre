function [rr,vv] = parorb2rv(a,e,i,OM,om,teta,mu)

%transformation from orbital elements to Cartesiane state
%
%
%input
%---------------------------------------------------------------------

% a     [1x1]  semi-major axis            [Km]
% e     [1x1]  eccentricity
% i     [1x1]  inclination                [rad]
% OM    [1x1]  RAAN                       [rad]
% om    [1x1]  pericenter anomaly         [rad]
% theta [1x1]  true anomaly               [rad]
% mu    [1x1] gravitational parameter    [Km^3/s^2]

%output
%-----------------------------------------------------------------------
% rr  [3x1]    posizion vector           [Km]
% vv  [3x1]     velocity vector          [Km/s]



R_OM =[cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];
R_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
R_om = [cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];

T = (R_om)*(R_i)*(R_OM);

x = (a*(1-e^2))/(1 + e*cos(teta))*cos(teta);
y = (a*(1-e^2))/(1 + e*cos(teta))*sin(teta);
z = 0.*teta;
rr = [x,y,z]';

vr = (sqrt(mu/(a*(1-e^2))))*e*sin(teta);
vteta = (sqrt(mu/(a*(1-e^2))))*(1+e*cos(teta));

vx = vr*cos(teta) - vteta*sin(teta);
vy = vr*sin(teta) + vteta*cos(teta);
vz = 0.*teta;
vv = [vx,vy,vz]';

rr = T'*rr;
vv = T'*vv;
end




