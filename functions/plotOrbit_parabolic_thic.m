function plotOrbit_parabolic_thic(rp,i,OM,om,th_i,th_f,dt,mu)
% plot of the orbit given the parameters with Linewidth 2 
%input
%---------------------------------------------------------------------
% rp    [1x1]  latus rectum               [Km]
% i     [1x1]  inclination                [rad]
% OM    [1x1]  RAAN                       [rad]
% om    [1x1]  pericenter anomaly         [rad]
% th_i  [1x1]  initial true anomaly       [rad]
% th_f  [1x1]  final true anomaly       [rad]
% dt    [1x1]  step
% mu    [1x1]  gravitational parameter    [Km^3/s^2]
% r     'r'     plot color

if ( th_i > th_f )
     th_f = th_f + 2*pi;
end

theta = (th_i:dt:th_f);
x = [];
y = [];
z = [];

x =( (2*rp)./(1 + cos(theta)) ).*cos(theta);
y =( (2*rp)./(1 + cos(theta)) ).*sin(theta);
z = 0.*theta;

R_OM =[cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];
R_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
R_om = [cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];

T = (R_om)*(R_i)*(R_OM);

rr = [x;y;z];
rr = T'*rr;

plot3(rr(1,:), rr(2,:), rr(3,:),'m', 'LineWidth',2);

end