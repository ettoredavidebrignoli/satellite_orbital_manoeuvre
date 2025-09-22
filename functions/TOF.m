function deltat = TOF(a,e,th1,th2,mu)

% Time of flight
%
% deltat = TOF( a,e,th1,th2,mu)
%
%-----------------------------------------------------------------------
% Input arguments:
% a     [1x1]   semi-major axis             [km]
% e     [1x1]   eccentricity                [-]
% th1   [1x1]   initial true anomaly        [rad]
% th2   [1x1]   final true anomaly          [rad]
% mu    [1x1]   gravitational parameter     [km^3/s^2]
%
%-----------------------------------------------------------------------
%
% Output argument:
% deltat [1x1] time of flight [s]
%
%-----------------------------------------------------------------------


if (th1 > 2*pi)
    th1= th1 - 2*pi;
end
if (th2 > 2*pi)
    th2= th2 - 2*pi;
end
if (th1 < 0 )
    th1 = th1 + 2*pi;
end
if (th2 < 0)
    th2 = th2 +2*pi;
end

E1 = 2*atan( sqrt ( (1-e)/(1+e)) * tan (th1 /2) );
E2 = 2*atan( sqrt ( (1-e)/(1+e)) * tan (th2 /2) );

if E1 <0
    E1 = 2*pi +E1;
end

if E2 <0
    E2 = 2*pi +E2;
end

deltat = sqrt( a^3 / mu)*( (E2 - E1 ) - e*( sin(E2) - sin(E1)));

if th2 < th1 
    deltat = deltat + ( 2*pi * sqrt ( a^3 / mu));
end

end




