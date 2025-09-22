function [dV2 , theta3, theta4] = CambioAnPericentro( a_i, e_i, om_i, om_f,mu)

% ruota l'orbita sullo stesso piano senza modificare gli altri parametri geometrici
% ci sono due punti in cui è possibile effettuare la manovra : punto A , punto B
%
% Input
%
%---------------------------------------------------------------------------
% a_i   [1x1] initial semimajor axis     [Km]
% e_i   [1x1] initial eccentricity     
% om_i  [1x1] initial pericenter anomaly [rad]
% om_f  [1x1] final pericenter anomaly [rad]
% mu    [1x1] gravitational parameter


% s = scelta 
% s = A --> manovra in A
% s = B --> manovra in B

% nota : vedi lab 2 se om_f = om finale oppure om_f = om finale + pi

delta_om = om_f - om_i;

s = input(' punto di manovra cambio om (A:1  o B:2):');

switch s
    case 1
        theta3 = (delta_om) / 2 ;
    case 2
        theta3 = pi + (delta_om)/2 ;
end


    theta4 = 2*pi - theta3;

% trovo il costo di manovra

Vr = sqrt (  mu / ( a_i * ( 1 - e_i^ 2)) ) * e_i * sin (( om_f - om_i)/2);
dV2= 2*abs(Vr);

end