function [u1,u2,dV1,om_2,theta2] = CambioPiano(a_i,e_i,i_i,OM_i,om_i,i_f,OM_f,mu)

% funzione che cambia il piano orbitale
% cambio di inclinazione e di ascensione retta {i_i ,   OM_i } --> { i_f , OM_f }
% mantenendo invariati { a , e }
% e variando anomalia del pericentro come conseguenza { om_i } --> { om_f }

% input
%
%--------------------------------------------------------------------------
% a_i   [1x1] initial semimajor axis     [Km]
% e_i   [1x1] initial eccentricity     
% i_i   [1x1] initial inclination        [rad]
% OM_i  [1x1] initial RAAN               [rad]
% om_i  [1x1] initial pericenter anomaly [rad]
% i_f   [1x1] final inclination          [rad]
% OM_f  [1x1] final RAAN                 [rad]
% mu    [1x1] gravitational parameter    [Km^3/s^2]
%
% output
%
%------------------------------------------------------------------------
% dV1    [1x1]  costo manovra              [Km/s]
% om_2   [1x1]  final pericenter anomaly   [rad]
% theta2 [1x1]  anomalia vera di manovra   [rad]



% s = scelta del punto di manovra ( puo essere effettuato sia in theta che
% in theta + pi )

s = input(' punto di manovra cambio piano ( 1:theta 2:theta + pi ): ');

if  OM_i == OM_f    %caso di rotazione di h attorno ad N 
    om_2 = om_i;
    theta2= 2*pi - om_i;
    alpha= (i_1 - i_2);
else

    % trovo alpha tra le velocità 

    alpha = acos( sin(i_i)*sin(i_f)*cos(OM_f-OM_i)+cos(i_i)*cos(i_f));

    % trovo u1 

    cos_u1 = ( ( cos(i_i)*cos(alpha)) - cos(i_f) ) / ( sin(i_i)*sin(alpha)) ;

    sin_u1 = ( sin (OM_f-OM_i)* sin(i_f)) / (sin(alpha)) ;

     % trovo u2

    cos_u2 = ( -( cos(i_f)*cos(alpha))  + cos(i_i) ) / ( sin(i_f)*sin(alpha)) ;
 
    sin_u2 = ( sin (OM_f-OM_i)* sin(i_i)) / (sin(alpha)) ;


    if( i_f - i_i) < 0
        cos_u1 = -cos_u1;
        cos_u2 = -cos_u2;
    end

    u2 = atan2( sin_u2 , cos_u2);
    u1 = atan2( sin_u1 , cos_u1);

    % trovo la posizione di manovra theta2 e l'anomalia al pericentro om_2

    if( ( i_f - i_i )* (OM_f - OM_i) > 0)   %delta I delta OM segno concorde

        theta2= u1 - om_i;
        om_2 =  u2 - theta2;
    else
        theta2 = 2*pi - om_i - u1;
        om_2   = 2*pi - u2 - theta2;

    end
    
if s == 2
    theta2= theta2 + pi;
   
end



if theta2 < 0 
    theta2 = 2*pi+theta2;
end

    %trovo il costo della manovra 
    V_t = sqrt( mu / (a_i * (1- e_i^2) ) )* ( 1 + e_i*cos(theta2));
    dV1= 2 * V_t * sin( (alpha)/2);

end