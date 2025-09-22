function [dV3 , dV4, dT4, a_t, e_t, theta_man, theta_arr] = CambioFormaOrbita(a_i, e_i, a_f, e_f, mu, om, i, OM, s)

% funzione per cambiare la forma dell'orbita attraverso una manovra bitangente
%
% input
%
%-------------------------------------------------------------------------
% a_i   [1x1] initial semimajor axis     [Km]
% e_i   [1x1] initial eccentricity 
% a_i   [1x1] final semimajor axis       [Km]
% e_i   [1x1] final eccentricity 
% mu    [1x1] gravitational parameter
% per il plot dell'orbita di trasferimento:
% om    [1x1] initial pericenter anomaly 
% i     [1x1] inclination                [rad]
% OM    [1x1] RAAN                       [rad]
% s = scelta della manovra 
% 1: pericentro --> apocentro
% 2: apocentro --> pericentro
% 3: apocentro --> apocentro
% 4: pericentro --> pericentro
%
% Output
%
%-------------------------------------------------------------------------
% dV3   [1x1] costo primo impulso                           [Km/s]
% dV4   [1x1] costo secondo impulso                         [Km/s]
% dT4   [1x1] tempo trascorso sull'orbita di trasferimento  [s]
% a_t   [1x1] transfer orbit semimajor axis                 [Km]
% e_t   [1x1] transfer orbit eccentricity 
% theta_man [1x1] anomalia vera prima manovra sull'orbita iniziale [rad]
% theta_arr [1x1] anomalia vera seconda manovra sull'orbita finale [rad]



rp_i = a_i * (1 - e_i);
ra_i = a_i * (1 + e_i);
rp_f = a_f * (1 - e_f);
ra_f = a_f * (1 + e_f);
vp_i = sqrt( 2 * mu * ( (1/rp_i) - (1/(2*a_i))) );
va_i = sqrt( 2 * mu * ( (1/ra_i) - (1/(2*a_i))) );
vp_f = sqrt( 2 * mu * ( (1/rp_f) - (1/(2*a_f))) );
va_f = sqrt( 2 * mu * ( (1/ra_f) - (1/(2*a_f))) );


switch s

    case 1  % P->A
        rp_t = rp_i;
        ra_t = ra_f;
        a_t= (rp_t + ra_t) / 2;
        e_t = (ra_t - rp_t)/(ra_t + rp_t);
        vp_t = sqrt( 2*mu*( (1/rp_t) - (1 / (2*a_t))));
        va_t = sqrt( 2*mu*( (1/ra_t) - (1 / (2*a_t))));
        dV3 = abs ( vp_i - vp_t);
        dV4 = abs ( va_t - va_f);
        dT4 = pi * sqrt( a_t^3 / mu );
        theta_man = 0;
        theta_arr = pi;
        om_t = om;
        theta_oi = 0;
        theta_of = pi;

    

    case 2  % A->P
        rp_t = rp_f;
        ra_t = ra_i;
        a_t= (rp_t + ra_t) / 2;
        e_t = (ra_t - rp_t)/(ra_t + rp_t);
        vp_t = sqrt( 2*mu*( (1/rp_t) - (1 / (2*a_t))));
        va_t = sqrt( 2*mu*( (1/ra_t) - (1 / (2*a_t))));
        dV3 = abs ( va_i - va_t);
        dV4 = abs ( vp_t - vp_f);
        dT4 = pi * sqrt( a_t^3 / mu );
        theta_man = pi;
        theta_arr= 0;
        % nel nostro caso rp_f < ra_i
        om_t=om;
        theta_oi = pi;
        theta_of = 2*pi;
   
   

    case 3 % P->P
        rp_t = rp_i;
        ra_t = rp_f;
        a_t= (rp_t + ra_t) / 2;
        e_t = (ra_t - rp_t)/(ra_t + rp_t);
        vp_t = sqrt( 2*mu*( (1/rp_t) - (1 / (2*a_t))));
        va_t = sqrt( 2*mu*( (1/ra_t) - (1 / (2*a_t))));
        dV3 = abs ( vp_i - vp_t);
        dV4 = abs ( va_t - vp_f);
        dT4 = pi * sqrt( a_t^3 / mu );
        theta_man = 0;
        theta_arr = 0;
        % considerando che l'orbita iniziale nel nostro problema è piu
        % piccola dell'orbita finale 
        om_t = om;
        theta_oi= 0;
        theta_of= pi;
         
    case 4 % A->A
        rp_t = ra_i;
        ra_t = ra_f;
        a_t= (rp_t + ra_t) / 2;
        e_t = (ra_t - rp_t)/(ra_t + rp_t);
        vp_t = sqrt( 2*mu*( (1/rp_t) - (1 / (2*a_t))));
        va_t = sqrt( 2*mu*( (1/ra_t) - (1 / (2*a_t))));
        dV3 = abs ( va_i - vp_t);
        dV4 = abs ( va_t - va_f);
        dT4 = pi * sqrt( a_t^3 / mu );
        theta_man = pi;
        theta_arr= pi;
        om_t = om+pi;
        theta_oi= 0;
        theta_of= pi;
end


 
 plotOrbit_thic(a_t,e_t,i,OM,om_t,theta_oi,theta_of,0.01,398600)


end

