
close all
clear all
clc


%% dati

mu = 398600;

% PUNTO INIZIALE

rr_i = [ -4091.2866 , 6563.7099 , 3713.1316 ];
vv_i = [ -6.0730 , -3.8110 , 0.1904 ];

% PUNTO FINALE

a_f = 13960;
e_f = 0.3339;
i_f = 0.7923;
OM_f = 2.1220;
om_f = 1.0030;
theta_f = 1.3000;
[rr_f,vv_f] = parorb2rv(a_f,e_f,i_f,OM_f,om_f,theta_f,mu);


% conversione da rr vv a parametri orbitali

[a_i,e_i,i_i,OM_i,om_i,theta_i] = rv2parorb(rr_i,vv_i, mu);

% ricavo u1 e u2 dal cambio piano da i_i OM_i a i_f OM_f dell'orbita 

[u1,u2,dV_prova,om_prova,theta_prova] = CambioPiano(a_f,e_f,i_i,OM_i,om_i,i_f,OM_f,mu);

% plot orbita iniziale e finale

  theta_oi = 0;    
  theta_of = 2*pi;

  Terra_3D(6371);
  hold on
  
  plotOrbit(a_i,e_i,i_i,OM_i,om_i,theta_oi,theta_of,0.01,mu,'m')

  plotOrbit(a_f,e_f,i_f,OM_f,om_f,theta_oi,theta_of,0.01,mu,'k')

  plot3(rr_i(1),rr_i(2),rr_i(3),'xg','LineWidth',2)
 
  plot3(rr_f(1),rr_f(2),rr_f(3),'xr','LineWidth',2)


%% STRATEGIA ALTERNATIVA 1 
%1- CAMBIO ANOMALIA NEL PRIMO PUNTO DISPONIBILE (A)  (1)
%2- CAMBIO FORMA CON MANOVRA PERICENTRO -> APOCENTRO (1)
%3- CAMBIO PIANO CON MANOVRA NEL PUNTO THETA + PI    (2)

%% 1
% cambio anomalia del pericentro per portare om_i --> om_t = om_f - d_om_parassita = om_f - (u2-u1)

om_t = om_f + u1 - u2;  

[dV0 , theta_01, theta_02] = CambioAnPericentro( a_i , e_i, om_i, om_t,mu);

% tempo per arrivare al punto di manovra cambio anomalia
deltat_0 = TOF(a_i, e_i, theta_i, theta_01, mu);

plotOrbit_thic(a_i, e_i, i_i, OM_i, om_i, theta_i, theta_01, 0.01, mu);

%% 2
% cambio forma 

s = input ( ' tipologia di manovra cambio forma ( 1:P-A  2:A-P 3:P-P  4:A-A ) : ');

[dV1 , dV2, dT2, a_t, e_t, theta_man, theta_arr] = CambioFormaOrbita(a_i, e_i, a_f, e_f, mu, om_t, i_i, OM_i, s);

% tempo per arrivare al punto di manovra cambio forma
plotOrbit_thic(a_i,e_i,i_i,OM_i,om_t,theta_02,theta_man,0.01,mu)

deltat_1 = TOF(a_i,e_i,theta_02,theta_man,mu);

% tempo sull'orbita di trasferimento 
deltat_2 = dT2;

%% 3 
% cambio piano

[u1,u2,dV3,om_3,theta_3] = CambioPiano(a_f,e_f,i_i,OM_i,om_t,i_f,OM_f,mu);


% tempo per arrivare al punto di manovra cambio piano
deltat_3 = TOF(a_f,e_f,theta_arr,theta_3,mu);

plotOrbit_thic(a_f,e_f,i_i,OM_i,om_t,theta_arr,theta_3,0.01,mu)

% tempo trascorso sull'orbita finale per arrivare al puntp finale
deltat_5 = TOF(a_f, e_f, theta_3, theta_f, mu);

plotOrbit_thic(a_f,e_f,i_f,OM_f,om_f,theta_3,theta_f,0.01,mu)

%% CONCLUSIONE

% costo tempo totale
deltat_TOT = deltat_0 + deltat_1 + deltat_2 + deltat_3+deltat_5
% costo energia totale
dVTOT = dV0+dV1 + dV2 + dV3

legend( '', ' Initial orbit ', 'Final orbit', 'Initial point', 'Final point', 'Initial Orbit','Transfer 2',' Transfer 1', 'Transfer 3', 'Final Orbit');

%% se necessari : plot completi delle orbite intermedie

%plotOrbit(a_f,e_f,i_i,OM_i,om_i,theta_oi,theta_of,0.01,mu)
%plotOrbit(a_f,e_f,i_i,OM_i,om_i+pi,theta_oi,theta_of,0.01,mu)
%plotOrbit(a_f,e_f,i_f,OM_f,om_3,theta_oi,theta_of,0.01,mu)
%plotOrbit(a_f,e_f,i_f,OM_f,om_f,theta_oi,theta_of,0.01,mu)



