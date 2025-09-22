close all
clear all
clc

% dati

mu = 398600;

% PUNTO INIZIALE

rr_i = [ -4091.2866 , 6563.7099 , 3713.1316 ];
vv_i = [ -6.0730 , -3.8110 , 0.1904 ];

% PUNTO FINALE

a_f = 13960;
E_f = 0.3339;
i_f = 0.7923;
OM_f = 2.1220;
om_f = 1.0030;
theta_f = 1.3000;
[rr_f,vv_f] = parorb2rv(a_f,E_f,i_f,OM_f,om_f,theta_f,mu);


% conversione da rr vv a parametri orbitali

[a_i,E_i,i_i,OM_i,om_i,theta_i] = rv2parorb(rr_i,vv_i, mu);

% plot orbita iniziale e finale
  theta_oi = 0;    
  theta_of = 2*pi;

  Terra_3D(6371);
  hold on
  
  plotOrbit(a_i,E_i,i_i,OM_i,om_i,theta_oi,theta_of,0.01,mu,'k-')

  plotOrbit(a_f,E_f,i_f,OM_f,om_f,theta_oi,theta_of,0.01,mu,'m-')

  plot3(rr_i(1),rr_i(2),rr_i(3),'xg','LineWidth',2)
 
  plot3(rr_f(1),rr_f(2),rr_f(3),'xr','LineWidth',2)


%% procedura standard CAMBIO FORMA -> CAMBIO PIANO -> CAMBIO ANOMALIA
%1- CAMBIO FORMA CON MANOVRA PERICENTRO-> APOCENTRO  (1)
%2- CAMBIO PIANO CON MANOVRA NEL PUNTO THETA + PI    (2)
%3- CAMBIO ANOMALIA CON MANOVRA NEL PUNTO A          (1)

%% 1
% cambio forma 

s = input ( ' tipologia di manovra cambio forma ( 1:P-A  2:A-P 3:P-P  4:A-A ) : ');

[dV1 , dV2, dT2, a_t, e_t, theta_man, theta_arr] = CambioFormaOrbita(a_i, E_i, a_f, E_f, mu, om_i, i_i, OM_i, s);

% tempo per arrivare al punto di manovra cambio forma
deltat_1 = TOF(a_i,E_i,theta_i,theta_man,mu);

plotOrbit_thic(a_i,E_i,i_i,OM_i,om_i,theta_i,theta_man,0.01,mu)

% tempo sull'orbita di trasferimento 
deltat_2 = dT2;

if s==1 || s==2
disp( ' la manovra scelta è A->P o P-> A  :  om = om_i ');
% anomalia pericentro nuova orbita 
om= om_i;


end

if s==3 || s==4

disp( ' la manovra scelta è P->P o A-> A  :  om = om_i + pi');
% anomalia pericentro nuova orbita
om = om_i + pi;

end

%% 2
% cambio piano

[u1,u2,dV3,om_3,theta_3] = CambioPiano(a_f,E_f,i_i,OM_i,om,i_f,OM_f,mu);

% tempo per arrivare al punto di manovra cambio piano
deltat_3 = TOF(a_f,E_f,theta_arr,theta_3,mu);

plotOrbit_thic(a_f,E_f,i_i,OM_i,om,theta_arr,theta_3,0.01,mu)


%% 3
% cambio anomalia del pericentro

[dV4 , theta_4, theta_5] = CambioAnPericentro( a_f , E_f, om_3, om_f,mu);

% tempo per arrivare al punto di manovra cambio anomalia
deltat_4 = TOF(a_f, E_f, theta_3, theta_4, mu);

plotOrbit_thic(a_f,E_f,i_f,OM_f,om_3,theta_3,theta_4,0.01,mu)


% tempo per arrivare al punto finale
deltat_5 = TOF(a_f, E_f, theta_5, theta_f, mu);

plotOrbit_thic(a_f,E_f,i_f,OM_f,om_f,theta_5,theta_f,0.01,mu)


%% CONCLUSIONE

% costo temporale totale
deltat_TOT = deltat_1 + deltat_2 + deltat_3 + deltat_4 + deltat_5
% costo energetico totale
dVTOT = dV1 + dV2 + dV3 + dV4


legend( '', ' Initial orbit ', 'Final orbit', 'Initial point', 'Final point', 'Transfer 2', ' Transfer 1', 'Transfer 3', 'Transfer 4', ' Final orbit');

%% se necessari : plot completi delle orbite intermedie

%plotOrbit(a_f,E_f,i_i,OM_i,om_i,theta_oi,theta_of,0.01,mu)
%plotOrbit(a_f,E_f,i_i,OM_i,om_i+pi,theta_oi,theta_of,0.01,mu)
%plotOrbit(a_f,E_f,i_f,OM_f,om_3,theta_oi,theta_of,0.01,mu)
%plotOrbit(a_f,E_f,i_f,OM_f,om_f,theta_oi,theta_of,0.01,mu)




