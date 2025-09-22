close all
clear all
clc

% dati

mu = 398600;
k = [0 0 1];
dt = 0.01;

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


ra_i = a_i * ( 1 + E_i);
rp_i = a_i * ( 1 - E_i);
p_i =  a_i * (1-E_i^2);
ra_f = a_f * ( 1 + E_f);
rp_f = a_f * ( 1 - E_f);
p_f =  a_f * (1 - E_f^2);


% plot orbita iniziale e finale
  theta_oi = 0;    
  theta_of = 2*pi;

  
  %% TROVO IL PIANO DELL'ORBITA DI TRASFERIMENTO

  % vettore normale al piano su cui giace l'orbita di trasferimento

  h_t = (cross(rr_i,rr_f) / norm(cross(rr_i,rr_f)));

  % inclinazione dell'orbita di trasferimento
  i_t = acos(h_t(3));

  % RAAN dell'orbita di trasferimento
  n_t = ( cross(k,h_t)/ norm(cross(k,h_t)));

  if n_t(2)>=0

      OM_t = acos(n_t(1));
  else
   
      OM_t = 2*pi - acos(n_t(1));
  end

  % la differenza di anomalia tra il punto iniziale e il punto finale
  % sull'orbita di trasferimento e':

  delta_theta = acos( dot(rr_i,rr_f) / (norm(rr_i)*norm(rr_f)));

  % trovo i raggi del punto iniziale e finale

  r_f = ( a_f*(1-E_f^2))/(1+E_f*cos(theta_f));
  r_i = ( a_i*(1-E_i^2))/(1+E_i*cos(theta_i));



  %% faccio variare l'eccentricità come parametro dal primo risultato reale fino a 1 (escluso)
  
  % per trovare om_t di ciascuna orbita trovo alpha 
   alpha = acos(dot(rr_i,n_t)/norm(rr_i));
   %{
    if rr_i(3) >= 0
      om_t = alpha - theta_ti
      else
      om_t = 2*pi - alpha - theta_ti
      end
   %}
   % dove theta_ti è l'anomalia vera del punto iniziale sull'orbita di trasferimento 

   E_t = 0.16:0.0001:0.17;

   LUN= 100;

for x= 1: LUN

   % ricavo theta_ti
   syms TH_1
   solve ( (r_i / r_f ) == ( ( 1 + E_t(x)*cos(TH_1 + delta_theta)) / (1 + E_t(x)*cos(TH_1) )), TH_1 );
   th_1 = double(ans);

   % trovo due valori di theta_ti per ciascun valore di eccentricità

   theta_t1 =  real(th_1(1));
   theta_t2 =  real(th_1(2));

   if( theta_t1 < 0 )
       theta_t1 = 2*pi + theta_t1;
   end
   if(theta_t2 < 0 )
       theta_t2 = 2*pi + theta_t2;
   end

   if(rr_i(3)>=0)

      om_t1 = alpha - theta_t1;
      om_t2 = alpha - theta_t2;

  else
      om_t1 = 2*pi - alpha - theta_t1;
      om_t2 = 2*pi - alpha - theta_t2;
   end
   
   % semiasse maggiore
   a_t1 = ( r_i * (1 + E_t(x) *cos(theta_t1)))/ (1 - E_t(x)^2);
   a_t2 = ( r_i * (1 + E_t(x) *cos(theta_t2)))/ (1 - E_t(x)^2);

   % dV manovra nel caso theta_ti = theta_t1

   [rr_a1,vv_a1] = parorb2rv(a_t1,E_t(x),i_t,OM_t,om_t1,theta_t1,mu);
   [rr_b1,vv_b1] = parorb2rv(a_t1,E_t(x),i_t,OM_t,om_t1,theta_t1 + delta_theta,mu);
   beta1 = acos( ( dot(vv_i,vv_a1)) / ( ( norm(vv_i )*norm(vv_a1) ) ) );
   delta1 = acos( ( dot(vv_b1,vv_f)) / ( ( norm(vv_b1 )*norm(vv_f) ) ) );
   % dV orbita iniziale -> orbita intermedia con formula di Crnot
   dVa1 = sqrt ( norm(vv_i)^2 + norm(vv_a1)^2 - 2* norm(vv_i)*norm(vv_a1)*cos(beta1));
   % dV orbita intermedia -> orbita finale
   dVb1 = sqrt ( norm(vv_f)^2 + norm(vv_b1)^2 - 2* norm(vv_f)*norm(vv_b1)*cos(delta1));
   % dV totale
   dVTOT1 = dVa1 + dVb1;


   % dV manovra nel caso theta_ti = theta_t1
   [rr_a2,vv_a2] = parorb2rv(a_t2,E_t(x),i_t,OM_t,om_t2,theta_t2,mu);
   [rr_b2,vv_b2] = parorb2rv(a_t2,E_t(x),i_t,OM_t,om_t2,theta_t2 + delta_theta,mu);
   beta2 = acos( ( dot(vv_i,vv_a2)) / ( ( norm(vv_i )*norm(vv_a2) ) ) );
   delta2 = acos( ( dot(vv_b2,vv_f)) / ( ( norm(vv_b2 )*norm(vv_f) ) ) );
   % dV orbita iniziale -> orbita intermedia con formula di Carnot
   dVa2 = sqrt ( norm(vv_i)^2 + norm(vv_a2)^2 - 2* norm(vv_i)*norm(vv_a2)*cos(beta2));
   % dV orbita intermedia -> orbita finale
   dVb2 = sqrt ( norm(vv_f)^2 + norm(vv_b2)^2 - 2* norm(vv_f)*norm(vv_b2)*cos(delta2));
   % dV totale
   dVTOT2 = dVa2 + dVb2;
   
   % vettore costo temporale ed energetivo per theta_t1, theta_t2

   delta_t1(x) = TOF(a_t1,E_t(x),theta_t1,theta_t1 + delta_theta, mu);
   delta_t2(x) = TOF(a_t2,E_t(x),theta_t2,theta_t2 + delta_theta, mu);
   delta_v1(x) = dVTOT1;
   delta_v2(x) = dVTOT2;

%{ 
plot di alcune orbite arbitrarie per osservare l'andamento
if (x == 1 || x == 4 || x == 7) && delta_t1(x)< 15000
   plotOrbit(a_t1, E_t(x), i_t, OM_t, om_t1, theta_t1, theta_t1 + delta_theta, 0.01, mu, 'g');
end
if (x == 1 || x == 4 || x == 7) && delta_t2(x)< 15000
   plotOrbit(a_t2, E_t(x), i_t, OM_t, om_t2, theta_t2, theta_t2 + delta_theta, 0.01, mu, 'g'); 
end
%}


   end

% riordino i vettori costi temporale ed energetico per un plot più ordinato

   for i= 1:LUN
     if delta_t1(i) < delta_t2(i)
         temp = delta_t1(i);
         delta_t1(i)=delta_t2(i);
         delta_t2(i)= temp;
     end
     if delta_v1(i) < delta_v2(i)
         temp = delta_v1(i);
         delta_v1(i)=delta_v2(i);
         delta_v2(i)= temp;
     end
     
   end

%% trasferimento con orbita parabolica ( e_t = 1)

% trovo theta_ti

syms th_t

solve (  (r_i/r_f) == ( (1+cos(th_t + delta_theta) ) / ( 1 + cos(th_t ))) , th_t );

TH_t1 = double(ans(1));
TH_t2 = double(ans(2));

% trovo i due valori di theta_ti per e_t = 1
theta_t1 = real(TH_t1);
theta_t2 =  real(TH_t2);

% trovo anomalia del pericentro per entrambi i casi

if(rr_i(3)>=0)
      om_t1 = alpha - theta_t1;
      om_t2 = alpha - theta_t2;

  else
      om_t1 = 2*pi - alpha - theta_t1;
      om_t2 = 2*pi - alpha - theta_t2;
end


% trovo il semilato retto 

rp_1 = ( r_i *  ( 1 + cos(theta_t1) ) ) /2;
rp_2 = ( r_i *  ( 1 + cos(theta_t2) ) ) /2;


hold on



% usiamo la parabola definita da theta_t1 , om_t1, rp_1 

% calcolo costo temporale

D_1 = tan(theta_t1 /2);
t_1 = 0.5 * ( sqrt(( 2*rp_1)^3 / mu)* ( D_1 + ( D_1 ^ 3 / 3 )));

D_2 = tan(( theta_t1 + delta_theta )  /2);
t_2 = 0.5 * ( sqrt(( 2*rp_1)^3 / mu)* ( D_2 + ( D_2 ^ 3 / 3 )));

delta_t_parabola = t_2 - t_1

% calcolo costo energetico

[rr_p1,vv_p1] = parorb2rv_parab(rp_1,i_t,OM_t,om_t1,theta_t1,mu);
[rr_p2,vv_p2] = parorb2rv_parab(rp_1,i_t,OM_t,om_t1,theta_t1 + delta_theta,mu);

beta = acos( ( dot(vv_p1,vv_i)) / ( ( norm(vv_p1 )*norm(vv_i) ) ) );
delta = acos( ( dot(vv_p2,vv_f)) / ( ( norm(vv_p2 )*norm(vv_f) ) ) );
dV1_parabola = sqrt ( norm(vv_i)^2 + norm(vv_p1)^2 - 2* norm(vv_i)*norm(vv_p1)*cos(beta));
dV2_parabola = sqrt ( norm(vv_f)^2 + norm(vv_p2)^2 - 2* norm(vv_f)*norm(vv_p2)*cos(delta));

dVTOT_parabola = dV1_parabola + dV2_parabola


   %% plotto l'orbita con dV minimo 

   % abbiamo individuato grafucamente la zona in cui è compreso il minimo
   % del costo energetico per la manovra di trasferimento diretto,
   % aumentando la precisione di e_t a 0.0001 abbiamo ricavato il valore di
   % eccentricità a cui corrisponde il minimo dV

   % con lo stesso procedimento di prima ricavo tutti i parametri corrispondenti a e_t2 = 0.2391

   e_t2 = 0.2391;

   % trovo theta_ti

   syms TH_1
   solve ( (r_i / r_f ) == ( ( 1 + e_t2*cos(TH_1 + delta_theta)) / (1 + e_t2*cos(TH_1) )), TH_1 );
   th_1 = double(ans);

   theta_t1 =  real(th_1(1));
   theta_t2 =  real(th_1(2));

   if( theta_t1 < 0 )
       theta_t1 = 2*pi + theta_t1;
   end
   if(theta_t2 < 0 )
       theta_t2 = 2*pi + theta_t2;
   end

   % anomalia del pericentro 

   if(rr_i(3)>=0)

      om_t1 = alpha - theta_t1;
      om_t2 = alpha - theta_t2;

  else
      om_t1 = 2*pi - alpha - theta_t1;
      om_t2 = 2*pi - alpha - theta_t2;
   end
   
   % semiasse maggiore
   a_t1 = ( r_i * (1 + e_t2 *cos(theta_t1)))/ (1 - e_t2^2);
   a_t2 = ( r_i * (1 + e_t2 *cos(theta_t2)))/ (1 - e_t2^2);

   % trovo il costo energetico totale per entrambi i casi poi seleziono il
   % caso con dV minimo

   [rr_a1,vv_a1] = parorb2rv(a_t1,e_t2,i_t,OM_t,om_t1,theta_t1,mu);
   [rr_b1,vv_b1] = parorb2rv(a_t1,e_t2,i_t,OM_t,om_t1,theta_t1 + delta_theta,mu);
   beta1 = acos( ( dot(vv_i,vv_a1)) / ( ( norm(vv_i )*norm(vv_a1) ) ) );
   dVa1 = sqrt ( norm(vv_i)^2 + norm(vv_a1)^2 - 2* norm(vv_i)*norm(vv_a1)*cos(beta1));
   delta1 = acos( ( dot(vv_b1,vv_f)) / ( ( norm(vv_b1 )*norm(vv_f) ) ) );
   dVb1 = sqrt ( norm(vv_f)^2 + norm(vv_b1)^2 - 2* norm(vv_f)*norm(vv_b1)*cos(delta1));
   dVTOT1 = dVa1 + dVb1;
   [rr_a2,vv_a2] = parorb2rv(a_t2,e_t2,i_t,OM_t,om_t2,theta_t2,mu);
   [rr_b2,vv_b2] = parorb2rv(a_t2,e_t2,i_t,OM_t,om_t2,theta_t2 + delta_theta,mu);
   beta2 = acos( ( dot(vv_i,vv_a2)) / ( ( norm(vv_i )*norm(vv_a2) ) ) );
   dVa2 = sqrt ( norm(vv_i)^2 + norm(vv_a2)^2 - 2* norm(vv_i)*norm(vv_a2)*cos(beta2));
   delta2 = acos( ( dot(vv_b2,vv_f)) / ( ( norm(vv_b2 )*norm(vv_f) ) ) );
   dVb2 = sqrt ( norm(vv_f)^2 + norm(vv_b2)^2 - 2* norm(vv_f)*norm(vv_b2)*cos(delta2));
   dVTOT2 = dVa2 + dVb2;
   
   if dVTOT1 < dVTOT2
   a_t_dvmin = a_t1;
   theta_i_dvmin = theta_t1;
   theta_f_dvmin = theta_t1 + delta_theta;
   dVa1_mintrasferimentodiretto = dVa1  ;
   dVb1_mintrasferimentodiretto = dVb1;
   dV1_min_trasferimentodiretto = dVTOT1
   plotOrbit_thic(a_t1, e_t2, i_t, OM_t, om_t1, theta_t1, theta_t1 + delta_theta, 0.01, mu);
   delta_t_dvmin = TOF(a_t1,e_t2,theta_t1,theta_t1 + delta_theta, mu)
   else
   a_t_dvmin = a_t2
   theta_i_dvmin = theta_t2;
   theta_f_dvmin = theta_t2 + delta_theta;
   dVa2_mintrasferimentodiretto = dVa2 ;
   dVb2_mintrasferimentodiretto = dVb2;
   dV2_min_trasferimentodiretto = dVTOT2
   plotOrbit_thic(a_t2, e_t2, i_t, OM_t, om_t2, theta_t2, theta_t2 + delta_theta, 0.01, mu);
   delta_t_dvmin = TOF(a_t2,e_t2,theta_t2,theta_t2 + delta_theta, mu)
   end 











%% plot dell'andamento di dV e dT al variare dell'eccentricità 

figure

 plot(E_t, delta_t1,'k', E_t, delta_t2, 'k', 1, delta_t_parabola, 'xr');
 ylim([0 10000])

legend(' delta t elliptic transfer ', ' ','delta t parabolic transfer')

figure

 plot(E_t, delta_v1,'k', E_t, delta_v2, 'k', 1, dVTOT_parabola, 'xr');
 legend(' delta v elliptic transfer ', ' ','delta v parabolic transfer')


