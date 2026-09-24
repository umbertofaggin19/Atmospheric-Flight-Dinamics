%6DOF in sistema di riferimento BODY
%per semplificare suppongo cI = 0 , cN = 0

%SI NOTI CHE IL MODO DI SPIRALE è INSTABILE (LENTA, ~ 150sec) , PERTANTO VA CORRETTA NEL TEMPO
%CON INPUT DEL PILOTA, da tenere conto nella formula di cI. (mode initiated
%by a rolling r heading disturbance, i.e. p or r =! 0)


%instead fugoid, short mode and dutch mode seems stable


%roll damping should be stable but seems unstable. in realtà roll damping
%funziona bene, è solo che eccitando p suscito sia roll damping che spiral
%mode. Roll damping funziona bene (p smorzato nei primi 5s), poi subentra
%spiral mode, quindi tutto ok



clc;
close all;
clear;

%importo dati di xflr5
T = readtable('dati23012.csv');
alpha_tab = T{5:65,1};
cL_tab = T{5:65,3};
cM_tab = T{5:65,9};
cD_tab = T{5:65,6};
cI_tab = T{5:65,8};
cN_tab = T{5:65,10};
cY_tab = T{5:65,7};


%PARAMETRI
m = 0.75;
g = 9.81;
b = 1.4;
S = 0.205;
Ar = b^2/S;
e = 0.75;
cD0 = 0.02;
k = 1/(e*pi*Ar);

%(dati relativi all'inerzia, supp Ixz = 0)
Ixx = 0.0481;
Iyy = 0.03591;
Izz = 0.08375;
Ixz = 0.11;
c_media = 0.151;

rho0 = 1.225;
z0 = 11000;
Z0 = -z0;

%STATO INIZIALE
alpha = deg2rad(0.74891 + 3);
beta = deg2rad(0);

phi = deg2rad(0);     %rollio
theta = alpha;   %beccheggio
psi = deg2rad(0);     %imbardata

p = deg2rad(-2);    %roll rate
q = deg2rad(0);    %pitch rate
r = deg2rad(0);    % yaw rate

p_dot = deg2rad(1);    
q_dot = deg2rad(0);    
r_dot = deg2rad(0);

V = 25;
u = V*cos(alpha);
v = 0;
w = V*sin(alpha);


x_ned = 0;
y_ned = 0;
z_ned = 0;
Z = 0;



dt = 0.001;
tmax = 100;
time = 0:dt:tmax;

errore_cumulato0 = 0;
V_old = 0;

for i = 1:length(time)

    t = time(i);
    time_eff(i) = t;

    rho = rho0*exp(z_ned/z0);
    rho_eff(i) = rho;

    qp = 0.5*rho*(V^2);
    

    alpha = atan2(w,u);
    alpha_deg = rad2deg(alpha);
    alpha_eff(i) = alpha_deg;
    beta = asin(v/V);  
    beta_eff(i) = rad2deg(beta);

    cL = interp1(alpha_tab,cL_tab,alpha_deg,'linear','extrap'); 
    cL_eff(i) = cL;
    
    if cL > 2 || isnan(cL)
        fprintf('stallo %.2f\n', time(i));
        break
    end

    cD = interp1(alpha_tab,cD_tab,alpha_deg,'linear','extrap');
    cD_eff(i) = cD;




    cY_beta = interp1(alpha_tab,cY_tab,alpha_deg,'linear','extrap');
    cY_p = 0.02272;
    cY_r = 0.1395;


    cY = cY_beta*beta + cY_p * (p*b)/(2*V) + cY_r * (r*b)/(2*V);

    Tx = 1.0153;
    Ty = 0;
    Tz = 0; 
    T = sqrt((Tx)^2 + (Ty)^2 + (Tz)^2);
    L = 0.5*rho*(V^2)*S*cL;
    D = 0.5*rho*(V^2)*S*cD;
    S_lat = qp*S*cY;
    
    Xa = - D*cos(alpha)*cos(beta) - S_lat*cos(alpha)*sin(beta) + L*sin(alpha);
    Ya = - D*sin(beta) + S_lat*cos(beta);
    Za = - D*sin(alpha)*cos(beta) - S_lat*sin(alpha)*sin(beta) -L*cos(alpha);
   
    %1 equazione Newton (TRASLAZIONE BARICENTRO)

    fx = -m*g*sin(theta) + Tx + Xa;
    fy = m*g*sin(phi)*cos(theta) + Ty + Ya;
    fz = m*g*cos(phi)*cos(theta) + Tz + Za;


    u_dot = (r*v - q*w) + fx/m;
    v_dot = (p*w - u*r) + fy/m;
    w_dot = (q*u - p*v) + fz/m;


    %2 equazione Newton (ROTAZIONE CORPO RIGIDO ATTORNO AL BARICENTRO)
    %prima costruisco i coefficienti dei momenti di rollio, beccheggio ed
    %imbardata
   
    cI_b = -0.0044919;
    cI_p = -0.55674;
    cI_r = 0.01956;

    cI =  cI_b * beta + cI_p * (p*b)/(2*V) +cI_r * (r*b)/(2*V);


    cI_1(i) = cI_b * beta;
    cI_2(i) = cI_p * (p*b)/(2*V);
    cI_3(i) = cI_r * (r*b)/(2*V);



    cN_b = 0.065829;
    cN_p = -0.023906;
    cN_r = -0.059857;
    
    cN = cN_b * beta + cN_p * (p*b)/(2*V) + cN_r * (r*b)/(2*V);


    cM_a = interp1(alpha_tab,cM_tab,alpha_deg,'linear','extrap'); 
    cM_q = -20.5;

    cM = cM_a + cM_q * (q*c_media)/(2*V);

    cM_eff(i) = cM;
    cN_eff(i) = cN;
    cI_eff(i) = cI;

    L_roll = qp*S*b*cI;
    M_pitch = qp*S*c_media*cM;
    N_yaw = qp*S*b*cN;

    p_dot = (L_roll - q*r*(Izz - Iyy) - Ixz*(r_dot + p*q))/Ixx;
    q_dot = (M_pitch - p*r*(Ixx - Izz))/Iyy;
    r_dot = (N_yaw - q*p*(Iyy - Ixx))/Izz;
    
    %tasso di variazione angoli rollio, beccheggio e imbardata

    phi_dot = p + tan(theta)*(q*sin(phi) + r*cos(phi));
    theta_dot = q*cos(phi) - r*sin(phi);
    psi_dot = (q*sin(phi) + r*cos(phi))/cos(theta);


    %aggiorno componenti
    
    u = u + u_dot*dt;
    v = v + v_dot*dt;
    w = w + w_dot*dt;

    p = p + p_dot*dt;
    q = q + q_dot*dt;    
    r = r + r_dot*dt;

    phi = phi + phi_dot*dt;
    theta = theta + theta_dot*dt;
    psi = psi + psi_dot*dt;


    V = sqrt(u^2 + v^2 + w^2);
    V_eff(i) = V;




    x_dot_ned = u*cos(theta)*cos(psi) + v*(sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi)) + w*(cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi));
    y_dot_ned = u*cos(theta)*sin(psi) + v*(sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi)) + w*(cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi));
    z_dot_ned = -u*sin(theta) + v*sin(phi)*cos(theta) + w*cos(phi)*cos(theta);

    x_ned = x_ned + x_dot_ned*dt;
    y_ned = y_ned + y_dot_ned*dt;
    z_ned = z_ned + z_dot_ned*dt;
    
    X(i) = x_ned;
    Y(i) = y_ned;
    Z(i) = z_ned;
    

    u_eff(i) = u;
    v_eff(i) = v;
    w_eff(i) = w;
  
    p_eff(i) = p;
    q_eff(i) = q;
    r_eff(i) = r;

    theta_eff(i) = theta;
    phi_eff(i) = phi;
    psi_eff(i) = psi;

    psi_dot_eff(i) = psi_dot;

    p_dot_eff(i) = p_dot;
    q_dot_eff(i) = q_dot;
    M_eff(i) = M_pitch;
    w_dot_eff(i) = w_dot;
    qp_eff(i) = qp;
    y_dot_ned_eff(i) = y_dot_ned;


    L_roll_eff(i) = L_roll;
    N_yaw_eff(i) = N_yaw;
end


figure(2)
plot(time_eff, alpha_eff)
grid on
title('grafico tempo-alpha')

figure(3)
plot(time_eff, u_eff)
grid on
title('grafico tempo-u')

figure(4)
plot(time_eff, v_eff)
grid on
title('grafico tempo-v')

figure(5)
plot(time_eff, r_eff)
grid on
title('grafico tempo-p')

figure(6)
plot(time_eff, r_eff)
grid on
title('grafico tempo-q')

figure(7)
plot(time_eff, M_eff)
grid on
title('grafico tempo- M/cM/q dot')
hold on 
plot(time_eff, cM_eff)
hold on
plot(time_eff,q_dot_eff)
legend('M','cM','q dot')        

figure(8)
plot(time_eff, beta_eff)
grid on
title('grafico tempo - beta')

figure(9)
plot(time_eff, theta_eff)
grid on
title('grafico tempo - theta')

figure(10)
plot(time_eff, phi_eff)
grid on
title('grafico tempo - phi')

figure(1)
plot3(Y,X,-Z)
xlim([-3000 3000])
ylim([-3000 3000])
zlim([-3000 3000])
grid on
title('Grafico traiettoria')