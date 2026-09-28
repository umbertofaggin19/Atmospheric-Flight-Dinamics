%6DOF in sistema di riferimento BODY

%SI NOTI CHE IL MODO DI SPIRALE è INSTABILE (LENTA, ~ 200/300 sec(autovalori circa lambda = 0.002,quindi normale che sia lento) , PERTANTO VA CORRETTA NEL TEMPO
%CON INPUT DEL PILOTA, tramite i vari delta----> implementare PID per poter
%regolare nel tempo alettoni/elevatori/timone
%by a rolling r heading disturbance, i.e. p or r =! 0)


%instead fugoid, short mode and dutch mode seems stable


%roll damping should be stable but seems unstable. in realtà roll damping
%funziona bene, è solo che eccitando p suscito sia roll damping che spiral
%mode. Roll damping funziona bene (p smorzato nei primi 5s), poi subentra
%spiral mode (si vede nel lungo periodo), quindi tutto ok

clc;
close all;
clear;



%PARAMETRI
m = 0.75;
g = 9.81;
b = 1.4;
S = 0.205;

%(dati relativi all'inerzia, supp Ixz = 0)
Ixx = 0.04814;
Iyy = 0.03591;
Izz = 0.08371;
Ixz = 0.0011;
c_media = 0.151;

rho0 = 1.225;
z0 = 11000;
Z0 = -z0;

%STATO INIZIALE
alpha = deg2rad(0);
beta = deg2rad(0);

phi = 0;     %rollio
gamma = deg2rad(0);
theta = alpha + gamma;   %beccheggio
psi = deg2rad(0);     %imbardata

p = deg2rad(0);    %roll rate
q = deg2rad(0);    %pitch rate
r = deg2rad(0);    % yaw rate


V = 50;
u = V*cos(alpha)*cos(beta);
v = V*sin(beta);
w = V*sin(alpha)*cos(beta);


x_ned = 0;
y_ned = 0;
z_ned = 0;
Z = 0;



dt = 0.001;
tmax = 1000;
time = 0:dt:tmax;





errore_cumulato1 = 0;
V_old = 19;

errore_cumulato2 = 0;
errore_precedente2 = 0;

errore_cumulato3 = 0;
errore_precedente3 = 0;

errore_cumulato4 = 0;
errore_precedente4 = 0;

errore_cumulato5 = 0;
errore_precedente5 = 0;

errore_cumulato6 = 0;
errore_precedente6 = 0;

for i = 1:length(time)
    
    t = time(i);
    time_eff(i) = t;
    
    %ATMOSFERA STANDARD
    
    rho0 = 1.225;
    h0 = 11000;
    h = z_ned;
    rho = rho0*exp(h/h0);
    rho_eff(i) = rho;

    qp = 0.5*rho*(V^2);
    

    alpha = atan2(w,u);
    alpha_eff(i) = alpha;
    alpha_trim = 0.0131;
    alpha_trim_eff(i) = alpha_trim;

    

    %5° PID per regolare quota mediante theta mediante elevatori
    
    if t < 36.5
        h_setpoint = 0 ; %segno + = verso l'alto
    elseif t > 36.5 && t < 100
        h_setpoint = 100;
    elseif t > 100 && t < 300
        h_setpoint = 600;
    elseif t > 300
        h_setpoint = 200;
    end

    errore5 = h_setpoint + z_ned;
    
    kP5 = 0.001;
    kI5 = 0;
    kD5 = 0;

    P5 = kP5 * errore5;

    errore_cumulato5 = errore_cumulato5 + errore5*dt;
    I5 = kI5 * errore_cumulato5;
    
    D5 = kD5*(errore5 - errore_precedente5)/dt;


    theta_comando = P5 + I5 + D5;

    errore2 = theta_comando - theta;

    kP2 = 10;
    kI2 = 2;
    kD2 = 0;

    P2 = kP2 * errore2;
    
    errore_cumulato2 = errore_cumulato2 + errore2*dt;
    I2 = kI2 * errore_cumulato2;

    D2 = kD2*(errore2 - errore_precedente2)/dt;
    errore_precedente2 = errore2;

    delta_e = -(P2 + I2 - D2);
    

    beta = atan2(v,sqrt(u^2 + w^2));  
    beta_eff(i) = beta;

    % 1° PID per regolare V mediante Tx.
    if t < 36.5
        V_setpoint = 40;
    elseif t > 36.5 && t < 100
        V_setpoint = 50;
    elseif t > 100 && t < 137.5
        V_setpoint = 40;
    elseif t > 137.5
        V_setpoint = 40;
    end
    
    
    
    
    V_setpoint_eff(i) = V_setpoint;
    errore1 = V_setpoint - V;

    kP1 = 5;
    kI1 = 3;
    kD1 = 0.1;

    P1 = kP1 * errore1;

    errore_cumulato1 = errore_cumulato1 + errore1*dt;
    I1 = kI1 * errore_cumulato1;

    D1 = kD1*(V-V_old)/dt;
    V_old = V;

    Tx = P1 + I1 - D1;

    Ty = 0;
    Tz = 0; 
    T = sqrt((Tx)^2 + (Ty)^2 + (Tz)^2);
    Tx_eff(i) = Tx;
   
    % 3° PID per regolare Y mediante PSI mediante PHI mediante DELTA_A (alettoni)
    
    if t<100
        y_setpoint = 200;
    elseif t>100 && t<400
        y_setpoint = 600;
    else 
        y_setpoint = 0;
    end

    errore6 = y_setpoint - y_ned; %RICONTROLLARE QUESTA RIGA
    errore6_eff(i) = errore6;

    kP6 = 0.001;
    kI6 = 0;
    kD6 = 0;

    P6 = kP6 * errore6;

    errore_cumulato6 = errore_cumulato6 + errore6*dt;
    I6 = kI6 * errore_cumulato6;

    D6 = kD6 * (errore6 - errore_precedente6)/dt;
    errore_precedente6 = errore6;

    psi_comando = P6 + I6 + D6; %PHI NECESSARIO PER VIRARE PER RAGGIUNGERE QUELLA Y
    
    errore7 = psi_comando - psi;
    kP7 = 2;
    phi_comando = kP7 * errore7;

    errore3 = phi_comando - phi;

    kP3 = 1;
    kI3 = 0;
    kD3 = 0;

    P3 = kP3*errore3;

    errore_cumulato3 = errore_cumulato3 + errore3*dt;
    I3 = kI3*errore_cumulato3;

    D3 = kD3*(errore3 - errore_precedente3)/dt;
    errore_precedente3 = errore3;

    delta_a = P3 + I3 - D3;


    % 4° PID per regolare psi(o beta) mediante timone (rudder) (delta_r)
    
    beta_setpoint = deg2rad(0);
    beta_setpoint_eff(i) = beta_setpoint;
    errore4 = beta_setpoint - beta;

    kP4 = 5;
    kI4 = 1;
    kD4 = 0;

    P4 = kP4*errore4;

    errore_cumulato4 = errore_cumulato4 + errore4*dt;
    I4 = kI4*errore_cumulato4;

    D4 = kD4*(errore4 - errore_precedente4)/dt;
    errore_precedente4 = errore4;

    delta_r = P4 + I4 - D4;
    
    %calcolo coefficienti aerodinamici cX, cY, cZ;

    cX0 = -0.009180;
    cXa = 0.052058;
    cXdelta_e = -0.0033842;
    cXdelta_a = -3.1 * 10^-4;
    cXdelta_r = -3.0426 * 10^-5;
    cX = cX0 + cXa * alpha  + cXdelta_e * delta_e + cXdelta_a * delta_a + cXdelta_r * delta_r;

    cY_b = -0.14908;
    cY_p = 0.02272;
    cY_r = 0.1395;
    cYdelta_e = 1.6143 * 10^-5;
    cYdelta_a = -2.94 * 10^-2;
    cYdelta_r = 0.08468;
    cY = cY_b*beta + cY_p * (p*b)/(2*V) + cY_r * (r*b)/(2*V) + cYdelta_e * delta_e + cYdelta_a * delta_a + cYdelta_r * delta_r;
    
    cZ0 = -0.011228;
    cZa = -5.507; 
    cZq = -10.21;
    cZdelta_e = -0.46773;
    cZdelta_a = 3.32 * 10^-5;
    cZdelta_r = 4.6428 * 10^-6;
    cZ = cZ0 + cZa * alpha  + cZq * (q*c_media)/(2*V) + cZdelta_e * delta_e + cZdelta_a * delta_a + cZdelta_r * delta_r;

    
    Xa = qp*S*cX;
    Ya = qp*S*cY;
    Za = qp*S*cZ;
    

    %1 equazione Newton (TRASLAZIONE BARICENTRO)
   
    %simulo vento che
    %spinge indietro, di conseguenza è richiesta più spinta (5 N in più, in
    %effetti...)
    Vx = -10;
    Vy = 0;
    Vz = 3;

    fx = -m*g*sin(theta) + Tx + Xa + Vx;
    fy = m*g*sin(phi)*cos(theta) + Ty + Ya + Vy;
    fz = m*g*cos(phi)*cos(theta) + Tz + Za + Vz;


    u_dot = (r*v - q*w) + fx/m;
    v_dot = (p*w - u*r) + fy/m;
    w_dot = (q*u - p*v) + fz/m;


    %2 equazione Newton (ROTAZIONE CORPO RIGIDO ATTORNO AL BARICENTRO)
    %prima costruisco i coefficienti dei momenti di rollio, beccheggio ed imbardata
   
    cI_b = -0.0044919;
    cI_p = -0.55674;
    cI_r = 0.01956;
    cIdelta_e = 0;
    cIdelta_a = 0.56655;
    cIdelta_r = 0.0029683 ;
    
    cI =  cI_b * beta + cI_p * (p*b)/(2*V) +cI_r * (r*b)/(2*V) + cIdelta_e*delta_e + cIdelta_a*delta_a + cIdelta_r*delta_r;

    cM_0 = 0.021068;
    cM_a = -1.6138;
    cM_u = -0.00048906;
    cM_q = -20.5;
    cMdelta_e = -1.7221;
    cMdelta_a = 5.29 * 10^-5;
    cMdelta_r = -2.7444 * 10^-5;

    cM = cM_0 + cM_a*alpha  + cM_q * (q*c_media)/(2*V) + cMdelta_e*delta_e + cMdelta_a*delta_a + cMdelta_r*delta_r;


    cN_b = 0.065829;
    cN_p = -0.023906;
    cN_r = -0.059857;
    cNdelta_e = 0;
    cNdelta_a = 0.014587;
    cNdelta_r = -0.038961;

    cN = cN_b * beta + cN_p * (p*b)/(2*V) + cN_r * (r*b)/(2*V) + cNdelta_e*delta_e + cNdelta_a*delta_a + cNdelta_r*delta_r;

    cM_eff(i) = cM;
    cN_eff(i) = cN;
    cI_eff(i) = cI;

    L_roll = qp*S*b*cI;
    M_pitch = qp*S*c_media*cM;
    N_yaw = qp*S*b*cN;

    p_dot = (L_roll - q*r*(Izz - Iyy))/Ixx;
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

figure(1)
plot3(Y,X,-Z)
xlim([-3000 3000])
ylim([-3000 3000])
zlim([-3000 3000])
grid on
title('Grafico traiettoria')

figure(2)
hold on
grid on
plot(time_eff,errore6_eff)
xlabel('tempo')
ylabel('errore 6')
title('tempo-errore 6')



