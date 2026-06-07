close all;
clear;
clc;

%PARAMETRI
m = 300;
b = 8;
S = 10;
Ar = b^2/S;
e = 0.75;
rho0 = 1.225;
g = 9.81;
cD0 = 0.02;
k = 1/(e*pi*Ar);

%FASI VOLO
T_manovra1 = 50;
T_manovra2 = 20;
T_manovra3 = 60;
T_manovra4 = 20;
T_manovra5 = 30;
T_manovra6 = 30;
T_manovra7 = 20;
T_manovra8 = 70;
T_manovra9 = 50;
T_manovra = [T_manovra1, T_manovra2, T_manovra3, T_manovra4, T_manovra5, T_manovra6, T_manovra7, T_manovra8, T_manovra9];


t_eff = zeros(1, 9); 
t_eff(1) = T_manovra(1);

for j = 2:9
    t_eff(j) = t_eff(j-1) + T_manovra(j);
end











dt = 0.001;
tmax = 350;
time = 0:dt:tmax;


%INIZIALIZZO I VETTORI
gamma = 0;
csi = 0;
alpha = 0;
mhu = deg2rad(0);
v = 1;
vx_ned = 0;
vy_ned = 0;
vz_ned = 0;
x_ned = 0; 
y_ned = 0;
z_ned = 0;
T = 0;
z0 = 11000;
cD = 0;
cL = 0;


for i = 1:length(time)

    t = time(i);
    time_eff(i) = t;
    rho = rho0*exp(z_ned/z0);
    rho_eff(i) = rho;

    D = 0.5 * rho * v^2 * S * cD;
    L = 0.5 * rho * v^2 * S * cL;
    W = m*g;
    
    
    nzb = L/W;
    nzb_eff(i) = nzb;
    
    %calcolo v_dot, gamma_dot, csi_dot
    v_dot = (T*cos(alpha)-D-W*sin(gamma))/(m);
    csi_dot = (L*sin(mhu))/(m*v*cos(gamma));
    gamma_dot = (T*sin(alpha)+L*cos(mhu)-W*cos(gamma))/(m*v);
    


    %aggiorno angoli + velocità
    gamma = gamma + gamma_dot*dt;
    csi = csi + csi_dot*dt;

    v = v + v_dot*dt;
    

    %simulo virata corretta nel piano orizzontale 
    %(assumo per semplicità virata di R = 500m)
    
    if t < t_eff(1) 
        v = 150/3.6;
        R = 500;
        
        cL = 2*m*sqrt(v^4/R^2 + g^2)/(rho*(v^2)*S);
        cD = cD0+k*cL^2;
        mhu = atan(v^2/(g*R));
        T = 0.5*rho*(v^2)*S*cD;

    end
    
    %simulo volo rettilineo nel piano orizzontale
    if t > t_eff(1) && t < t_eff(2) || (t > t_eff(3) && t < t_eff(4)) || (t > t_eff(8))
        
        v = 150/3.6;
        
        cL = 2*m*g/(rho*(v^2)*S);
        cD = cD0 + k*cL^2;
        T = 0.5*rho*(v^2)*S*cD;
        mhu = 0;
    
    end
    

    %simulo loop nel piano verticale
    if t > t_eff(2) && t < t_eff(3)
        v = 150/3.6;
        
        R = v*(T_manovra3)/(2*pi);
        cL = 2*(m*g*cos(gamma)+m*(v^2)/R)/(rho*(v^2)*S);
        cD = cD0 + k*cL^2;
        T = 0.5*rho*(v^2)*S*cD + m*g*sin(gamma);
        if T < 0
            T = 0;
        
        %noto infatti che T < 0, il che non è possibile chiaramente.
        %allora decido di regolare cD(per esempio attivando degli aerofreni
        %in modo tale che risulti T = 0;
        %cioè non calcolo più cD dalla polare, ma da: (praticamente
        %sovrascrivo cD)
        cD = (-2*m*g*sin(gamma))/(rho*(v^2)*S);
        end
      
    
    end

    %provo a simulare manovra di salita nel piano verticale 
    %con gamma_dot tale che alla fine della manovra risulti gamma = 45.
    if t > t_eff(4) && t < t_eff(5)
        csi_dot = 0;
        v = 150/3.6;  
        gamma_dot = (deg2rad(45) - (gamma_eff(t_eff(4)/dt)-2*pi))/(T_manovra5);
        
        cL = 2*(m*v*gamma_dot + m*g*cos(gamma))/(rho*(v^2)*S);
        cD = cD0 + k*cL^2;
        T = 0.5*rho*(v^2)*S*cD + m*g*sin(gamma);
   
    
    end
    
    %simulo spirale con asse verticale (facendo un giro completo in 10s)
    if t > t_eff(5) && t < t_eff(6)
        v = 150/3.6;
        n_giri = 3;
        gamma_dot = 0;
        csi_dot = 2*pi*n_giri/T_manovra6;
        
        mhu = atan(v*csi_dot/g);
        cL = (2*m*g*cos(gamma))/(rho*(v^2)*S*cos(mhu));
        cD = cD0 + k*cL^2;
        T = 0.5*rho*(v^2)*S*cD + m*g*sin(gamma);
    end

    %riporto in moto rettilineo nel piano orizzontale in 10s (praticamente
    %è una discesa nel piano
    %verticale progettata in modo tale che gamma(10s) = 0, così poi posso
    %fare volo rettilineo orizzontale

   
    if t >= t_eff(6) && t < t_eff(7)
       
        mhu = 0;
        gamma_dot = -deg2rad(45)/T_manovra7; 
        %csi_dot = 0;
        v = 150/3.6;
        
        cL = 2*(m*g*cos(gamma)+m*v*gamma_dot)/(rho*(v^2)*S);
        cD = cD0+k*cL^2;
        T = 0.5*rho*(v^2)*S*cD + m*g*sin(gamma);

    end
    
    %ora simulo spirale con asse orizzontale, in cui csi_dot e gamma_dot
    %variano sinusoidalmente
    
    if t > t_eff(7) && t < t_eff(8)
        n_giri2 = 7;
        v = 150/3.6;
        A = deg2rad(30);
        w = 2*pi*n_giri2/T_manovra8;
        
        csi_dot = A*w*cos(w*(t-210));
        gamma_dot = -A*w*sin(w*(t-210));
        mhu = atan((v*cos(gamma)*csi_dot)/(v*gamma_dot + g*cos(gamma)));
        cL = (2*m*v*csi_dot*cos(gamma))/(rho*(v^2)*S*sin(mhu));
        
        
        
        cD = cD0 + k*cL^2;
        T = 0.5*rho*(v^2)*S*cD + m*g*sin(gamma);
        if T < 0
            T = 0;
       
        %noto infatti che T < 0, il che non è possibile chiaramente.
        %allora decido di regolare cD(per esempio attivando degli aerofreni
        %in modo tale che risulti T = 0;
        %cioè non calcolo più cD dalla polare, ma da: (praticamente
        %sovrascrivo cD)
        cD = (-2*m*g*sin(gamma))/(rho*(v^2)*S);
        end

    
    end







    gamma_dot_eff(i) = gamma_dot;
    csi_dot_eff(i) = csi_dot;
    
    csi_eff(i) = csi;
  
    
    cl_eff (i) = cL;
    cd_eff (i) = cD;
    mhu_eff (i) = mhu;
    gamma_eff (i) = gamma;
    
    
    T_eff(i) = T;
    v_eff(i) = v;
    P = T*v;
    P_eff(i) = P;

    Efficienza = abs(cL/cD);
    efficienza_min (i) = 1/(2*sqrt(cD0*k));
    eff(i) = Efficienza;
    
    if cL > 1.5
        fprintf('stallo %.2f\n', time(i));
        
        break
    end


    
    
    %velocità in ned 
    vx_ned = v*cos(gamma)*cos(csi);
    vy_ned = v*cos(gamma)*sin(csi);
    vz_ned = -v*sin(gamma);

    %posizioni in ned 
    x_ned = x_ned + vx_ned*dt;
    y_ned = y_ned + vy_ned*dt;
    z_ned = z_ned + vz_ned*dt;
    
    %salvo i vettori
    X(i) = x_ned;
    Y(i) = y_ned;
    Z(i) = z_ned;




end


figure(1)
plot3(Y,X,-Z)
xlim([-7000 3000])
ylim([-5000 3000])
zlim([-1000 3000])
grid on
title('Grafico traiettoria')

figure(2)
plot(time_eff,rho_eff)
grid on
title('Grafico tempo-densità atmosferica')
xlabel('Tempo [s]')
ylabel('Densità atmosferica [Kg/m^3]')

figure(3)
plot(time_eff,nzb_eff)
grid on
title('Grafico tempo-fattore di carico')
xlabel('Tempo [s]')
ylabel('Fattore di carico')

figure(4)
plot(time_eff,T_eff)
grid on
title('Grafico tempo-spinta')
xlabel('Tempo [s]')
ylabel('Spinta richiesta [N]')

figure(5)
plot(time_eff,P_eff)
grid on
title('Grafico tempo-potenza')
xlabel('Tempo [s]')
ylabel('Potenza richiesta [W]')

figure(6)
plot(time_eff,cl_eff)
grid on
title('Grafico tempo-cL')
xlabel('Tempo [s]')
ylabel('Coefficiente di Lift')

figure(7)
plot(time_eff,cd_eff)
grid on
title('Grafico tempo-cD')
xlabel('Tempo [s]')
ylabel('Coefficiente di Drag')

figure(8)
plot(time_eff,mhu_eff)
grid on
title('Grafico tempo-angolo di rollio')
xlabel('Tempo [s]')
ylabel('Angolo di rollio [rad]')

figure(9)
plot(time_eff,eff)
grid on
title('Grafico tempo-efficienza')
xlabel('Tempo [s]')
ylabel('Efficienza')
hold on 
plot(time_eff,efficienza_min)
legend ('Efficienza reale','Efficienza massima')