% Alek Rudy
%    AAE 439 HW3 Prob 3.5
%    2-D Ballistic Missile Trajectory
%    Date: 9/18/2022
clc;clearvars
close all
%% Initialize

g = 32.2; %f/s^2
Cd = 0.3;
Cl = 0;
Dref = 3.1/12; %ft
t = linspace(0,2,41); %sec
F_t = g*(18 - 9*t); %lb
mo = 1.1; %lb
Isp = 165; %s
A_ref = Dref^2 * pi /4; %ft^2
p = 0.0765; %lb/ft^3 

%% Part a
g = 32.2; %f/s^2
mo = 1.1; %lb
Isp = 165; %s
t = linspace(0,2,41); %sec
F_t = g*(18 - 9*t); %lb


m_t = flip(mo - (F_t/(g*Isp))) ;  


F_t(41:400) = F_t(41) ;
%from element 1->41 force decreases, 41 onwards the equation F = 0
m_t(41:400) = m_t(41);
%from element 1->41 mass decreases, 41 onwards is final mass 
%Propellant burnt out Mp = 0 at t=2

%% Part b
t = linspace(0,2,400); %sec
clear t
% Initial Cond
Lrod = 5; %Length of rod (ft)
Vw = 14.67;  %velocity of wind ft/s

%% For theta = 70 deg
    theta(1) = deg2rad(70);
    V(1) = 0;
    t(1) = 0;
    dt = 0.05; %Step input
    x(1) = 0;
    z(1) = 0;

for i=1:1:399
    F(i) = F_t(i);
    m(i) = m_t(i);
    Vrel(i) = sqrt(  (V(i)*sin(theta(i)))^2 +   (Vw + V(i)*sin(theta(i)))^2    );
    D(i) = 0.5*Cd*(A_ref)*p*Vrel(i)^2;
    phi(i) = atan(  V(i)*sin(theta(i))/ (V(i)*cos(theta(i)) + Vw)  );
    phi(1) = 0;
    G(i) = ((F(i) - D(i))./m(i)) .* cos(phi(i)-theta(i)) - g.*sin(theta(i));
    r(i) = sqrt(x(i)^2 + z(i)^2);
    if r(i) < Lrod
        H(i) = 0;
    else
        H(i) = (F(i)-D(i))./(m(i).*V(i)) .* sin(phi(i)-theta(i)) - g./V(i).*cos(theta(i));
    end
    
    %Time to get *values (labeled _n for now)
    t(i+1) = t(i) + dt;
    F_n(i) = F_t(i+1);
    m_n(i) = m_t(i+1);
    V_n(i) = V(i) + dt*G(i);
    theta_n(i) = theta(i) + dt*H(i);
   
    Vrel_n(i) = Vrel(i) + dt*sqrt(  (V(i)*sin(theta(i)))^2 +   (Vw + V(i)*sin(theta(i)))^2    );
    D_n(i) = 0.5*Cd*(A_ref)*p*Vrel_n(i)^2;
    phi_n(i) = atan(  V_n(i)*sin(theta_n(i))/ (V_n(i)*cos(theta_n(i)) + Vw)  );
    G_n(i) = ((F_n(i) - D_n(i))./m_n(i)) .* cos(phi_n(i)-theta_n(i)) - g.*sin(theta_n(i));
    if r(i) < Lrod
        H_n(i) = 0;
    else
        H_n(i) = (F_n(i)-D_n(i))./(m_n(i).*V_n(i)) .* sin(phi_n(i)-theta_n(i)) - g./V_n(i).*cos(theta_n(i));
    end
        
    % Update V, theta, z, x
    V(i+1) = V(i) + (dt/2)*(G(i) + G_n(i));
    theta(i+1) = theta(i) + (dt/2)*(H(i) + H_n(i));
    z(i+1) = z(i) + (dt/2)*(V(i)*sin(theta(i)) + V(i+1)*sin(theta(i+1)));
    x(i+1) = x(i) + (dt/2)*(V(i)*cos(theta(i)) + V(i+1)*cos(theta(i+1)));
end


%% For theta = 75 deg
    theta_75 =  deg2rad(75);
    V_75(1) = 0;
    t(1) = 0;
    dt = 0.05; %Step input
    x_75(1) = 0;
    z_75(1) = 0;

for i=1:1:399
    F_75(i) = F_t(i);
    m_75(i) = m_t(i);
    Vrel_75(i) = sqrt(  (V_75(i)*sin(theta_75(i)))^2 +   (Vw + V_75(i)*sin(theta_75(i)))^2    );
    D_75(i) = 0.5*Cd*(A_ref)*p*Vrel_75(i)^2;
    phi_75(i) = atan(  V_75(i)*sin(theta_75(i))/ (V(i)*cos(theta_75(i)) + Vw)  );
    phi_75(1) = 0;
    G_75(i) = ((F_75(i) - D_75(i))./m_75(i)) .* cos(phi_75(i)-theta_75(i)) - g.*sin(theta_75(i));
    r_75(i) = sqrt(x_75(i)^2 + z_75(i)^2);
    if r_75(i) < Lrod
        H_75(i) = 0;
    else
        H_75(i) = (F_75(i)-D_75(i))./(m_75(i).*V_75(i)) .* sin(phi_75(i)-theta_75(i)) - g./V_75(i).*cos(theta_75(i));
    end
    
    %Time to get *values (labeled _n for now)
    t(i+1) = t(i) + dt;
    F_n_75(i) = F_t(i+1);
    m_n_75(i) = m_t(i+1);
    V_n_75(i) = V_75(i) + dt*G_75(i);
    theta_n_75(i) = theta_75(i) + dt*H_75(i);
    
    Vrel_n_75(i) = Vrel_75(i) + dt*sqrt(  (V_75(i)*sin(theta_75(i)))^2 +   (Vw + V_75(i)*sin(theta_75(i)))^2    );
    D_n_75(i) = 0.5*Cd*(A_ref)*p*Vrel_n_75(i)^2;
    phi_n_75(i) = atan(  V_n_75(i)*sin(theta_n_75(i))/ (V_n_75(i)*cos(theta_n_75(i)) + Vw)  );
    G_n_75(i) = ((F_n_75(i) - D_n_75(i))./m_n_75(i)) .* cos(phi_n_75(i)-theta_n_75(i)) - g.*sin(theta_n_75(i));
    if r(i) < Lrod
        H_n_75(i) = 0;
    else
        H_n_75(i) = (F_n_75(i)-D_n_75(i))./(m_n_75(i).*V_n_75(i)) .* sin(phi_n_75(i)-theta_n_75(i)) - g./V_n_75(i).*cos(theta_n_75(i));
    end
        
    % Update V, theta, z, x
    V_75(i+1) = V_75(i) + (dt/2)*(G_75(i) + G_n_75(i));
    theta_75(i+1) = theta_75(i) + (dt/2)*(H_75(i) + H_n_75(i));
    z_75(i+1) = z_75(i) + (dt/2)*(V_75(i)*sin(theta_75(i)) + V_75(i+1)*sin(theta_75(i+1)));
    x_75(i+1) = x_75(i) + (dt/2)*(V_75(i)*cos(theta_75(i)) + V_75(i+1)*cos(theta_75(i+1)));
end

%% For theta = 80 deg
    theta_80 =  deg2rad(80);
    V_80(1) = 0;
    t(1) = 0;
    dt = 0.05; %Step input
    x_80(1) = 0;
    z_80(1) = 0;

for i=1:1:399
    F_80(i) = F_t(i);
    m_80(i) = m_t(i);
    Vrel_80(i) = sqrt(  (V_80(i)*sin(theta_80(i)))^2 +   (Vw + V_80(i)*sin(theta_80(i)))^2    );
    D_80(i) = 0.5*Cd*(A_ref)*p*Vrel_80(i)^2;
    phi_80(i) = atan(  V_80(i)*sin(theta_80(i))/ (V_80(i)*cos(theta_80(i)) + Vw)  );
    phi_80(1) = 0;
    G_80(i) = ((F_80(i) - D_80(i))./m_80(i)) .* cos(phi_80(i)-theta_80(i)) - g.*sin(theta_80(i));
    r_80(i) = sqrt(x_80(i)^2 + z_80(i)^2);
    if r_80(i) < Lrod
        H_80(i) = 0;
    else
        H_80(i) = (F_80(i)-D_80(i))./(m_80(i).*V_80(i)) .* sin(phi_80(i)-theta_80(i)) - g./V_80(i).*cos(theta_80(i));
    end
    
    %Time to get *values (labeled _n for now)
    t(i+1) = t(i) + dt;
    F_n_80(i) = F_t(i+1);
    m_n_80(i) = m_t(i+1);
    V_n_80(i) = V_80(i) + dt*G_80(i);
    theta_n_80(i) = theta_80(i) + dt*H_80(i);
   
    Vrel_n_80(i) = Vrel_80(i) + dt*sqrt(  (V_80(i)*sin(theta_80(i)))^2 +   (Vw + V_80(i)*sin(theta_80(i)))^2    );
    D_n_80(i) = 0.5*Cd*(A_ref)*p*Vrel_n_80(i)^2;
    phi_n_80(i) = atan(  V_n_80(i)*sin(theta_n_80(i))/ (V_n_80(i)*cos(theta_n_80(i)) + Vw)  );
    G_n_80(i) = ((F_n_80(i) - D_n_80(i))./m_n_80(i)) .* cos(phi_n_80(i)-theta_n_80(i)) - g.*sin(theta_n_80(i));
    if r_80(i) < Lrod
        H_n_80(i) = 0;
    else
        H_n_80(i) = (F_n_80(i)-D_n_80(i))./(m_n_80(i).*V_n_80(i)) .* sin(phi_n_80(i)-theta_n_80(i)) - g./V_n_80(i).*cos(theta_n_80(i));
    end
        
    % Update V, theta, z, x
    V_80(i+1) = V_80(i) + (dt/2)*(G_80(i) + G_n_80(i));
    theta_80(i+1) = theta_80(i) + (dt/2)*(H_80(i) + H_n_80(i));
    z_80(i+1) = z_80(i) + (dt/2)*(V_80(i)*sin(theta_80(i)) + V_80(i+1)*sin(theta_80(i+1)));
    x_80(i+1) = x_80(i) + (dt/2)*(V_80(i)*cos(theta_80(i)) + V_80(i+1)*cos(theta_80(i+1)));
end


%% Table & Plots

% Table below
time_tab = [t]';
X_tab = [x_75]';
Z_tab = [z_75]';
G_75(length(X_tab)) = 0;
accel_tab = [G_75]';
Vel_tab = [V_75]';
theta_tab = [rad2deg(theta_75)]';
Tables = table(time_tab, X_tab, Z_tab, accel_tab, Vel_tab, theta_tab)



figure(1)
plot(x,z,'r')
grid on
hold on
plot(x_75,z_75,'g--')
plot(x_80,z_80,'b-.')
xlabel('$x$ (ft)','Interpreter','latex');
ylabel('$z$  (ft)','Interpreter','latex');
title('$X$(t) vs $Z$(t)  [Alek Rudy]','Interpreter','latex');
set(gca,'FontSize',16);
legend('70 deg','75 deg','80 deg')
axis([0 1600 0 1600])


figure(2)
plot(t,z,'r')
grid on
hold on
plot(t,z_75,'g--')
plot(t,z_80,'b-.')
xlabel('$t$ (s)','Interpreter','latex');
ylabel('$z$  (ft)','Interpreter','latex');
title('$z$(t)  [Alek Rudy]','Interpreter','latex');
set(gca,'FontSize',16);
legend('70 deg','75 deg','80 deg')
axis([0 20 0 1200])



figure(3)
plot(t,V,'r')
grid on
hold on
plot(t,V_75,'g--')
plot(t,V_80,'b-.')
xlabel('$t$ (s)','Interpreter','latex');
ylabel('$V$  (ft/s)','Interpreter','latex');
title('$t$(t) vs $V$(t)  [Alek Rudy]','Interpreter','latex');
set(gca,'FontSize',16);
legend('70 deg','75 deg','80 deg')
axis([0 20 0 400])



