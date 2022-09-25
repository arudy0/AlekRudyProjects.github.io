% Alek Rudy
%    AAE 340 HW8
%    Problem 2: Plots of Binary Star System
%    Date:  Mar. 10, 2022
clc;clearvars
close all

%define constants
m1 = 1.987E30; %kg    (Star 1 mass)
m2 = m1/2; %kg        (Star 2 mass)
G = 6.673E-20; %km^3/kg*s^2    (Grav parameter) 
%% Part c
t = 0:1E6:1E9;  %to to tf (sec)

% I.C.
u1 = 0; %km         (R 1)
u2 = 0;
u3 = 0;  %km        (R 2)
u4 = 1E9;
u5 = 2; %km/s    (Rdot 1)
u6 = 0;
u7 = 12; %km/s   (Rdot 2)
u8 = 4;
IC = [u1, u2, u3, u4, u5, u6, u7, u8];

% Numerically Integrate (ode45)
options = odeset('RelTol',1E-12,'AbsTol',1E-12);
[T,State] = ode45(@BinaryStars,t,IC,options,G,m1,m2);

x1 = State(:,1);
y1 = State(:,2);
x2 = State(:,3);
y2 = State(:,4);
xdot1 = State(:,5);
ydot1 = State(:,6);
xdot2 = State(:,7);
ydot2 = State(:,8);

% Center of mass
CM_x = (2*x1 + x2)/3;
CM_y = (2*y1 + y2)/3;

CM_xdot = (2*xdot1 + xdot2)/3;
CM_ydot = (2*ydot1 + ydot2)/3;

CM_P1_x = x1 - CM_x;
CM_P1_y = y1 - CM_y;
CM_P2_x = x2 - CM_x;
CM_P2_y = y2 - CM_y;


%% Plot for part c

figure(1)
plot(x1,y1,'.k') %Star 1
hold on
grid on
plot(x2,y2,'.b') %Star 1
legend('Star 1 Trajectory','Star 2 Trajectory')
title('Binary Star system Motions [Alek Rudy]')
axis equal
xlabel('x [km]')
ylabel('y [km]')

%% Plot for part d

figure(2)
plot(CM_x,CM_y,'m', 'LineWidth', 1.5)
grid on
title('Motion of Center of Mass [Alek Rudy]')
axis equal
xlabel('x^{OC} [km]')
ylabel('y^{OC} [km]')

%% Plot for part e

figure(3)
plot(CM_xdot,CM_ydot,'*g', 'LineWidth', 1.25)
grid on
title('Velocity of the Center of Mass [Alek Rudy]')
axis([0 6 0 6])  %make axis equal
xlabel('V_x^{OC} [km/s]')
ylabel('V_y^{OC} [km/s]')

%% Plot for part f

figure(4)
plot(CM_P1_x,CM_P1_y,'b-.', 'LineWidth', 1.5)
hold on
plot(CM_P2_x,CM_P2_y,'r-.', 'LineWidth', 1.5)
grid on
title('Binary Stars Motion with Respect to Center of Mass  [Alek Rudy]')
axis equal  
xlabel('x - x^{OC} [km]')
ylabel('y - y^{OC} [km]')
legend('Star 1','Star 2')


    function udot = BinaryStars(t,u,G,m1,m2)
    % EOM's in state variable form

    udot(1) = u(5);
    udot(2) = u(6);
    udot(3) = u(7);
    udot(4) = u(8);
    udot(5) = G*m2*(u(3)-u(1))/ (((u(1)-u(3))^2 + (u(2)-u(4))^2)^(3/2));
    udot(6) = G*m2*(u(4)-u(2))/ (((u(1)-u(3))^2 + (u(2)-u(4))^2)^(3/2));
    udot(7) = G*m1*(u(1)-u(3))/ (((u(1)-u(3))^2 + (u(2)-u(4))^2)^(3/2)); 
    udot(8) = G*m1*(u(2)-u(4))/ (((u(1)-u(3))^2 + (u(2)-u(4))^2)^(3/2)); 

    udot = udot';
    return
    end
