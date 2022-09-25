% Alek Rudy
%    AAE 340 HW6
%    Problem 1: Part 3: Plots of Galileo Orbiter & probe around Jupiter
%    Date: Feb. 22, 2022
clc;clearvars
close all

%define constants
Rj = 71398; %km    (jupiter radius)
V_inf = 5.455; %km/s   (orbiter and probe)
mu = 1.267E8; %km^3/s^2    (Jupiter's grav parameter) 

%% Jupiter dimensions
circ = linspace(0,2*pi,1000);
Xj = Rj*cos(circ);
Yj = Rj*sin(circ);

%% Orbiter Elliptic Trajectory
t = [0, 198*24*3600];  %to to tf (sec)

% I.C.
r_0 = 285592; %km
rdot_0 = 0;  %km/s
theta_0 = 0; %rad
thetadot_0 = 1.03536E-4; %rad/s
IC = [r_0, rdot_0, theta_0, thetadot_0];

% Numerically Integrate (ode45)
options = odeset('RelTol',1E-12,'AbsTol',1E-12);
[T1,State1] = ode45(@GalileoTraj,t,IC,options,mu);
r_1 = State1(:,1);
rdot_1 = State1(:,2);
theta_1 = State1(:,3);
thetadot_1 = State1(:,4);

%Polar --> cartesian (for plotting)
x_1 = r_1.*cos(theta_1);
y_1 = r_1.*sin(theta_1);

%% Orbiter Hyperbolic Trajectory
t = [0, 32*24*3600];  %to to tf (sec for 32 days)

% I.C.
r_0 = 300*Rj; %km
rdot_0 = -6.4362;  %km/s
theta_0 = -2.7173; %rad
thetadot_0 = 1.8851E-8; %rad/s
IC = [r_0, rdot_0, theta_0, thetadot_0];

% Numerically Integrate (ode45)
options = odeset('RelTol',1E-12,'AbsTol',1E-12);
[T2,State2] = ode45(@GalileoTraj,t,IC,options,mu);
r_2 = State2(:,1);
rdot_2 = State2(:,2);
theta_2 = State2(:,3);
thetadot_2 = State2(:,4);

%Polar --> cartesian (for plotting)
x_2 = r_2.*cos(theta_2);
y_2 = r_2.*sin(theta_2);

%% Probe Hyperbolic Trajectory
t = [0, 32*24*3600];  %to to tf (sec for 32 days)

% I.C.
r_0 = 300*Rj; %km
rdot_0 = -6.4557;  %km/s
theta_0 = -2.9262; %rad
thetadot_0 = 9.3099E-9; %rad/s
IC = [r_0, rdot_0, theta_0, thetadot_0];

% Numerically Integrate (ode45)
options = odeset('RelTol',1E-12,'AbsTol',1E-12);
[T3,State3] = ode45(@GalileoTraj,t,IC,options,mu);
r_3 = State3(:,1);
rdot_3 = State3(:,2);
theta_3 = State3(:,3);
thetadot_3 = State3(:,4);

%Polar --> cartesian (for plotting)
x_3 = r_3.*cos(theta_3);
y_3 = r_3.*sin(theta_3);


%% Plot
figure(1)
plot(Xj,Yj,'k')
hold on
grid on
plot(x_1,y_1,'r-')  %Ellipse 
plot(x_2,y_2,'b:','Linewidth', 1.5) %Orbiter Hyperbolic
plot(x_3,y_3,'g-.','Linewidth', 1.5) %Probe Hyperbolic
legend('Jupiter','Orbiter Ellipse','Orbiter Hyperbolic Trajectory','Probe Hyperbolic Trajectory')
title('Galileo Satellite and Probe entering Jupiters Orbit [Alek Rudy]')
axis equal
xlabel('x [km]')
ylabel('y [km]')



    function xdot = GalileoTraj(t,x,mu)
    % EOM's in state variable form


    xdot(1) = x(2);
    xdot(2) = x(1)*x(4)^2 - (mu/x(1)^2);
    xdot(3) = x(4);
    xdot(4) = -2*x(2)*x(4)/x(1);

    xdot = xdot';
    return
    end



