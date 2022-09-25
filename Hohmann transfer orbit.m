% Alek Rudy
%    AAE 340 HW5
%    Problem 3: Transfer orbit from LEO to GEO
%    Date: Feb. 17, 2022
%% Part c
clc;clearvars
close all

%constants
R_e = 6379;  %km
u = 3.986E5;  %km^3/s^2
R_leo = 300;

%Initial and final conditions
r_0 = R_e + R_leo;
rdot_0 = 0;
theta_0 = 0;
thetadot_0 = 0.001521;
r_f = 42241.08;
theta_circ = linspace(0,2*pi,1000); %define circle of theta = 0 to 2pi

%time in state space
t_0 = 0;
t_f = 19034.3; 
t_state = [t_0 t_f];

%Initial state space
x_state = [r_0; rdot_0; theta_0; thetadot_0];  %also called x_0

%Integrate with ode113
options = odeset('RelTol',1E-12,'AbsTol',1E-12);
[T,X] = ode113(@OrbMech_EOMS, t_state, x_state,options, u);

%Found that ode113 is better than ode45 bc its a non-stiff d.e. solver in
%variable order method. It is the best for solving computational instensive
%ODEs like orbital mechanic problems. It works the same for the Transfer
%orbit plot but not for the energy comparison plot

%Numeric solution
r = X(:,1);
rdot = X(:,2);
theta = X(:,3);
thetadot = X(:,4);

%cartesian numerical
x_num = r.*cos(theta);
y_num = r.*sin(theta);

%LEO orbit (300km) in cartesian coord.
LEO_x = r_0*cos(theta_circ);
LEO_y = r_0*sin(theta_circ);

%GEO orbit (42241km) in cartesian coord.
GEO_x = r_f*cos(theta_circ);
GEO_y = r_f*sin(theta_circ);

%Transfer orbit plot
figure(1)
plot(0,0,'.k') %earth
hold on
plot(LEO_x, LEO_y,'g','linewidth', 2) %LEO orbit
plot(GEO_x, GEO_y,'b','linewidth', 2) %GEO orbit
plot(x_num,y_num,'r','linewidth', 1.5)
axis([-5.5E4 5.5E4 -4.3E4 4.3E4])
grid on
legend('Earth', 'LEO','GEO','Transfer orbit')
xlabel('$x$ [$km$]','Interpreter','latex')
ylabel('$y$ [$km$]','Interpreter','latex')
title('Transfer orbit from LEO to GEO [Alek Rudy]')


%Energy Eq.
e_0 = (0.5)*(rdot_0^2 + (r_0^2 * thetadot_0^2)) - u/r_0;
e_num = (0.5).*(rdot.^2 + (r.^2 .* thetadot.^2)) - u./r;
energy = e_0 - e_num;


%Plot of energy computed vs actual energy
figure(2) 
plot(T,energy,'Linewidth', 1.5)
xlabel('t [sec]')
ylabel('E_{error}    ')
title('Specific Energy Computed vs Actual Specific Energy   [Alek Rudy]')
grid on


    function xdot = OrbMech_EOMS(t,x,u)
    % EOM's in state variable form

    xdot(1,1) = x(2,1);
    xdot(2,1) = x(1,1)*x(4,1)^2 - (u/x(1,1)^2);
    xdot(3,1) = x(4,1);
    xdot(4,1) = -2*x(2,1)*x(4,1)/x(1,1);
    
    end
