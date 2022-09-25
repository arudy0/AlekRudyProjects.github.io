%    PSP-Liquids: Traj
%    Finalized Radial Heat Transfer

clc
clear all
close all

innerRadius = .14986;
outerRadius = .1524;
thermal1 = .175;
thermal2 = 20;
convectiveHeatK = 180.1; %52
tankLength = 1.218;
tInf1 = 71;
tInf2 = 298;
Cp = 1.676;
mass = 140;
t = 0;
dt = 0.001;

timeVals = [t];
tInf1Vals = [tInf1];
deltaTVals = [0];

while tInf1 <= 295

    rConv1 = 1 / (2 * pi * innerRadius * tankLength * thermal1);
    rConv2 = 1 / (2 * pi * outerRadius * tankLength * thermal2);
    rCyl = (log(outerRadius / innerRadius)) / (2 * pi * tankLength * convectiveHeatK);
    rTotal = rConv1 + rConv2 + rCyl;
    qDot = (tInf2 - tInf1) / rTotal;

    deltaT = qDot / (mass * Cp) * dt;
    tInf1 = tInf1 + deltaT;
    t = t + dt;

    timeVals(end + 1) = t;
    tInf1Vals(end + 1) = tInf1;
    deltaTVals(end + 1) = deltaT;
        
end

figure
plot(timeVals, tInf1Vals)
title('LOX temperature vs. Time')
xlabel('Time [s]')
ylabel('Temperature [K]')
axes('pos',[.1 .75 .3 .1])
imshow('PSP2.png')

figure
plot(timeVals, deltaTVals)
title('Temperature Change vs. Time')
xlabel('Time [s]')
ylabel('Temperature Change [K]')
axes('pos',[.1 .75 .3 .1])
imshow('PSP2.png')
