%% Prepare workspace
clear all %#ok<CLALL>
close all
clc

%% Input parameters
vehicle(1).m    = 1031.9; % kg
vehicle(1).Iz   = 1850; % kg-m^2
vehicle(1).a    = 0.9271; % Distance from front axle to CG, in meters
vehicle(1).b    = 1.5621; % Distance from rear axle to CG, in meters
vehicle(1).Caf  = -77500; % N/rad;
vehicle(1).Car  = -116250; % N/rad;
U = 20; % U is forward velocity of vehicle in longitudinal direction, [m/s] (rule of thumb: mph ~= 2* m/s)

% Items used to define steering inputs
steering_amplitude_degrees = 2; % 2 degrees of steering amplitude for input sinewave
Period = 3; % Units are seconds. A typical lane change is about 3 to 4 seconds based on experimental highway measurements

% Define items used to determine how long to run sim, number of time
% points, etc.
total_duration = 1.5*Period; % This is how long the simulation will run. Usually 1.5 times the period is enough.
delta_t        = 0.01; % This is the time step of the simulation. 
number_of_simulation_steps = floor(total_duration/delta_t)+1; % This is the number of time steps we should have

% Initialize variables
V = nan(1,number_of_simulation_steps);
r = nan(1,number_of_simulation_steps);
X = nan(1,number_of_simulation_steps);
Y = nan(1,number_of_simulation_steps);
Phi = nan(1,number_of_simulation_steps);
V(1) = 0;
r(1) = 0;
X(1) = 0;
Y(1) = 0;
Phi(1) = 0;
counter = 1;
for time = delta_t:delta_t:total_duration
    delta_f = (1-1*(0<time-Period))*...
              (pi/180)*steering_amplitude_degrees*sin((2*pi/Period)*time);
    alpha = fcn_slipAngles(U, V(counter), r(counter), delta_f, vehicle);
    Fy = fcn_lateralForces(alpha, vehicle);
    [~, y]  = ode45(@(t,y) fcn_lateralDynamics(t, y, U, Fy, vehicle), ...
                    [0 delta_t], [V(counter), r(counter)]);
    counter    = counter+1;
    V(counter) = y(end,1);
    r(counter) = y(end,2);
    
    [~, y]  = ode45(@(t,y) fcn_body2globalCoordinates(t, y, U, V(counter), r(counter)), ...
                    [0 delta_t], [X(counter-1), Y(counter-1), Phi(counter-1)]);
    X(counter) = y(end,1);
    Y(counter) = y(end,2);
    Phi(counter) = y(end,3);
end

figure(12345)
clf
plot(X,Y)