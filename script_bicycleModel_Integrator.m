% script_bicycleModel_Integrator.m
% Plots the results of a outputs of a simulink model: the all-integrator 
% form and ODE45 implemented as MATLAB functions.
% Note that both are same (which is the point). It also tests this across 
% 3 vehicles.
%
% This script was originally written on 2007_06_25 by S. Brennan
% Modified by Satya Prasad on 2021/04/29
% Questions or comments? szm888@psu.edu

%% Prep the workspace
close all;  % close all the plots
clear all %#ok<CLALL>
clc

%% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _       
%  |_   _|                 | |      
%    | |  _ __  _ __  _   _| |_ ___ 
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |                  
%              |_| 
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set core simulation inputs: the speed, and steering amplitude
U = 20;  % U is forward velocity of vehicle in longitudinal direction, [m/s] (rule of thumb: mph ~= 2* m/s)

% Items used to define steering inputs
steering_amplitude_degrees = 2; % 2 degrees of steering amplitude for input sinewave
Period = 3; % Units are seconds. A typical lane change is about 3 to 4 seconds based on experimental highway measurements

% Define items used to determine how long to run sim, number of time
% points, etc.
TotalTime = 1.5*Period;  % This is how long the simulation will run. Usually 1.5 times the period is enough.
deltaT = 0.01; % This is the time step of the simulation. See the "Model Settings" submenu in Simulink to see where this variable is used.
N_timeSteps = floor(TotalTime/deltaT) + 1; % This is the number of time steps we should have



%% Fill in the vehicle parameters. 
% Use a structure array so we can have several vehicles

% Approximately a Ford Taurus
vehicle(1).m          = 1031.9; % kg
vehicle(1).Iz         = 1850; % kg-m^2
vehicle(1).a          = 0.9271; % Distance from front axle to CG, in meters
vehicle(1).b          = 1.5621;  % Distance from rear axle to CG, in meters
vehicle(1).Caf        = -77500; % N/rad;
vehicle(1).Car	      = -116250; % N/rad;

% A random vehicle (race car tires)
vehicle(2).m          = 1670; % kg
vehicle(2).Iz         = 2100; % kg-m^2
vehicle(2).a          = 0.99; % meters
vehicle(2).b          = 1.7;  % meters
vehicle(2).Caf        = -123200; % N/rad;
vehicle(2).Car	      = -104200; % N/rad;

% A variation on the random vehicle (go-cart like)
vehicle(3).m          = 167; % kg
vehicle(3).Iz         = 210; % kg-m^2
vehicle(3).a          = 0.99; % meters
vehicle(3).b          = 1.7;  % meters
vehicle(3).Caf        = -70200; % N/rad;
vehicle(3).Car	      = -80200; % N/rad;

% The number of vehicles is the length of the structure array
N_vehicles = length(vehicle);

%% Main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This runs the codes to predict core behavior within SIMULINK and MATLAB,
% looping through each vehicle.
% For this simulation set, the core behavior is
% V: lateral speed
% r: rotational rate of vehicle around z-axis

% Initialize all the arrays with Not-a-Number (nan). When plotting, any
% values that are nan will be left empty.
all_X_sim   = nan(N_timeSteps,N_vehicles);
all_Y_sim   = nan(N_timeSteps,N_vehicles);
all_phi_sim = nan(N_timeSteps,N_vehicles);
all_r_sim   = nan(N_timeSteps,N_vehicles);
all_V_sim   = nan(N_timeSteps,N_vehicles);

all_X_mat   = nan(N_timeSteps,N_vehicles);
all_Y_mat   = nan(N_timeSteps,N_vehicles);
all_phi_mat = nan(N_timeSteps,N_vehicles);
all_r_mat   = nan(N_timeSteps,N_vehicles);
all_V_mat   = nan(N_timeSteps,N_vehicles);

% Loop through all the vehicles, simulating the trajectory of each within
% the for loop
for vehicle_i=1:N_vehicles
    % Print to the console which vehicle we are working on
    fprintf(1,'Working on vehicle: %d\n', vehicle_i);
    
    %% Using Simulink model
    % Fill in the parameters needed by the simulation for it to run
    m = vehicle(vehicle_i).m;
    Iz = vehicle(vehicle_i).Iz;
    a = vehicle(vehicle_i).a;
    b = vehicle(vehicle_i).b;
    L = a+b;
    Caf = vehicle(vehicle_i).Caf;
    Car = vehicle(vehicle_i).Car;
    
    % Run the simulation in SIMULINK
    sim('mdl_bicycleModel_Integrator.slx', TotalTime);
    
    % Save the results in a big array (for plotting in next part) 
    % Before saving, we need to check if the full vector is shorter than
    % expected length of N_timeSteps
    if length(t) ~= N_timeSteps
        warning('More time was spent than expected in the simulation. Keeping only the expected time portion.')        
    end    
    
    % Keep the shorter of either the actual length, or expected length:
    shorter_index = min(N_timeSteps,length(t));
    
    % Fill in the data arrays
    all_X_sim(1:shorter_index,vehicle_i)   = X_int(1:shorter_index);
    all_Y_sim(1:shorter_index,vehicle_i)   = Y_int(1:shorter_index);
    all_phi_sim(1:shorter_index,vehicle_i) = phi_int(1:shorter_index);
    all_r_sim(1:shorter_index,vehicle_i)   = r_int(1:shorter_index);
    all_V_sim(1:shorter_index,vehicle_i)   = V_int(1:shorter_index);
    
    %% ODE45 in MATLAB Script
    % Set initial conditions
    all_X_mat(1,vehicle_i)   = 0;
    all_Y_mat(1,vehicle_i)   = 0;
    all_phi_mat(1,vehicle_i) = 0;
    all_r_mat(1,vehicle_i)   = 0;
    all_V_mat(1,vehicle_i)   = 0;
    
    counter = 1; % Initialize the counter
    for time = deltaT:deltaT:TotalTime
        delta_f = (1-1*(0<time-Period))*...
                  (pi/180)*steering_amplitude_degrees*sin((2*pi/Period)*time);
        
        input_states = [all_V_mat(counter,vehicle_i), ...
                        all_r_mat(counter,vehicle_i)];
        % estimate slip-angles
        alpha = fcn_slipAngles(U, input_states(1), input_states(2), ...
                               delta_f, vehicle(vehicle_i));
        
        % estimate lateral forces
        Fy = fcn_lateralForces(alpha, vehicle(vehicle_i));
        
        % estimate lateral velocity and yaw rate
        [~, y] = ode45(@(t,y) fcn_lateralDynamics(t, y, U, Fy, ...
                       vehicle(vehicle_i)), [0 deltaT], input_states);
        counter = counter+1;
        all_V_mat(counter,vehicle_i) = y(end,1);
        all_r_mat(counter,vehicle_i) = y(end,2);
        
        input_pose = [all_X_mat(counter-1,vehicle_i), ...
                      all_Y_mat(counter-1,vehicle_i), ...
                      all_phi_mat(counter-1,vehicle_i)];
        states = [all_V_mat(counter,vehicle_i), ...
                  all_r_mat(counter,vehicle_i)];
        % estimate global position
        [~, y] = ode45(@(t,y) fcn_body2globalCoordinates(t, y, U, ...
                       states(1), states(2)), [0 deltaT], input_pose);
        all_X_mat(counter,vehicle_i) = y(end,1);
        all_Y_mat(counter,vehicle_i) = y(end,2);
        all_phi_mat(counter,vehicle_i) = y(end,3);
    end
end

%% Plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____  _       _   _   _             
%  |  __ \| |     | | | | (_)            
%  | |__) | | ___ | |_| |_ _ _ __   __ _ 
%  |  ___/| |/ _ \| __| __| | '_ \ / _` |
%  | |    | | (_) | |_| |_| | | | | (_| |
%  |_|    |_|\___/ \__|\__|_|_| |_|\__, |
%                                   __/ |
%                                  |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start the plotting. Both the plots for each situation should agree with 
% each other (otherwise the equations are not consistent)

ind = (2==t);
% Loop through each of the vehicles
for vehicle_i=1:N_vehicles
    % Plot the yaw rate
    h1 = figure(99);
    hold on;
    set(h1,'Name','Yawrate')
    plot(t,all_r_mat(:,vehicle_i),'ro',...
         t,all_r_sim(:,vehicle_i),'b');
    grid on;
    text(2, all_r_sim(ind,vehicle_i), ['vehicle' num2str(vehicle_i)]);
    legend('MATLAB','Simulink');
    xlabel('Time (sec)'); 
    ylabel('Yawrate (rad/sec)');
    title('Yawrate');
    
    % Plot the lateral velocity
    h2 = figure(88);
    hold on;
    set(h2,'Name','LatVel')
    plot(t,all_V_mat(:,vehicle_i),'ro',...
         t,all_V_sim(:,vehicle_i),'b');
    grid on;
    text(2, all_V_sim(ind,vehicle_i), ['vehicle' num2str(vehicle_i)]);
    legend('MATLAB','Simulink');
    xlabel('Time (sec)'); 
    ylabel('Lateral Velocity (m/sec)');
    title('Lateral Velocity');
    
    % The XY Plots
    h3 = figure(77);
    hold on;
    set(h2,'Name','XYposition')
    plot(all_X_mat(:,vehicle_i),all_Y_mat(:,vehicle_i),'ro',...
         all_X_sim(:,vehicle_i),all_Y_sim(:,vehicle_i),'b');
    grid on;
    text(all_X_sim(ind,vehicle_i), all_Y_sim(ind,vehicle_i), ...
        ['vehicle' num2str(vehicle_i)]);
    legend('MATLAB','Simulink');
    xlabel('X position [m]');
    ylabel('Y position [m]');
    title('Position');
end

% The Steering Input
h4 = figure(66);
hold on;
set(h4,'Name','Steering Input')
plot(t,df,'b');
grid on;
xlabel('Time (sec)'); 
ylabel('Steering input (rad)');
title('Steering Input');