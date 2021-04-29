function dydt = fcn_lateralDynamics(~, y, U, Fy, vehicle)
%% fcn_lateralDynamics
%   This function calculates lateral acceleration and yaw acceleration.
%
% FORMAT:
%
%   dy = fcn_lateralDynamics(~, y, U, Fy, vehicle)
%
% INPUTS:
%
%   y: A 2x1 vector of lateral velocity and yaw rate
%   U: Longitudinal velocity
%   Fy: A 2x1 vector containing lateral-forces at front and rear wheels
%   vehicle: MATLAB structure containing vehicle properties
%
% OUTPUTS:
%
%   dydt: A 2x1 vector of lateral acceleration and yaw acceleration
%
% This function was written on 2021_04_28 by Satya Prasad
% Questions or comments? szm888@psu.edu
%
% TODO:
% 1. Add examples

flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1, 'STARTING function: %s, in file: %s\n', st(1).name, st(1).file);
end

%% Check input arguments
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
if flag_check_inputs
    % Are there the right number of inputs?
    if 5 ~= nargin
        error('Incorrect number of input arguments.')
    end
    
    % Check the 'y' input
    if (2~=length(y)) || ~isnumeric(y)
        error('States must be a 2x1 vector.');
    end
    
    % Check the 'U' input
    if (1~=length(U)) || ~isnumeric(U) || (0>U)
        error('Longitudinal velocity (U) must be a non-negative number.');
    end
    
    % Check the 'Fy' input
    if (2~=size(Fy,1)) || ~isnumeric(Fy)
        error('Lateral forces (Fy) must be a 2x1 vector.');
    end
end

%% Calculate Lateral acceleration and Yaw acceleration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = y(2);
dVdt = (1/vehicle.m)*(Fy(1)+Fy(2))-r*U;
drdt = (1/vehicle.Iz)*(vehicle.a*Fy(1)-vehicle.b*Fy(2));
dydt = [dVdt; drdt];

%% Any debugging?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _                 
%  |  __ \     | |                
%  | |  | | ___| |__  _   _  __ _ 
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_debug
    fprintf(1, 'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end

end