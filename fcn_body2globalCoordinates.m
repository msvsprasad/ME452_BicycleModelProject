function dydt = fcn_body2globalCoordinates(~, y, U, V, r)
%% fcn_body2globalCoordinates
%   This function calculates velocites in global coordinates.
%
% FORMAT:
%
%   dydt = fcn_body2globalCoordinates(~, y, U, V, r)
%
% INPUTS:
%
%   y: A 3x1 vector of global pose
%   U: Longitudinal velocity
%   V: Lateral velocity
%   r: Yaw rate
%
% OUTPUTS:
%
%   dydt: A 2x1 vector of velocities in global coordinates
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
    if (3~=length(y)) || ~isnumeric(y)
        error('Global pose must be a 3x1 vector.');
    end
    
    % Check the 'U' input
    if (1~=length(U)) || ~isnumeric(U) || (0>U)
        error('Longitudinal velocity (U) must be a non-negative number.');
    end
    
    % Check the 'V' input
    if (1~=length(V)) || ~isnumeric(V)
        error('Lateral velocity (V) must be a scalar.');
    end
    
    % Check the 'r' input
    if (1~=length(r)) || ~isnumeric(r)
        error('Yaw rate (r) must be a scalar.');
    end
end

%% Calculate velocities in Global coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Phi = y(3);

dXdt   = U*cos(Phi)-V*sin(Phi);
dYdt   = U*sin(Phi)+V*cos(Phi);
dPhidt = r;
dydt   = [dXdt; dYdt; dPhidt];

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