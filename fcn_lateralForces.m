function Fy = fcn_lateralForces(alpha, vehicle)
%% fcn_lateralForces
%   This function computes lateral-forces at front and rear wheels. It
%   assumes a linear model.
%
% FORMAT:
%
%   Fy = fcn_lateralForces(alpha, vehicle)
%
% INPUTS:
%
%   alpha: A 2x1 vector containing slip-angles of front and rear wheels
%   vehicle: MATLAB structure containing vehicle properties
%
% OUTPUTS:
%
%   Fy: A 2x1 vector containing lateral-forces at front and rear wheels
%
% This function was written on 2021_04_28 by Satya Prasad
% Questions or comments? szm888@psu.edu
%
% TODO:
% 1. Add examples
% 2. Add brush model

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
    if 2 ~= nargin
        error('Incorrect number of input arguments.')
    end
    
    % Check the 'alpha' input
    if (2~=size(alpha,1)) || ~isnumeric(alpha)
        error('Slip-angles (alphs) must be a 2x1 vector.');
    end
end

%% Calculate Lateral Tire Forces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fyf = vehicle.Caf*alpha(1);
Fyr = vehicle.Car*alpha(2);
Fy  = [Fyf; Fyr];

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