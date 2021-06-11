function alpha = fcn_slipAngles(U, V, r, delta_f, vehicle)
%% fcn_slipAngles
%   This function computes slip-angles for front and rear wheels. It
%   assumes a bicycle model.
%
% FORMAT:
%
%   alpha = fcn_slipAngles(U, V, r, delta_f, vehicle)
%
% INPUTS:
%
%   U: Longitudinal velocity [m/s]
%   V: Lateral velocity [m/s]
%   r: Yaw rate [rad/s]
%   delta_f: Steering angle (FRONT) [rad]
%   vehicle: MATLAB structure containing vehicle properties
%
% OUTPUTS:
%
%   alpha: A 2x1 vector containing slip-angles of front and rear wheels [rad]
%
% This function was written on 2021_04_28 by Satya Prasad
% Questions or comments? szm888@psu.edu
%

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
    
    % Check the 'U' input
    if ~isreal(U) || ~isnumeric(U) || 1~=numel(U) || 0>U
        error('Longitudinal velocity (U) must be a non-negative number.');
    end
    
    % Check the 'V' input
    if ~isreal(V) || ~isnumeric(V) || 1~=numel(V)
        error('Lateral velocity (V) must be a real number.');
    end
    
    % Check the 'r' input
    if ~isreal(r) || ~isnumeric(r) || 1~=numel(r)
        error('Yaw rate (r) must be a real number.');
    end
    
    % Check the 'delta_f' input
    if ~isreal(delta_f) || ~isnumeric(delta_f) || 1~=numel(delta_f)
        error('Front steering angle (delta_f) must be a real number.');
    end
end

%% Calculate Slip-Angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_f = (V+vehicle.a*r)/U - delta_f;
alpha_r = (V-vehicle.b*r)/U;
alpha   = [alpha_f; alpha_r];

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