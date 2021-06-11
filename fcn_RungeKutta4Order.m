function [final_time, final_states] = fcn_RungeKutta4Order(input_function, ...
                                      initial_states, initial_time, ...
                                      time_interval)
%% fcn_RungeKutta4Order
%   This function calculates 'final_states' by integrating the
%   'input_function' using Runge-Kutta 4th Order and using
%   'initial_states', 'initial_time', and 'time_interval'.
%
% FORMAT:
%
%   [final_time, final_states] = fcn_RungeKutta4Order(input_function, ...
%                                initial_states, initial_time, ...
%                                time_interval)
%
% INPUTS:
%
%   input_function: A function handle.
%   initial_states: Input states. The size of output states will be same as
%   the input states.
%   initial_time: Time related to 'initial_states'.
%   time_interval: Time interval.
%
% OUTPUTS:
%
%   final_time: Final time is sum of initial time and time interval.
%   final_states: Output states.
%
% This function was written on 2021_06_10 by Satya Prasad
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
    if 4 ~= nargin
        error('Incorrect number of input arguments.')
    end
    
    % Check the 'input_function' input
    if ~isa(input_function,'function_handle')
        error('Input function must be a function handle.');
    end
    
    % Check the 'initial_time' input
    if ~isreal(initial_time) || ~isnumeric(initial_time) || ...
            1~=numel(initial_time) || 0>initial_time
        error('Initial time must be a non-negative number.');
    end
    
    % Check the 'time_interval' input
    if ~isreal(time_interval) || ~isnumeric(time_interval) || ...
            1~=numel(time_interval) || 0>=time_interval
        error('Time interval must be a positive number.');
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
final_time = initial_time+time_interval;  % find the final time

k1 = input_function(initial_time, initial_states);
k2 = input_function(initial_time+time_interval/2, ...
                    initial_states+(time_interval/2)*k1);
k3 = input_function(initial_time+time_interval/2, ...
                    initial_states+(time_interval/2)*k2);
k4 = input_function(final_time, initial_states+time_interval*k3);
final_states = initial_states + (time_interval/6)*(k1+2*k2+2*k3+k4);

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