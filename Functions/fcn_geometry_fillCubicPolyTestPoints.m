function [test_points, true_points] = fcn_geometry_fillCubicPolyTestPoints (coeff_x3, coeff_x2, coeff_x1, coeff_x0, x_range, M, sigma, varargin)
%% fcn_geometry_fillCubicPolyTestPoints
% 
% Given coefficients, range of x (x_range), M (no. of points per unit
% distance), sigma (standard deviation of points), this function generates
% a set of points around the perimeter of polynomial (max: cubic
% polynomial) with randomly distributed with variance sigma.
%
% FORMAT: 
%
% [test_points, true_points] = fcn_geometry_fillCubicPolyTestPoints (coeff_x3, coeff_x2, coeff_x1, coeff_x0, x_range, num_points, sigma, (fig_num)) 
% 
% INPUTS:
%
%      coeff_x3: the coefficient for x^3 term of the cubic equation
%
%      coeff_x2: the coefficient for x^2 term of the cubic equation
%
%      coeff_x: the coefficient for x term of the cubic equation
%
%      coeff_x0: the coefficient for x^ term of the cubic equation or
%      simply "constant" 
%
%      x_range: the range of x values, for example: [-5 5] or [0 2]
%
%      M: the number of test points to generate per unit
%      distance.
%
%      sigma: the standard deviation in points
%
%
% (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      test_points: a list of test points used to test regression fitting.
%      It is a Nx2 vector where N is the number of points.
%
%      true_points: the true values of the cubic polynomial equations used
%      to fill the points
%
%
% DEPENDENCIES:
%
%   (NONE)
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_fillCubicPolyTestPoints
% for a full test suite.
%
% This function was written on 2024_05_08 by Aneesh Batchu
% Questions or comments? abb6486@psu.edu or sbrennan@psu.edu

% Revision History
% 2024_05_08 - Aneesh Batchu
% -- wrote the code 



%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==8 && isequal(varargin{end},-1))
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS");
    MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG = getenv("MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS);
    end
end

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 34838; %#ok<NASGU>
else
    debug_fig_num = []; %#ok<NASGU>
end

%% check input arguments
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

if 0==flag_max_speed
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(7,8);

    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (8<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

%% Main Code starts from here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Generate x values within the specified range
x_values = (x_range(1):1/M:x_range(2))'; 

% Evaluate the cubic polynomial function to get true y values
true_y_values = coeff_x3 * x_values.^3 + coeff_x2 * x_values.^2 + coeff_x1 * x_values + coeff_x0;

% Generate random deviations (noise) with normal distribution
noise = sigma * randn(length(x_values(:,1)), 1);

% Add noise to true y values to get test points
test_y_values = true_y_values + noise;

% Combine x and y values to form test points
test_points = [x_values, test_y_values];

% Combine x and true y values to form true points
true_points = [x_values, true_y_values];

%% Plot the results (for debugging)?
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
if flag_do_plots

    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end

    hold on;
    grid on;
    % axis equal
    xlabel('X [m]')
    ylabel('Y [m]')


    % plot the true points
    plot(true_points(:,1), true_points(:,2), 'k-', 'MarkerSize',10)

    % plot the test points
    plot(test_points(:,1), test_points(:,2), 'b.', 'MarkerSize',10)

    % plot the first point of the test points
    plot(true_points(1,1), true_points(1,2), 'r.', 'MarkerSize',20)

    % plot the last point of the test points
    plot(true_points(end,1), true_points(end,2), 'r.', 'MarkerSize',20)

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end

%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§


