function test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma, varargin)
% fcn_geometry_fillCircleTestPoints
% given N points, with N>=2, creates a set of M points per unit distance
% aroundt he perimeter of a circle with points randomly distributed
% radially, with variance sigma.
%
% [test_points] = fcn_geometry_fillCircleTestPoints(seed_points, M, sigma)
%
% INPUTS:
%      circle_center: a 1x2 vector denoting the [X Y] location of the
%      center of the circle
%
%      circle_radius: a 1x1 vector denoting the radius of the circle
%
%      M: the number of test points to generate per unit
%      distance around the circumference
%
%      sigma: the standard deviation in points in the radial direction
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.     
%
% OUTPUTS:
%
%      test_points: a list of test points used to test regression fitting
% 
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_fillCircleTestPoints
% for a full test suite.
%
% This function was written on 2024_01_08 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2024_01_08 by S. Brennan
% -- wrote the code using fcn_geometry_fillArcTestPoints as starter

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==5 && isequal(varargin{end},-1))
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
    debug_fig_num = 34838;
    figure(debug_fig_num);
    clf;
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

if (0==flag_max_speed)
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(4,5);

        % Check the circle_center input to be [1x2] vector
        fcn_DebugTools_checkInputsToFunctions(...
            circle_center, '2column_of_numbers',[1 1]);

        % Check the circle_radius input to be [1x1] vector
        fcn_DebugTools_checkInputsToFunctions(...
            circle_radius, '1column_of_numbers',[1 1]);

        % Check the M input to be [1x1] vector
        fcn_DebugTools_checkInputsToFunctions(...
            M, '1column_of_numbers',[1 1]);

        % Check the sigma input to be [1x1] vector
        fcn_DebugTools_checkInputsToFunctions(...
            sigma, '1column_of_numbers',[1 1]);
    end
end

% Does user want to show the plots?
flag_do_plots = 0;
if (5 == nargin) && (0==flag_max_speed)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end


%% Solve for the circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the number points
projection_distances = (0:(1/M):2*pi*circle_radius)';
angle_increment = 2*pi/(length(projection_distances)-1);



% Convert to angles, and then to points
angles = (0:angle_increment:(2*pi-angle_increment))';
N_points = length(angles);

unit_radial_vectors = [cos(angles) sin(angles)];
radii = randn(N_points,1)*sigma + ones(N_points,1)*circle_radius;
perturbed_points = ones(N_points,1)*circle_center + radii.*unit_radial_vectors;

% Convert to test points
test_points = perturbed_points;




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
    figure(fig_num);
    axis equal;
    hold on;
    grid on;

    % Plot the input circle
    plot(circle_center(1,1), circle_center(1,2), 'r+','MarkerSize',30);
    fcn_geometry_plotCircle(circle_center,circle_radius,'b-', fig_num);

    % Plot the results
    plot(test_points(:,1), test_points(:,2), 'k.','MarkerSize',20);

    % Make axis slightly larger
    temp = axis;
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

    
end % Ends check if plotting

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends main function

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







