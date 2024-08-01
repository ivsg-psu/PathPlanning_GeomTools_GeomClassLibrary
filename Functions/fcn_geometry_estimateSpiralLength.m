function [estimated_spiralLength, estimated_spiralStartAngle, estimated_spiralEndAngle, angle_larger_to_smaller] = ...
    fcn_geometry_estimateSpiralLength(circle1_parameters, circle2_parameters, varargin)
%% fcn_geometry_estimateSpiralLength
% Calculates an estimated spiral length that joins a second circle to a
% unit circle.
%
% Format:
% estimated_length = fcn_geometry_estimateSpiralLength(circle2_radius, offset, (fig_num))
%
% INPUTS:
%
%      circle2_radius: the normalized radius of circle2, e.g. r2/r1
%
%      offset: the distance bween the circle edges at their closest point.
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      estimated_length: the estimated length of the spiral arc. 
%
%      estimated_spiralStartAngle, estimated_spiralEndAngle: the estimated
%      start and end angles of the connecting arc, relative to the vector
%      connecting the larger circle center to the smaller circle center
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_estimateSpiralLength
% for a full test suite.
%
% This function was written on 2024_07_30 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% 2024_07_30 - S. Brennan
% -- wrote the code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==3 && isequal(varargin{end},-1))
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS");
    MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG = getenv("MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG); 
        flag_check_inputs  = str2double(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

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
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(2,3);

        % % Check the points input to be length greater than or equal to 2
        % fcn_DebugTools_checkInputsToFunctions(...
        %     points, '2column_of_numbers',[2 3]);
        %
        % % Check the transverse_tolerance input is a positive single number
        % fcn_DebugTools_checkInputsToFunctions(transverse_tolerance, 'positive_1column_of_numbers',1);
        %
        % % Check the station_tolerance input is a positive single number
        % if ~isempty(station_tolerance)
        %     fcn_DebugTools_checkInputsToFunctions(station_tolerance, 'positive_1column_of_numbers',1);
        % end
    end
end

% Does user want to specify fig_num?
flag_do_plots = 0;
if 0==flag_max_speed && 3<=nargin
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp; %#ok<NASGU>
        flag_do_plots = 1;
    end
end


%% Solve for the line fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
center1 = circle1_parameters(1,1:2);
radius1 = circle1_parameters(1,3);

center2 = circle2_parameters(1,1:2);
radius2 = circle2_parameters(1,3);

vector_center1_to_center2 = center2-center1;
distance_center1_to_center2 = real(sum(vector_center1_to_center2.^2,2).^0.5);

% Determine which circle is the larger one
if radius1>radius2
    larger_center = center1;
    normalizing_distance = radius1;
    circle2_radius = radius2/normalizing_distance;
    offset = (radius1 - radius2 - distance_center1_to_center2)/normalizing_distance;
    larger_to_smaller_vector = vector_center1_to_center2;
else
    larger_center = center2;
    normalizing_distance = radius2;
    circle2_radius = radius1/normalizing_distance;
    offset = (radius2 - radius1 - distance_center1_to_center2)/normalizing_distance;
    larger_to_smaller_vector = -vector_center1_to_center2;
end

assert(offset>0)

%% Make sure offsets and radii are within the fit range
flag_inputs_are_in_good_range = 1;

min_offset = 0.01/1000;
max_offset = 0.1;
if offset<min_offset || offset>max_offset
    warning('Input offset is out of range');
    flag_inputs_are_in_good_range = 0;
end

min_fraction_of_radial_space_for_r2 = 1 - 0.2;
max_fraction_of_radial_space_for_r2 = 1 - 0.01/1000;
max_circle2_radii = (1-offset)*(max_fraction_of_radial_space_for_r2);
min_circle2_radii = (1-offset)*(min_fraction_of_radial_space_for_r2);

if circle2_radius<min_circle2_radii || circle2_radius>max_circle2_radii
    warning('Input smaller radii is out of range');
    flag_inputs_are_in_good_range = 0;
end

%% Perform calculations

if 1==flag_inputs_are_in_good_range
    estimated_spiralLength     = fcn_INTERNAL_estimateSpiralLength(circle2_radius, offset)*normalizing_distance;
    estimated_spiralStartAngle = -fcn_INTERNAL_estimateSpiralStartAngle(circle2_radius, offset);
    estimated_spiralEndAngle   = estimated_spiralStartAngle + ((1/radius1 + 1/radius2)/2)*estimated_spiralLength;
else
    % fill the values with NaN
    estimated_spiralStartAngle   = nan;
    estimated_spiralLength       = nan;
    estimated_spiralEndAngle     = nan;
end
angle_larger_to_smaller    = atan2(larger_to_smaller_vector(2),larger_to_smaller_vector(1));

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
    axis equal;

    title('Joining circles with spiral');
    xlabel('X [meters]');
    ylabel('Y [meters]')

    % Calculate the correct values
    spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, [], -1);
    h0           = spiral_join_parameters(1,3);
    spiralLength = spiral_join_parameters(1,4);
    % x0           = spiral_join_parameters(1,3);
    % y0           = spiral_join_parameters(1,4);
    K0           = spiral_join_parameters(1,5);
    Kf           = spiral_join_parameters(1,6);
    analytical_end_angle   = h0 + (Kf+K0)/2*spiralLength;


    % Plot the inputs
    fcn_geometry_plotGeometry('circle',circle1_parameters, [], sprintf(' ''LineWidth'',4'));
    fcn_geometry_plotGeometry('circle',circle2_parameters, [], sprintf(' ''LineWidth'',4'));

    % Plot the spiral result
    fcn_geometry_plotGeometry('spiral',spiral_join_parameters, [], sprintf(' ''LineWidth'',3'));

    % Plot the join angle
    points_join_angle = [larger_center; larger_center + normalizing_distance*[cos(angle_larger_to_smaller) sin(angle_larger_to_smaller)]];
    plot(points_join_angle(:,1),points_join_angle(:,2),'k-');

    % Plot the true and estimated start angle
    points_true_start_angle = [larger_center; larger_center + normalizing_distance*[cos(h0-pi/2) sin(h0-pi/2)]];
    plot(points_true_start_angle(:,1),points_true_start_angle(:,2),'-','LineWidth',5,'Color',[0 1 0]);
    points_est_start_angle = [larger_center; larger_center + normalizing_distance*[cos(estimated_spiralStartAngle-pi/2) sin(estimated_spiralStartAngle-pi/2)]];
    plot(points_est_start_angle(:,1),points_est_start_angle(:,2),'-','LineWidth',2,'Color',[0 0.5 0]);
    
    % Plot the true and estimated end angle
    points_true_end_angle = [larger_center; larger_center + normalizing_distance*[cos(analytical_end_angle-pi/2) sin(analytical_end_angle-pi/2)]];
    plot(points_true_end_angle(:,1),points_true_end_angle(:,2),'-','LineWidth',5,'Color',[1 0 0]);
    points_est_end_angle = [larger_center; larger_center + normalizing_distance*[cos(estimated_spiralEndAngle-pi/2) sin(estimated_spiralEndAngle-pi/2)]];
    plot(points_est_end_angle(:,1),points_est_end_angle(:,2),'-','LineWidth',2,'Color',[0.5 0 0.5]);
    




    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end

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

%% fcn_INTERNAL_estimateLengthSlopeFromOffset
function slope = fcn_INTERNAL_estimateLengthSlopeFromOffset(offset)
% This function estimates the slope from the offset

% Define the x-data used for fitting
xdata = log(-log(offset));

% xdata = log(-log(offsets_to_try));
% ydata = log(0.5-slopes);


% Uses equation of form
% y = x2*a + x*b + c --> Y = [X^2 X 1]*[a b c]'
X_vec = [xdata.^3 xdata.^2 xdata ones(length(xdata(:,1)),1)];

fit_params = [-0.023713568015240
   0.205326944457505
  -0.596356120454393
   0.588196591187278];

y_pred = X_vec*fit_params;

slope = 0.5-exp(y_pred);
end % Ends fcn_INTERNAL_estimateLengthSlopeFromOffset

%% fcn_INTERNAL_estimateLengthInterceptFromOffset
function lengthIntercept = fcn_INTERNAL_estimateLengthInterceptFromOffset(offset)
% Define the x-data used for fitting
xdata = log(offset);

% xdata = log(offsets_to_try);
% ydata = intercepts;


% Uses equation of form
% y = x2*a + x*b + c --> Y = [X^2 X 1]*[a b c]'
X_vec = [xdata.^3 xdata.^2 xdata ones(length(xdata(:,1)),1)];

% These were fit from earlier
fit_params = [
    -0.000648129505156
  -0.018419651621482
   0.315399921979639
   0.849677809445654
   ];

y_pred = X_vec*fit_params;

lengthIntercept = y_pred;
end % Ends fcn_INTERNAL_estimateLengthInterceptFromOffset


%% fcn_INTERNAL_estimateSpiralLength
function estimated_length = fcn_INTERNAL_estimateSpiralLength(circle2_radius, offset)

% x_data = log(1-circle2_radii);
% y_data = log(lengths);

x_data = log(1-circle2_radius);

% Get the slope and intercept
fit_slope = fcn_INTERNAL_estimateLengthSlopeFromOffset(offset);
fit_intercept = fcn_INTERNAL_estimateLengthInterceptFromOffset(offset);

estimated_length = exp(fit_slope*x_data + fit_intercept);
end % Ends fcn_INTERNAL_estimateSpiralLength


%% fcn_INTERNAL_estimateStartAngleSlopeFromOffset
function startAngleSlope = fcn_INTERNAL_estimateStartAngleSlopeFromOffset(offset)
% This function estimates the slope from the offset

    % Define the x-data used for fitting
xdata = log(-log(offset));

% Uses equation of form
% y = x2*a + x*b + c --> Y = [X^2 X 1]*[a b c]'
X_vec = [xdata.^3 xdata.^2 xdata ones(length(xdata(:,1)),1)];

fit_params = [ 
    0.027224755607716
  -0.232029761011433
   0.662948799520068
  -1.142906989473347];

y_pred = X_vec*fit_params;

startAngleSlope = y_pred;
end % Ends fcn_INTERNAL_estimateStartAngleSlopeFromOffset

%% fcn_INTERNAL_estimateStartAngleInterceptFromOffset
function startAngleIntercept = fcn_INTERNAL_estimateStartAngleInterceptFromOffset(offset)
% Define the x-data used for fitting
xdata = log(offset);

% xdata = log(offsets_to_try);
% ydata = startAngle_intercepts;


% Uses equation of form
% y = x2*a + x*b + c --> Y = [X^2 X 1]*[a b c]'
X_vec = [xdata.^3 xdata.^2 xdata ones(length(xdata(:,1)),1)];

% These were fit from earlier
fit_params = [
    -0.000805873132752
  -0.022273590488785
   0.285242206383266
   0.081271471353003
   ];

y_pred = X_vec*fit_params;

startAngleIntercept = y_pred;
end % Ends fcn_INTERNAL_estimateStartAngleInterceptFromOffset


%% fcn_INTERNAL_estimateSpiralStartAngle
function estimated_startAngle = fcn_INTERNAL_estimateSpiralStartAngle(circle2_radius, offset)

% x_data = log(1-circle2_radii);
% y_data = log(start_angles);

x_data = log(1-circle2_radius);

% Get the slope and intercept
fit_slope = fcn_INTERNAL_estimateStartAngleSlopeFromOffset(offset);
fit_intercept = fcn_INTERNAL_estimateStartAngleInterceptFromOffset(offset);

estimated_startAngle = exp(fit_slope*x_data + fit_intercept);
end % Ends fcn_INTERNAL_estimateSpiralStartAngle

