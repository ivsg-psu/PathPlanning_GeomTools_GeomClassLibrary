function [spiral_join_parameters, space_between_circles] = fcn_geometry_spiralFromCircleToCircle(circle1_radius, circle2_radius, circle2_center_xy, varargin)
%% fcn_geometry_spiralFromCircleToCircle
% Calculates the spiral parameters that join one circle to another. The
% spiral is assumed to leave the first circle in a counter-clockwise
% (positively-changing angle) direction at a heading of 0, and the
% initial position of the spiral is at the origin. Thus, the first circle
% is created at a posiion of [0 circle1_radius] in the XY plane. The second
% circle can have either a counter-clockwise or clockwise connection to the
% spiral.
%
% Format:
% spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_radius, circle2_radius, circle2_center_xy, (flag_circle2_is_counterclockwise),  (fig_num))
%
% INPUTS:
%
%      circle1_radius: the radius of the first circle.
%
%      circle2_radius: the radius of the second circle.
%
%      circle2_center_xy: the center XY point of the second circle in [X Y]
%      format 
%
%
%      (OPTIONAL INPUTS)
%
%      flag_circle2_is_counterclockwise: a flag to indicate arc direction.
%      Set flag_circle2_is_counterclockwise = 1 if the connection to circle 2
%      would cause a counter-clockwise path (default), or 
%      flag_circle2_is_counterclockwise = -1, if the circle 2 would cause
%      a clockwise path.
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      spiral_join_parameters: the parameter set describing the
%      spiral segment geometry that joins the circles. See
%      fcn_geometry_fillEmptyDomainStructure for details, specifically the
%      structure for 'spiral'. If a spiral is not possible, then nan values
%      are returned for all parameters.
%
%      space_between_circles: the amount of space between the circles if a
%      spiral is to be joined between them. Positive values will generate
%      feasible spirals, negative values are not feasible. This
%      space_between_circles is a useful measure to determine how "close"
%      in physical distance a configuration that does NOT work would be to
%      a configuration that does work.
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_spiralFromCircleToCircle
% for a full test suite.
%
% This function was written on 2024_04_24 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_04_24 - S. Brennan
% -- wrote the code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==5 && isequal(varargin{end},-1))
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
        narginchk(3,5);

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

% Does user want to specify best_fit_domain_box_projection_distance?
flag_circle2_is_counterclockwise = 1; 
if (4<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        flag_circle2_is_counterclockwise = temp; 
    end
end

% Does user want to specify fig_num?
flag_do_plots = 0;
if 5<= nargin && 0==flag_max_speed
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
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

circle1_center_xy = [0 circle1_radius];

% Calculate the space that the spiral must fit within for feasibility. The
% smaller circle MUST fit within the larger circle - otherwise, no feasible
% spiral is possble.
center_to_center_distance_between_circles = sum((circle2_center_xy - circle1_center_xy).^2,2).^0.5;

if circle1_radius>circle2_radius
    larger_radius = circle1_radius;
    smaller_radius = circle2_radius;
else
    larger_radius = circle2_radius;
    smaller_radius = circle1_radius;
end

if 1==flag_circle2_is_counterclockwise
    % Small circle must be completely inside the large circle
    space_between_circles = larger_radius - center_to_center_distance_between_circles - smaller_radius;
else
    % Small circle must be completely outside the large circle
    space_between_circles = center_to_center_distance_between_circles - larger_radius - smaller_radius;
end

if space_between_circles>0
    spiralLength = fcn_INTERNAL_findLengthFromOffset(circle1_radius, circle2_radius*flag_circle2_is_counterclockwise, center_to_center_distance_between_circles);

    % Check results?
    if 1==flag_do_debug
        % Set up station coordinates
        s  = (0:0.01:1)'*spiralLength;

        % Call the function fcn_geometry_extractXYfromSTSpiral to predict the
        % spiral and calculate the offsets, plotting the results in
        % figure 1234

        figure(1234);
        clf;
        hold on;
        grid on;
        axis equal;

        h0 = 0;
        x0 = 0;
        y0 = 0;
        K0 = 1/circle1_radius;
        Kf = 1/(circle2_radius*flag_circle2_is_counterclockwise);

        [x_spiral,y_spiral] = fcn_geometry_extractXYfromSTSpiral(s,spiralLength,h0,x0,y0,K0,Kf,(1234));


        % Find the center of the circle tangent at the end of the spiral
        % Find the unit vector (need to do this analytically!)
        s_tangent = [0.99999999 1]'*spiralLength;
        [x_tangent,y_tangent] = fcn_geometry_extractXYfromSTSpiral(s_tangent,spiralLength,h0,x0,y0,K0,Kf);
        unit_tangent = fcn_geometry_calcUnitVector([diff(x_tangent) diff(y_tangent)]);
        approximate_end_angle = atan2(unit_tangent(2),unit_tangent(1));
        analytical_end_angle   = h0 + (Kf-K0)*spiralLength/2 + K0*spiralLength;
        disp([approximate_end_angle analytical_end_angle]);
        unit_tangent           = [cos(analytical_end_angle) sin(analytical_end_angle)];
        unit_orthogonal = unit_tangent*[0 1; -1 0];
        calculated_circle2_center_xy = circle2_radius*flag_circle2_is_counterclockwise*unit_orthogonal + [x_spiral(end) y_spiral(end)];

        % Plot the circle's centers
        plot(0,circle1_radius,'b+');
        plot(circle2_center_xy(:,1),circle2_center_xy(:,2),'r+');
        plot(calculated_circle2_center_xy(:,1),calculated_circle2_center_xy(:,2),'m+');

        % Plot the circles
        fcn_geometry_plotCircle([0 circle1_radius], circle1_radius,'b-',(1234));
        fcn_geometry_plotCircle([0 circle1_radius], center_to_center_distance_between_circles,'b--',(1234));

        fcn_geometry_plotCircle(calculated_circle2_center_xy, circle2_radius,'r-',(1234));
        fcn_geometry_plotCircle(circle2_center_xy, circle2_radius,'r-',(1234));
        


    end % Ends plotting


    % Find the angles between the original centers
    vector_from_center1_to_center2 = circle2_center_xy - circle1_center_xy;
    unit_vector_from_center1_to_center2 = fcn_geometry_calcUnitVector(vector_from_center1_to_center2);
   
    % Find the angle to correct the spiral so that it lands exactly on top
    % of the original circle1. This requires us to recalculate the initial
    % heading the attachment point of the spiral [x0,y0] on circle1
    h0 = 0;
    x0 = 0;
    y0 = 0;
    K0 = 1/circle1_radius;
    Kf = 1/(circle2_radius*flag_circle2_is_counterclockwise);
    [x_spiral,y_spiral] = fcn_geometry_extractXYfromSTSpiral(spiralLength,spiralLength,h0,x0,y0,K0,Kf,(1234));

    analytical_end_angle   = h0 + (Kf-K0)*spiralLength/2 + K0*spiralLength;
    unit_tangent           = [cos(analytical_end_angle) sin(analytical_end_angle)];
    unit_orthogonal = unit_tangent*[0 1; -1 0];
    calculated_circle2_center_xy = circle2_radius*flag_circle2_is_counterclockwise*unit_orthogonal + [x_spiral(end) y_spiral(end)];
    vector_from_center1_to_calculated_center2 = calculated_circle2_center_xy - circle1_center_xy;
    unit_vector_from_center1_to_calculated_center2 = fcn_geometry_calcUnitVector(vector_from_center1_to_calculated_center2);
    cross_product = cross([unit_vector_from_center1_to_center2 0],[unit_vector_from_center1_to_calculated_center2 0]);
    dot_product   = dot(unit_vector_from_center1_to_center2,unit_vector_from_center1_to_calculated_center2);
    angle_magnitude = acos(dot_product);
    rotation_angle = sign(cross_product(3))*angle_magnitude;
    circle1_start_xy = circle1_center_xy + circle1_radius*[cos(-90*pi/180 - rotation_angle) sin(-90*pi/180 - rotation_angle)];

    spiral_join_parameters(1,1) = spiralLength;
    spiral_join_parameters(1,2) = -rotation_angle;
    spiral_join_parameters(1,3) = circle1_start_xy(1,1);
    spiral_join_parameters(1,4) = circle1_start_xy(1,2);
    spiral_join_parameters(1,5) = K0;
    spiral_join_parameters(1,6) = Kf;

else
    spiral_join_parameters = nan(1,6);
end % Ends if statement to check if spiral is possible

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

    % Plot the circles
    fcn_geometry_plotCircle(circle1_center_xy, circle1_radius,'b-');
    fcn_geometry_plotCircle(circle2_center_xy, circle2_radius,'r-');

    % Plot the spiral result
    fcn_geometry_plotGeometry('spiral',spiral_join_parameters);


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


%% fcn_INTERNAL_findLengthFromOffset
function sprialLength = fcn_INTERNAL_findLengthFromOffset(arc1_radius, arc2_radius, center_to_center_distance)
function_to_optimize = @(x)fcn_INTERNAL_calcSpiralOffsetError(x, arc1_radius, arc2_radius, center_to_center_distance);
% options = optimset('Display','iter','PlotFcns',@optimplotfval);
% sprialLength = fminsearch(function_to_optimize,1,options);
X0 = 1;
X_solution = fminsearch(function_to_optimize,X0);
sprialLength = X_solution(1,1);

end % Ends fcn_INTERNAL_findLengthFromOffset


%% fcn_INTERNAL_calcSpiralOffsetError
function error = fcn_INTERNAL_calcSpiralOffsetError(X, circle1_radius, circle2_radius, center_to_center_distance)

% Pull out the optimization inputs
spiralLength = X;

% Fix all the other values
h0 = 0;
x0 = 0;
y0 = 0;
K0 = 1/circle1_radius;
Kf = 1/circle2_radius;

% Call the function
[x_spiral,y_spiral] = fcn_geometry_extractXYfromSTSpiral(spiralLength,spiralLength,h0,x0,y0,K0,Kf,-1);

% This is the equation for the angle of the spiral relative to horizontal
% at each position along its length. (This was confirmed via a script).
analytical_end_angle   = h0 + (Kf-K0)*spiralLength/2 + K0*spiralLength;

unit_tangent_where_spiral_joins_circle2    = [cos(analytical_end_angle) sin(analytical_end_angle)];
unit_orthogonal_where_spiral_joins_circle2 = unit_tangent_where_spiral_joins_circle2*[0 1; -1 0];

circle1_center = [0 circle1_radius];
circle2_center_from_spiral = circle2_radius*unit_orthogonal_where_spiral_joins_circle2 + [x_spiral(end) y_spiral(end)];
spiral_predicted_center_to_center_distance = sum((circle1_center-circle2_center_from_spiral).^2,2).^0.5;

flag_do_debug = 0;
if 1==flag_do_debug
    % Set up station coordinates
    s  = (0:0.01:1)'*spiralLength;

    % Call the function fcn_geometry_extractXYfromSTSpiral to predict the
    % spiral and calculate the offsets, plotting the results in
    % figure 1234

    figure(1234);
    clf;
    hold on;
    grid on;
    axis equal;

    h0 = 0;
    x0 = 0;
    y0 = 0;
    K0 = 1/circle1_radius;
    Kf = 1/circle2_radius;

    fcn_geometry_extractXYfromSTSpiral(s,spiralLength,h0,x0,y0,K0,Kf,(1234));


    % Plot the circle's centers
    plot(0,circle1_radius,'b+','MarkerSize',20);
    plot(circle2_center_from_spiral(:,1),circle2_center_from_spiral(:,2),'m+','MarkerSize',20);

    % Plot the circles
    fcn_geometry_plotCircle([0 circle1_radius], circle1_radius,'b-');
    fcn_geometry_plotCircle(circle2_center_from_spiral, circle2_radius,'b-');

    fcn_geometry_plotCircle([0 circle1_radius], center_to_center_distance,'b--');

end % Ends plotting


error = abs(center_to_center_distance - spiral_predicted_center_to_center_distance);


end % Ends fcn_INTERNAL_calcSpiralOffsetError

