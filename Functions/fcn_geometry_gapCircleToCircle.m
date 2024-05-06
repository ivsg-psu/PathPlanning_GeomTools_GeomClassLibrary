function gap_between_circles = fcn_geometry_gapCircleToCircle(circle1_radius, circle2_radius, circle2_center_xy, flag_circle2_is_inside, varargin)
%% fcn_geometry_gapCircleToCircle
% Calculates the gap distance from one circle at a position of [0 radius1]
% to another circle of arbitrary XY center location. The gap distance
% is the distance from one edge of the circle to another and is defined
% depending on whether one circle is meant to be completely inside or
% completely outside. By default, the distance assumes completely inside.
%
% If one circle is meant to be completely inside, the gap is positive if
% there is no overlap of one circle to another, and the gap is a measure of
% the distance from the edge of the inner circle to the edge of the outer
% circle. If the circles overlap or one is outside the other, the gap is
% negative.
%
% If one circle is meant to be completely outside the other, the gap is
% positive if there is no overlap of one circle to another, and the gap is
% again the measure of the distance from the edge of one circle to the edge
% of another circle. If the circles overlap or one circle is inside the
% other, the gap is negative.
%
% Format:
% gap_between_circles = fcn_geometry_gapCircleToCircle(circle1_radius, circle2_radius, circle2_center_xy, flag_circle2_is_inside,  (fig_num))
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
%      flag_circle2_is_inside: a flag to indicate which gap to check. Set
%      flag_circle2_is_inside = 1 if one circle is meant to be inside the
%      other (default), or flag_circle2_is_inside = -1, if one circle is
%      meant to be outside the other.
%
%      (OPTIONAL INPUTS)
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
% See the script: script_test_fcn_geometry_gapCircleToCircle
% for a full test suite.
%
% This function was written on 2024_05_06 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_05_06 - S. Brennan
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
        narginchk(4,5);

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

if 1==flag_circle2_is_inside
    % Small circle must be completely inside the large circle
    gap_between_circles = larger_radius - center_to_center_distance_between_circles - smaller_radius;
else
    % Small circle must be completely outside the large circle
    gap_between_circles = center_to_center_distance_between_circles - larger_radius - smaller_radius;
end


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

    title(sprintf('Gap is: %.4f',gap_between_circles));
    xlabel('X [meters]');
    ylabel('Y [meters]')

    % Plot the circles
    fcn_geometry_plotCircle(circle1_center_xy, circle1_radius,'g-');
    fcn_geometry_plotCircle(circle2_center_xy, circle2_radius,'r-');


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


