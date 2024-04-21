function [revised_line_parameters, revised_arc_parameters, revised_spiral_join_parameters] = fcn_geometry_alignLineToArc(line_parameters, arc_parameters, flag_arc_is_first, varargin)
%% fcn_geometry_alignLineToArc
% Revises the geometric parameters of an arc and line such that they align
% where they join. It does this by checking the offset between the two
% objects at the join location. If the offset is less than a threshold
% (default is 0.1 meter), the second geometric object is shifted to force
% alignment.
%
% Format:
% [revised_line_parameters, revised_arc_parameters]  = ...
% fcn_geometry_alignLineToArc(line_parameters, arc_parameters, flag_arc_is_first, (threshold), (fig_num))
%
% INPUTS:
%
%      line_parameters: the parameter set describing the line segment
%      geometry. See fcn_geometry_fillEmptyDomainStructure for details,
%      specifically the structure for 'Vector regression segment fit'.
%
%      arc_parameters: the parameter set describing the arc geometry. See
%      fcn_geometry_fillEmptyDomainStructure for details, specifically the
%      structure for 'Regression arc'.
%
%      flag_arc_is_first: set to 1 if the arc precedes the line, or 0 if
%      the arc is after the line.
%
%      (OPTIONAL INPUTS)
%
%      threshold: the offset, in meters, between the arc and the line such
%      that this offset is removed by shifting. If the offset is larger
%      than this, then the outputs are set to empty.If this is
%      entered as a 2x1 or 1x2, then this specifies the threshold first in
%      the transverse direction, and then in the station direction. For
%      example, an entry of [0.02 3] would have 0.02 meters threshold in
%      the transverse direction, but 3 meters threshold in the station
%      direction.
%
%      continuity_level: the level of continuity desired in the alignment.
%      Input values include 0 for C0 continuity, 1 for C1 continuity
%      (default), or 2 for C2 continuity. For an explanation of continuity,
%      see fcn_geometry_alignGeometriesInSequence
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      revised_line_parameters: the parameter set describing the line
%      segment geometry that joins the geometries. See
%      fcn_geometry_fillEmptyDomainStructure for details, specifically the
%      structure for 'Vector regression segment fit'.
%
%      revised_arc_parameters: the parameter set describing the arc
%      segment geometry that joins the geometries. See
%      fcn_geometry_fillEmptyDomainStructure for details, specifically the
%      structure for 'Regression arc'.
%
%      revised_spiral_join_parameters: the parameter set describing the
%      spiral segment geometry that joins the line and arc geometries if C2
%      continuity is specified. See fcn_geometry_fillEmptyDomainStructure
%      for details, specifically the structure for 'spiral'.
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_alignLineToArc
% for a full test suite.
%
% This function was written on 2024_04_12 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_04_12 - S. Brennan
% -- wrote the code
% 2024_04_19
% -- renamed from fcn_geometry_joinLineToArc
% -- fixed bug where calculation still works if error larger than tolerance
% -- added continuity_level input
% 2024_04_20
% -- added St conversion functions
% -- added powerful debugging plots (VERY useful - caught lots of mistakes)
% -- finished functionalizing code
% -- added C0 and C1 continuity, confirmed via script testing they work
% -- added C2 continuity and revised_spiral_join_parameters output

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==6 && isequal(varargin{end},-1))
    flag_do_debug = 0; % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % Flag to plot the results for debugging
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
    debug_fig_num = 34838;
else
    debug_fig_num = [];
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
        narginchk(3,6);

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
threshold = 0.1;
if (4<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        threshold = temp;
    end
end

% Does user want to specify continuity_level?
continuity_level = 1;
if (4<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        continuity_level = temp;
        if ~any(continuity_level == [0 1 2])
            error('The continuity_level input must be 0, 1, or 2');
        end
    end
end

% Does user want to specify fig_num?
flag_do_plots = 0;
if 6<= nargin && 0==flag_max_speed
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

%% Plot inputs?
if flag_do_debug
    figure(debug_fig_num);
    clf;

    % Plot the inputs
    subplot(3,2,1);


    fcn_geometry_plotGeometry('line',line_parameters);
    fcn_geometry_plotGeometry('arc',arc_parameters);

    temp = axis;
    %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    new_axis = [temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y];
    axis([min(new_axis) max(new_axis) min(new_axis) max(new_axis)]);

    hold on;
    grid on;
    axis equal;
    xlabel('X [meters]');
    ylabel('Y [meters]')

    title('Inputs');

    debug_axis = axis;

end

%% Rearrange parameters so line is always the 1st input, arc is 2nd
% Fix the parameters to make the line segment first, arc second, and make
% sure the line and arc point into and then out of the junction
% respectively. For situations where arc is actually the first input, this
% is fixed in later steps using a flag.
[clean_line_parameters, clean_arc_parameters] = fcn_INTERNAL_fixOrientationAndOrdering(line_parameters, arc_parameters);

if flag_do_debug
    figure(debug_fig_num);

    % Plot the cleaned inputs
    subplot(3,2,2);

    fcn_geometry_plotGeometry('line',clean_line_parameters);
    fcn_geometry_plotGeometry('arc',clean_arc_parameters);
    title('Cleaned inputs');
    axis(debug_axis);
end

%% Rotate the geometries so that the line is oriented horizontally
% This is to make the debugging MUCH easier, as it reduces permutations.
% Again, this is fixed in later steps.
[st_line_parameters, st_arc_parameters, St_transform, rotation_angle] = fcn_INTERNAL_convertParametersToStOrientation(clean_line_parameters, clean_arc_parameters);

if flag_do_debug
    figure(debug_fig_num);

    % Plot the rotated inputs
    subplot(3,2,3);

    fcn_geometry_plotGeometry('line',st_line_parameters);
    fcn_geometry_plotGeometry('arc',st_arc_parameters);
    title('Rotated into St');
    axis(debug_axis);
end


%% Check to see if arc and line intersect
intersection_point = fcn_INTERNAL_findLineArcIntersection(st_line_parameters,st_arc_parameters);


if flag_do_debug
    % Plot the intersection
    figure(debug_fig_num);
    subplot(3,2,3);
    plot(intersection_point(:,1),intersection_point(:,2),'r.','MarkerSize',20);
    axis(debug_axis);
end


%% Check how much shift is needed to connect line to arc
[delta_transverse, delta_station, desired_arc_start_point, desired_line_end_point,spiral_join_parameters] = fcn_INTERNAL_findShiftToMatchArcToLine(st_line_parameters, st_arc_parameters,continuity_level, intersection_point);
% Deltas are from desired to actual

if flag_do_debug
    % Plot the bounding box
    figure(debug_fig_num);
    subplot(3,2,3);

    % bounding_box = [0 0; 1 0; 1 1; 0 1; 0 0]*[threshold(1) 0; 0 threshold(1)];
    % plot(bounding_box(:,1),bounding_box(:,2),'r-','LineWidth',1);

    fcn_geometry_plotCircle(desired_line_end_point, threshold(1),'r-');

    plot(desired_line_end_point(:,1),desired_line_end_point(:,2),'b.','MarkerSize',30);
    plot(desired_arc_start_point(:,1),desired_arc_start_point(:,2),'c.','MarkerSize',10);
    fcn_geometry_plotGeometry('spiral',spiral_join_parameters);
    
    axis(debug_axis);
end

%% Do we join the line to the arc?
shift_distance = sum(([delta_transverse delta_station]).^2,2).^0.5;
if abs(shift_distance)<threshold
    [revised_line_parameters_St,revised_arc_parameters_St] = fcn_INTERNAL_performShift(flag_arc_is_first, st_line_parameters, st_arc_parameters,continuity_level, delta_transverse,delta_station, desired_arc_start_point, desired_line_end_point);
    revised_spiral_join_parameters_St = spiral_join_parameters;
else
    % Not possible to shift
    revised_line_parameters_St = [];
    revised_arc_parameters_St  = [];
    revised_spiral_join_parameters_St = [];
end


if flag_do_debug
    % Plot the results
    figure(debug_fig_num);
    subplot(3,2,4);

    fcn_geometry_plotGeometry('line',revised_line_parameters_St);
    fcn_geometry_plotGeometry('arc',revised_arc_parameters_St);
    fcn_geometry_plotGeometry('spiral',revised_spiral_join_parameters_St);

    title('St outputs');
    axis(debug_axis);
end

%% Rotate results out of St

if ~isempty(revised_line_parameters_St)
    [revised_line_parameters,revised_arc_parameters, revised_spiral_join_parameters] = ...
        fcn_INTERNAL_convertParametersOutOfStOrientation(...
        revised_line_parameters_St, revised_arc_parameters_St, St_transform, rotation_angle, revised_spiral_join_parameters_St);
else
    % Not possible to shift
    revised_line_parameters = [];
    revised_arc_parameters  = [];
    revised_spiral_join_parameters = [];
end

if flag_do_debug
    % Plot the results
    figure(debug_fig_num);
    subplot(3,2,5);

    fcn_geometry_plotGeometry('line',revised_line_parameters);
    fcn_geometry_plotGeometry('arc',revised_arc_parameters);
    fcn_geometry_plotGeometry('spiral',revised_spiral_join_parameters);
    
    title('St outputs');
    axis(debug_axis);
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

    subplot(1,2,1);
    hold on;
    grid on;
    axis equal;
    sgtitle('Arc joining with line segment');
    title('Original');
    xlabel('X [meters]');
    ylabel('Y [meters]')

    % Plot the inputs
    fcn_geometry_plotGeometry('segment',line_parameters);
    fcn_geometry_plotGeometry('arc',arc_parameters);

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end

    temp_axis = axis;

    subplot(1,2,2);
    hold on;
    grid on;
    axis equal;
    title('Revised');
    xlabel('X [meters]');
    ylabel('Y [meters]')

    % Plot the inputs
    fcn_geometry_plotGeometry('segment',line_parameters);
    fcn_geometry_plotGeometry('arc',arc_parameters);


    % Plot the outputs
    fcn_geometry_plotGeometry('segment',revised_line_parameters);
    fcn_geometry_plotGeometry('arc',revised_arc_parameters);
    fcn_geometry_plotGeometry('spiral',revised_spiral_join_parameters);


    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        temp_axis2 = [temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y];
        good_axis = [min([temp_axis(1) temp_axis2(1)])  max([temp_axis(2) temp_axis2(2)])  min([temp_axis(3) temp_axis2(3)])  max([temp_axis(4) temp_axis2(4)]) ];
        axis(good_axis);
    end

    % Force the subplots to have matching axes
    subplot(1,2,2);
    good_axis = axis;

    subplot(1,2,1);
    axis(good_axis)

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

%% fcn_INTERNAL_fixOrientationAndOrdering
function [clean_line_parameters, clean_arc_parameters] = fcn_INTERNAL_fixOrientationAndOrdering(line_parameters, arc_parameters)
% This function takes the parameter inputs and produces parameter sets such
% that the line is first, it is oriented so that it ends at the junction
% with the arc, and the arc starts at the junction. It also forces the
% station of the line to start at 0 and the base point of the line to be
% the start point, if they do not already.

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_unit_tangent_vector = line_parameters(1,1:2);
line_base_point_xy       = line_parameters(1,3:4);
line_s_start             = line_parameters(1,5);
line_s_end               = line_parameters(1,6);
line_start_xy            = line_base_point_xy + line_unit_tangent_vector*line_s_start;
line_end_xy              = line_base_point_xy + line_unit_tangent_vector*line_s_end;

% Make sure the line segment is well-formed, e.g. the station at the end is
% larger than the station at the start. If not, need to correct
if line_s_end<line_s_start
    % % Flip the order
    % line_s_start             = line_parameters(1,6);
    % line_s_end               = line_parameters(1,5);

    % Flip the vector
    line_unit_tangent_vector = -line_unit_tangent_vector;

end

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy                = arc_parameters(1,1:2);
arc_radius                   = arc_parameters(1,3);
arc_start_angle_in_radians   = arc_parameters(1,4);
arc_end_angle_in_radians     = arc_parameters(1,5);
arc_start_xy                 = arc_center_xy + arc_radius*[cos(arc_start_angle_in_radians) sin(arc_start_angle_in_radians)];
arc_end_xy                   = arc_center_xy + arc_radius*[cos(arc_end_angle_in_radians) sin(arc_end_angle_in_radians)];
arc_is_counter_clockwise     = arc_parameters(1,7);

% Find the change in angle of the arc
arc_start_unit_vector        = [cos(arc_start_angle_in_radians) sin(arc_start_angle_in_radians)];
arc_end_unit_vector          = [cos(arc_end_angle_in_radians)   sin(arc_end_angle_in_radians)  ];
if arc_is_counter_clockwise
    cross_product_direction = 1;
else
    cross_product_direction = -1;
end
change_in_arc_angle = fcn_geometry_findAngleUsing2PointsOnCircle([0 0],1, arc_start_unit_vector, arc_end_unit_vector, cross_product_direction);


% Find line and arc's join points, e.g. where they meet. This can
% happen at either end
distances_to_check = sum((...
    [line_start_xy; line_start_xy; line_end_xy; line_end_xy] - [arc_start_xy; arc_end_xy; arc_start_xy; arc_end_xy]).^2,2).^0.5;
[~,closest_pair] = min(distances_to_check);

% Fix the line
switch closest_pair
    case 1 % Line start, arc start
        % Line segment is pointing away from junction, need to fix orientation,
        % base point, and station
        corrected_line_unit_tangent_vector = -line_unit_tangent_vector;
        corrected_line_base_point_xy       = line_end_xy;

        % The arc is leaving the junction at its start. This is
        % correct so just pass through the variables
        corrected_arc_is_counter_clockwise   = arc_is_counter_clockwise;
        corrected_arc_start_angle_in_radians = atan2(arc_start_unit_vector(2),arc_start_unit_vector(1));
        corrected_arc_end_angle_in_radians   = corrected_arc_start_angle_in_radians + change_in_arc_angle;

    case 2 % Line start, arc end
        % Line segment is pointing away from junction, need to fix orientation,
        % base point, and station
        corrected_line_unit_tangent_vector = -line_unit_tangent_vector;
        corrected_line_base_point_xy       = line_end_xy;

        % The arc is entering the junction at its end. This is
        % not correct. Need to "flip" the arc's orientation.
        if 1==arc_is_counter_clockwise
            corrected_arc_is_counter_clockwise = 0;
        else
            corrected_arc_is_counter_clockwise = 1;
        end
        corrected_arc_start_angle_in_radians = atan2(arc_end_unit_vector(2),arc_end_unit_vector(1));
        corrected_arc_end_angle_in_radians   = corrected_arc_start_angle_in_radians - change_in_arc_angle;

    case 3 % Line end, arc start
        % Line segment is pointing into junction, need to just fix base point and station
        corrected_line_unit_tangent_vector = line_unit_tangent_vector;
        corrected_line_base_point_xy       = line_start_xy;

        % The arc is leaving the junction at its start. This is
        % correct so just pass through the variables
        corrected_arc_is_counter_clockwise = arc_is_counter_clockwise;
        corrected_arc_start_angle_in_radians = atan2(arc_start_unit_vector(2),arc_start_unit_vector(1));
        corrected_arc_end_angle_in_radians   = corrected_arc_start_angle_in_radians + change_in_arc_angle;

    case 4 % Line end, arc end
        % Line segment is pointing into junction, need to just fix base point and station
        corrected_line_unit_tangent_vector = line_unit_tangent_vector;
        corrected_line_base_point_xy       = line_start_xy;

        % The arc is entering the junction at its end. This is
        % not correct. Need to "flip" the arc's orientation.
        if 1==arc_is_counter_clockwise
            corrected_arc_is_counter_clockwise = 0;
        else
            corrected_arc_is_counter_clockwise = 1;
        end
        corrected_arc_start_angle_in_radians = atan2(arc_end_unit_vector(2),arc_end_unit_vector(1));
        corrected_arc_end_angle_in_radians   = corrected_arc_start_angle_in_radians - change_in_arc_angle;

    otherwise
        error('Impossible case encountered - must stop!');
end

% Set the line start and end by distances only, starting at 0 and ending at
% total distance
corrected_line_s_start             = 0;
corrected_line_s_end               = sum((line_end_xy - line_start_xy).^2,2).^0.5;

clean_line_parameters(1,1:2)  = corrected_line_unit_tangent_vector;
clean_line_parameters(1,3:4)  = corrected_line_base_point_xy;
clean_line_parameters(1,5:6)  = [corrected_line_s_start corrected_line_s_end];

clean_arc_parameters(1,1:2)   = arc_parameters(1,1:2); % center of the arc does not change
clean_arc_parameters(1,3)     = arc_parameters(1,3);   % radius of the arc does not change
clean_arc_parameters(1,4:5)   = [corrected_arc_start_angle_in_radians corrected_arc_end_angle_in_radians];
clean_arc_parameters(1,6)     = arc_parameters(1,6);   % flag is circle
clean_arc_parameters(1,7)     = corrected_arc_is_counter_clockwise;


end % ends fcn_INTERNAL_fixOrientationAndOrdering


%% fcn_INTERNAL_findLineArcIntersection
function  intersection_point = fcn_INTERNAL_findLineArcIntersection(clean_line_parameters,clean_arc_parameters)

% Calculate needed values from parameter sets

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_unit_tangent_vector     = clean_line_parameters(1,1:2);
line_base_point_xy           = clean_line_parameters(1,3:4);
% line_s_start               = clean_line_parameters(1,5);
% line_s_end                 = clean_line_parameters(1,6);
% line_start_xy              = line_base_point_xy + line_unit_tangent_vector*line_s_start;
% line_end_xy                = line_base_point_xy + line_unit_tangent_vector*line_s_end;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy                = clean_arc_parameters(1,1:2);
arc_radius                   = clean_arc_parameters(1,3);
arc_start_angle_in_radians   = clean_arc_parameters(1,4);
arc_end_angle_in_radians     = clean_arc_parameters(1,5);
arc_start_xy                 = arc_center_xy + arc_radius*[cos(arc_start_angle_in_radians) sin(arc_start_angle_in_radians)];
arc_end_xy                   = arc_center_xy + arc_radius*[cos(arc_end_angle_in_radians) sin(arc_end_angle_in_radians)];
% arc_is_circle                = clean_arc_parameters(1,6);
arc_is_counter_clockwise     = clean_arc_parameters(1,7);
% change_in_arc_angle = arc_end_angle_in_radians-arc_start_angle_in_radians; % Find the change in angle of the arc


% With the cleaned parameters, the line vector always
% points toward the joint of the line and the arc.
% to_joint_line_unit_tangent_vector = line_unit_tangent_vector;
% to_joint_line_unit_ortho_vector= to_joint_line_unit_tangent_vector*[0 1; -1 0];


% flag_intersection_points_found = 0;
% Check if the line and circle intersect
if line_unit_tangent_vector(1)==0
    slope = inf;
    intercept = line_base_point_xy(1);
else
    slope = line_unit_tangent_vector(2)/line_unit_tangent_vector(1);
    intercept = line_base_point_xy(2) - slope*line_base_point_xy(1); % b = y - m*x
end

% Use MATLAB's linecirc algorithm to find intersections
[xout,yout] = linecirc(slope,intercept,arc_center_xy(1,1),arc_center_xy(1,2),arc_radius);

% Check results of above
if 1==0
    figure(233);
    clf;
    hold on;
    grid on;
    axis equal
    angles = (0:1*pi/180:2*pi)';
    test_circle_XY = ones(length(angles),1)*arc_center_xy + arc_radius*[cos(angles) sin(angles)];
    plot(test_circle_XY(:,1),test_circle_XY(:,2),'r-','LineWidth',3);

    if ~isinf(slope)
        line_xdata = (min(test_circle_XY(:,1)):(arc_radius/100):max(test_circle_XY(:,1)))';
        line_ydata = slope*line_xdata + intercept*ones(length(line_xdata),1);
    else
        line_ydata = (min(test_circle_XY(:,2)):(arc_radius/100):max(test_circle_XY(:,2)))';
        line_xdata = ones(length(line_ydata),1)*intercept;
    end
    plot(line_xdata,line_ydata,'b-','LineWidth',3);

    intersections = [xout', yout'];
    plot(intersections(:,1),intersections(:,2),'g.','MarkerSize',20);
end

if ~isnan(xout)
    % intersection points were found!

    % Which point(s) to keep?
    two_intersection_points = [xout', yout'];

    % Are the intersections within the arc range that we were given? To
    % check this, we use the three points on the arc - the start, the
    % intersection, and the end to calculate the arc direction. We then
    % check to see if it is the same as the given direction - if it is, the
    % point is on the arc.
    potential_intersection_points = nan(size(two_intersection_points));
    for ith_row = 1:length(two_intersection_points(:,1))
        intersection_is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(arc_start_xy, two_intersection_points(ith_row,:), arc_end_xy,-1);
        if arc_is_counter_clockwise == intersection_is_counterClockwise
            potential_intersection_points(ith_row,:) = two_intersection_points(ith_row,:);
        else
            potential_intersection_points(ith_row,:) = [nan nan];
        end
    end

    if isequal(potential_intersection_points(1,:),potential_intersection_points(2,:))
        intersection_point = potential_intersection_points(1,:);
    else
        % Find which point is closest to the line's start point
        vectors_from_line_base_point = potential_intersection_points - [1;1]*line_base_point_xy;

        % Distances are dot product with the line's vector
        distances = sum([1; 1]*line_unit_tangent_vector.*vectors_from_line_base_point,2).^0.5;
        if distances(1)<distances(2)
            intersection_point = potential_intersection_points(1,:);
        else
            intersection_point = potential_intersection_points(2,:);
        end
    end

else
    intersection_point = [nan nan];
end
end % Ends fcn_INTERNAL_findLineArcIntersection



%% fcn_INTERNAL_findShiftToMatchArcToLine
function [delta_transverse, delta_station, desired_closest_arc_point_to_joint, desired_closest_line_point_to_joint,spiral_join_parameters] = fcn_INTERNAL_findShiftToMatchArcToLine(clean_line_parameters, clean_arc_parameters,continuity_level, intersection_point)
% Calculates the delta amount to match the arc to the line. The delta
% values are measured FROM desired point TO actual point

spiral_join_parameters = []; % Initialize output. This is ONLY used for spiral fits.

% Calculate needed values from parameter sets

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_unit_tangent_vector     = clean_line_parameters(1,1:2);
line_base_point_xy           = clean_line_parameters(1,3:4);
% line_s_start               = clean_line_parameters(1,5);
line_s_end                 = clean_line_parameters(1,6);
% line_start_xy              = line_base_point_xy + line_unit_tangent_vector*line_s_start;
line_end_xy                = line_base_point_xy + line_unit_tangent_vector*line_s_end;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy                = clean_arc_parameters(1,1:2);
arc_radius                   = clean_arc_parameters(1,3);
arc_start_angle_in_radians   = clean_arc_parameters(1,4);
arc_end_angle_in_radians     = clean_arc_parameters(1,5);
arc_start_xy                 = arc_center_xy + arc_radius*[cos(arc_start_angle_in_radians) sin(arc_start_angle_in_radians)];
arc_end_xy                   = arc_center_xy + arc_radius*[cos(arc_end_angle_in_radians) sin(arc_end_angle_in_radians)];
% arc_is_circle                = clean_arc_parameters(1,6);
arc_is_counter_clockwise     = clean_arc_parameters(1,7);
% change_in_arc_angle = arc_end_angle_in_radians-arc_start_angle_in_radians; % Find the change in angle of the arc


% With the cleaned parameters, the line vector always
% points toward the joint of the line and the arc.
to_joint_unit_ortho_vector   = line_unit_tangent_vector*[0 1; -1 0];


if 0 == continuity_level
    % For C0 continuity of a line to an arc, the closest point of the
    % desired joint to the arc is either the intersection of the line and
    % the arc, or where the arc ends
    if any(isnan(intersection_point))
        desired_closest_arc_point_to_joint = line_end_xy;
        desired_closest_line_point_to_joint = line_end_xy;
    else
        desired_closest_arc_point_to_joint = intersection_point;
        desired_closest_line_point_to_joint = intersection_point;
    end

elseif 1 == continuity_level
    % For C1 continuity of a line to an arc, the closest point of the
    % desired joint to the arc is always going to be the point on the arc
    % where the circle is tangent to the line. This tangent point can be
    % found by using unit vectors that are aligned with the line segment,
    % and orthogonal to the segment.

    % Calculate the offset from the arc's circle to the line
    vector_from_line_anti_joint_to_arc_center = arc_center_xy - line_base_point_xy;
    % unit_vector_from_line_anti_joint_to_arc_center = fcn_geometry_calcUnitVector(vector_from_line_anti_joint_to_arc_center);
    signed_distance_along_line_to_joint = sum(line_unit_tangent_vector.*vector_from_line_anti_joint_to_arc_center,2);

    % signed_distance_ortho_line_to_arc_center = sum(to_joint_unit_ortho_vector.*vector_from_line_anti_joint_to_arc_center,2);
    % crossProduct_to_find_arc_sign = cross([to_joint_unit_tangent_vector 0],[unit_vector_from_line_anti_joint_to_arc_center 0]);
    % arc_direction_relative_to_line_to_joint = crossProduct_to_find_arc_sign(3);

    % If everything is done right, signed distance along the line to circle
    % is always positive
    assert(signed_distance_along_line_to_joint>=0);

    desired_closest_arc_point_to_joint = line_end_xy;
    desired_closest_line_point_to_joint = line_end_xy;

elseif 2 == continuity_level
    % For C2 continuity of a line to an arc, the arc must connect to the
    % line via a spiral. Calculation of the spiral is difficult and
    % requires numerical iteration, and the inputs require the calculation
    % of the offset of the outer part of the circle from the line. This
    % offset in the direction of the circle center MUST be strictly
    % positive or the spiral cannot work.

    % Calculate the offset from the arc's circle to the line
    vector_from_line_anti_joint_to_arc_center = arc_center_xy - line_base_point_xy;

    % unit_vector_from_line_anti_joint_to_arc_center = fcn_geometry_calcUnitVector(vector_from_line_anti_joint_to_arc_center);
    % signed_distance_along_line_to_joint = sum(line_unit_tangent_vector.*vector_from_line_anti_joint_to_arc_center,2);

    if 1==arc_is_counter_clockwise
        vector_from_joint_to_arc_center = to_joint_unit_ortho_vector;
    else
        vector_from_joint_to_arc_center = -to_joint_unit_ortho_vector;
    end
    unit_vector_from_joint_to_arc_center = fcn_geometry_calcUnitVector(vector_from_joint_to_arc_center);

    signed_distance_ortho_line_to_arc_center = sum(unit_vector_from_joint_to_arc_center.*vector_from_line_anti_joint_to_arc_center,2);
    offset_from_arc_edge_to_line = signed_distance_ortho_line_to_arc_center - arc_radius;

    if offset_from_arc_edge_to_line>0
        h0 = 0; % Initial heading
        x0 = line_base_point_xy(1,1);
        y0 = line_base_point_xy(1,2);
        K0 = 0; % Initial curvature
        Kf = 1/arc_radius;
        spiralLength = fcn_INTERNAL_findLengthFromOffset(offset_from_arc_edge_to_line, h0, x0, y0, K0, Kf);
        
        % Find the angle that the spiral ends at - it is always the spiral length/2
        [x_spiralEnd,y_spiralEnd] = fcn_geometry_extractXYfromSTSpiral(spiralLength,spiralLength,h0,x0,y0,K0,Kf,-1);

        spiral_end_angles_radians = spiralLength/2;

        % Find where the circle center is at based on this end angle
        unit_tangent_vector = [cos(spiral_end_angles_radians) sin(spiral_end_angles_radians)];
        unit_orthogonal_vector = unit_tangent_vector*[0 1; -1 0];
        spiral_predicted_circle_center = arc_radius*unit_orthogonal_vector + [x_spiralEnd(end) y_spiralEnd(end)];
        offset_error = arc_center_xy - spiral_predicted_circle_center;
        x0 = offset_error(1);
        [x_spiralEnd,y_spiralEnd] = fcn_geometry_extractXYfromSTSpiral(spiralLength,spiralLength,h0,x0,y0,K0,Kf,-1);

        if 1==0
            % Set up station coordinates            
            s  = (0:0.01:1)'*spiralLength;

            % Call the function fcn_geometry_extractXYfromSTSpiral to predict the
            % spiral and calculate the offsets, plotting the results
            fcn_geometry_extractXYfromSTSpiral(s,spiralLength,h0,x0,y0,K0,Kf,(1234));

            fcn_geometry_plotGeometry('arc',clean_arc_parameters);

            
            % Find the center of the circle tangent at the end of the spiral
            % Find the unit vector (need to do this analytically!)            
            s_tangent = [0.999999 1]'*spiralLength;
            [x_tangent,y_tangent] = fcn_geometry_extractXYfromSTSpiral(s_tangent,spiralLength,h0,x0,y0,K0,Kf);
            unit_tangent = fcn_geometry_calcUnitVector([diff(x_tangent) diff(y_tangent)]);
            unit_orthogonal = unit_tangent*[0 1; -1 0];
            calculated_circle_center = arc_radius*unit_orthogonal + [x_tangent(end) y_tangent(end)];

            % Plot the circle's center
            plot(calculated_circle_center(:,1),calculated_circle_center(:,2),'r+');

            % Plot the circle
            fcn_geometry_plotCircle(arc_center_xy, arc_radius,'r-',(1234));
        end % Ends plotting

        % Make sure x0 and final angles are within the line segment and the
        % arc respectively, 
        if 0>x0 || x0>line_s_end
            error('Spiral fit does not fit within x-range of line. Unable to continue.')
        end

        spiral_join_xy = [x_spiralEnd,y_spiralEnd];
        intersection_is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(arc_start_xy, spiral_join_xy, arc_end_xy,-1);
        if arc_is_counter_clockwise ~= intersection_is_counterClockwise
            error('Spiral fit does not fit within angle range of arc. Unable to continue.')
        end

        delta_transverse = 0;
        delta_station    = 0;
        desired_closest_arc_point_to_joint = spiral_join_xy;
        desired_closest_line_point_to_joint = [x0 0];
        spiral_join_parameters = [spiralLength,h0,x0,y0,K0,Kf];
        return;
    else
        error('Spirals cannot be formed between arc curves and intersecting lines. This is geometrically impossible');
    end % Ends if statement to check if spiral is possible

    desired_closest_arc_point_to_joint = line_end_xy;
    desired_closest_line_point_to_joint = line_end_xy;

else
    error('This continuity not possible yet')
end


% Depending on which side of the line the arc is located, either need to
% subtract or add on the circle radius. Note: the result is designed to be
% positive if the line doesn't quite meet the circle as a tangent creating
% a gap - in this case, a spiral is needed. The result is negative if the
% line would intersect the circle; in this case, a line shift offset is
% needed.

if 0 == continuity_level
    difference_vector_from_line_end_to_arc_start = arc_start_xy - desired_closest_line_point_to_joint;
    delta_transverse = sum(to_joint_unit_ortho_vector.*difference_vector_from_line_end_to_arc_start,2);
    delta_station    = sum(line_unit_tangent_vector.*difference_vector_from_line_end_to_arc_start,2);
elseif 1 == continuity_level
    St_vector_from_line_end_to_center = to_joint_unit_ortho_vector.*vector_from_line_anti_joint_to_arc_center;
    unit_St_vector_from_line_end_to_center = fcn_geometry_calcUnitVector(St_vector_from_line_end_to_center);

    desired_arc_center_xy = desired_closest_line_point_to_joint + arc_radius*unit_St_vector_from_line_end_to_center;
    difference_vector_from_desired_to_actual_circle_center = arc_center_xy - desired_arc_center_xy;
    delta_transverse = sum(to_joint_unit_ortho_vector.*difference_vector_from_desired_to_actual_circle_center,2);
    delta_station    = sum(line_unit_tangent_vector.*difference_vector_from_desired_to_actual_circle_center,2);
    % if signed_distance_ortho_line_to_arc_center < 0
    %     delta_transverse   = -1*delta_transverse;
    % end

elseif 2 == continuity_level
    error('Not coded yet');
else
    error('This continuity not possible yet')
end
end % Ends fcn_INTERNAL_findShiftToMatchArcToLine

%% fcn_INTERNAL_performShift
function [revised_line_parameters,revised_arc_parameters] = fcn_INTERNAL_performShift(flag_arc_is_first, clean_line_parameters, clean_arc_parameters,~, delta_transverse,delta_station, arc_start_point, line_end_point)

% Calculate needed values from parameter sets

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_unit_tangent_vector     = clean_line_parameters(1,1:2);
line_base_point_xy           = clean_line_parameters(1,3:4);
% line_s_start               = clean_line_parameters(1,5);
% line_s_end                 = clean_line_parameters(1,6);
% line_start_xy              = line_base_point_xy + line_unit_tangent_vector*line_s_start;
% line_end_xy                = line_base_point_xy + line_unit_tangent_vector*line_s_end;
assert(isequal(round(line_unit_tangent_vector,4),[1 0]));
assert(isequal(round(line_base_point_xy,4),[0 0]));

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy                = clean_arc_parameters(1,1:2);
arc_radius                   = clean_arc_parameters(1,3);
% arc_start_angle_in_radians   = clean_arc_parameters(1,4);
arc_end_angle_in_radians     = clean_arc_parameters(1,5);
% arc_start_xy                 = arc_center_xy + arc_radius*[cos(arc_start_angle_in_radians) sin(arc_start_angle_in_radians)];
arc_end_xy                     = arc_center_xy + arc_radius*[cos(arc_end_angle_in_radians) sin(arc_end_angle_in_radians)];
arc_is_circle                = clean_arc_parameters(1,6);
% arc_is_counter_clockwise     = clean_arc_parameters(1,7);
% change_in_arc_angle = arc_end_angle_in_radians-arc_start_angle_in_radians; % Find the change in angle of the arc


% % Need to find where the line and arc join
% join_point_xy = line_base_point_xy + signed_distance_along_line_to_join*to_joint_line_unit_tangent_vector;
%
%
% % Calculate how much shift is needed to connect the line exactly to the
% % arc
% vector_from_arc_center_to_join = join_point_xy - arc_center_xy;
% unit_vector_from_arc_center_to_join = fcn_geometry_calcUnitVector(vector_from_arc_center_to_join);
% point_shift_xy = -1*offset_dist_from_line_toward_circle*unit_vector_from_arc_center_to_join;
%

% Calculate the offset from the arc's circle to the line
vector_from_line_anti_joint_to_arc_center = arc_center_xy - line_base_point_xy;
unit_vector_from_line_anti_joint_to_arc_center = fcn_geometry_calcUnitVector(vector_from_line_anti_joint_to_arc_center);

crossProduct_to_find_arc_sign = cross([line_unit_tangent_vector 0],[unit_vector_from_line_anti_joint_to_arc_center 0]);
arc_direction_relative_to_line_to_joint = crossProduct_to_find_arc_sign(3);


point_shift_xy = [-delta_station -delta_transverse];


%%% Fix the line
if 0==flag_arc_is_first
    % The line is first, do not shift it
    new_line_unit_tangent_vector = line_unit_tangent_vector;
    new_line_base_point_xy       = line_base_point_xy;
    new_line_s_start             = 0;
    new_line_s_end               = sum((line_base_point_xy - line_end_point).^2,2).^0.5;
    if arc_direction_relative_to_line_to_joint>=0
        % Arc veers to the left relative to the line's vector
        new_arc_is_counter_clockwise = 1;
    else
        % Arc veers to the right relative to the line's vector
        new_arc_is_counter_clockwise = 0;
    end

    new_arc_center_xy = arc_center_xy + point_shift_xy;
    new_arc_end_xy    = arc_end_xy    + point_shift_xy;

    vector_from_arc_center_to_join = arc_start_point - new_arc_center_xy;
    arc_angle_at_join = atan2(vector_from_arc_center_to_join(1,2),vector_from_arc_center_to_join(1,1));

    % The angle the arc ends is where the arc stops
    vector_from_arc_center_to_end = new_arc_end_xy - new_arc_center_xy;
    arc_angle_at_end = atan2(vector_from_arc_center_to_end(1,2),vector_from_arc_center_to_end(1,1));

    new_arc_start_angle_in_radians = arc_angle_at_join;
    new_arc_end_angle_in_radians   = arc_angle_at_end;


else
    % The line is last, need to shift it
    new_line_unit_tangent_vector = -1*line_unit_tangent_vector;
    new_line_base_point_xy       = line_end_point - point_shift_xy;
    new_line_s_start             = 0;
    new_line_s_end               = sum((line_base_point_xy - line_end_point).^2,2).^0.5;
    if arc_direction_relative_to_line_to_joint>=0
        new_arc_is_counter_clockwise = 0;
    else
        new_arc_is_counter_clockwise = 1;
    end

    % The line is last, do not shift the arc
    new_arc_center_xy = arc_center_xy;

    % The angle the arc starts at the join point
    vector_from_arc_center_to_join = arc_start_point - new_arc_center_xy;
    arc_angle_at_join = atan2(vector_from_arc_center_to_join(1,2),vector_from_arc_center_to_join(1,1));

    % The angle the arc ends is where the arc stops
    new_arc_end_xy = arc_end_xy;
    vector_from_arc_center_to_end = new_arc_end_xy - new_arc_center_xy;
    arc_angle_at_end = atan2(vector_from_arc_center_to_end(1,2),vector_from_arc_center_to_end(1,1));

    new_arc_end_angle_in_radians   = arc_angle_at_join;
    new_arc_start_angle_in_radians = arc_angle_at_end;
end
revised_line_parameters(1,1:2) = new_line_unit_tangent_vector;
revised_line_parameters(1,3:4) = new_line_base_point_xy;
revised_line_parameters(1,5)   = new_line_s_start;
revised_line_parameters(1,6)   = new_line_s_end;


revised_arc_parameters(1,1:2) = new_arc_center_xy;
revised_arc_parameters(1,3)   = arc_radius;
revised_arc_parameters(1,4)   = new_arc_start_angle_in_radians;
revised_arc_parameters(1,5)   = new_arc_end_angle_in_radians;
revised_arc_parameters(1,6)   = arc_is_circle;
revised_arc_parameters(1,7)   = new_arc_is_counter_clockwise;

% Fix the arc angles to be between 0 and 2*pi
revised_arc_parameters(1,4:5) = mod(revised_arc_parameters(1,4:5),2*pi);
end % Ends fcn_INTERNAL_performShift

%% fcn_INTERNAL_convertParametersToStOrientation
function [st_line_parameters, st_arc_parameters, St_transform, rotation_angle] = fcn_INTERNAL_convertParametersToStOrientation(clean_line_parameters, clean_arc_parameters)
% Calculate St-equivalent values from parameter sets

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_unit_tangent_vector     = clean_line_parameters(1,1:2);
line_base_point_xy           = clean_line_parameters(1,3:4);
line_s_start                 = clean_line_parameters(1,5);
line_s_end                   = clean_line_parameters(1,6);

line_angle = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));
st_line_parameters(1,1:2)    = [1 0];
st_line_parameters(1,3:4)    = [0 0];
st_line_parameters(1,5)      = line_s_start;
st_line_parameters(1,6)      = line_s_end;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy                = clean_arc_parameters(1,1:2);
arc_radius                   = clean_arc_parameters(1,3);
arc_start_angle_in_radians   = clean_arc_parameters(1,4);
arc_end_angle_in_radians     = clean_arc_parameters(1,5);
arc_is_circle                = clean_arc_parameters(1,6);
arc_is_counter_clockwise     = clean_arc_parameters(1,7);

% Use the SE2 toolbox to transform. Unfortunately, we have to force a
% translation FIRST, then rotation. (There's no obvious way to do this in
% one step).
translation_to_St       = -line_base_point_xy;
rotation_angle          = line_angle;
transformMatrix_translation_into_St = se2(0,'theta',translation_to_St);
transformMatrix_rotation_into_St    = se2(rotation_angle,'theta',[0 0]);
transformMatrix_into_St             = transformMatrix_rotation_into_St*transformMatrix_translation_into_St;
arc_center_St                 = transform(transformMatrix_into_St,arc_center_xy);

arc_start_angle_in_radians_St = arc_start_angle_in_radians+rotation_angle;
arc_end_angle_in_radians_St   = arc_end_angle_in_radians  +rotation_angle;

st_arc_parameters(1,1:2)      = arc_center_St;
st_arc_parameters(1,3)        = arc_radius;
st_arc_parameters(1,4)        = arc_start_angle_in_radians_St;
st_arc_parameters(1,5)        = arc_end_angle_in_radians_St;
st_arc_parameters(1,6)        = arc_is_circle;
st_arc_parameters(1,7)        = arc_is_counter_clockwise;

St_transform = transformMatrix_into_St;
end % fcn_INTERNAL_convertParametersToStOrientation



%% fcn_INTERNAL_convertParametersOutOfStOrientation
function [revised_line_parameters, revised_arc_parameters, revised_spiral_join_parameters] = ...
    fcn_INTERNAL_convertParametersOutOfStOrientation(...
    revised_line_parameters_St, revised_arc_parameters_St, St_transform, rotation_angle, revised_spiral_join_parameters_St)
% Calculate XY parameters from St-equivalent values of parameter sets

transformMatrix_outOf_St = inv(St_transform);


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_unit_tangent_vector     = revised_line_parameters_St(1,1:2);
line_base_point_xy           = revised_line_parameters_St(1,3:4);
line_s_start                 = revised_line_parameters_St(1,5);
line_s_end                   = revised_line_parameters_St(1,6);

revised_vector_head               = transform(transformMatrix_outOf_St,line_unit_tangent_vector);
revised_vector_tail               = transform(transformMatrix_outOf_St,[0 0]);
revised_unit_tangent_vector       = revised_vector_head - revised_vector_tail;
revised_line_base_point_xy        = transform(transformMatrix_outOf_St,line_base_point_xy);

revised_line_parameters(1,1:2)    = revised_unit_tangent_vector;
revised_line_parameters(1,3:4)    = revised_line_base_point_xy;
revised_line_parameters(1,5)      = line_s_start;
revised_line_parameters(1,6)      = line_s_end;


% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy                = revised_arc_parameters_St(1,1:2);
arc_radius                   = revised_arc_parameters_St(1,3);
arc_start_angle_in_radians   = revised_arc_parameters_St(1,4);
arc_end_angle_in_radians     = revised_arc_parameters_St(1,5);
arc_is_circle                = revised_arc_parameters_St(1,6);
arc_is_counter_clockwise     = revised_arc_parameters_St(1,7);



arc_center_St                 = transform(transformMatrix_outOf_St,arc_center_xy);
arc_start_angle_in_radians_St = arc_start_angle_in_radians - rotation_angle;
arc_end_angle_in_radians_St   = arc_end_angle_in_radians   - rotation_angle;

revised_arc_parameters(1,1:2)      = arc_center_St;
revised_arc_parameters(1,3)        = arc_radius;
revised_arc_parameters(1,4)        = arc_start_angle_in_radians_St;
revised_arc_parameters(1,5)        = arc_end_angle_in_radians_St;
revised_arc_parameters(1,6)        = arc_is_circle;
revised_arc_parameters(1,7)        = arc_is_counter_clockwise;

% Get the spiral parameters
if ~isempty(revised_spiral_join_parameters_St)
    spiralLength      = revised_spiral_join_parameters_St(1,1);
    h0                = revised_spiral_join_parameters_St(1,2);
    x0                = revised_spiral_join_parameters_St(1,3);
    y0                = revised_spiral_join_parameters_St(1,4);
    K0                = revised_spiral_join_parameters_St(1,5);
    Kf                = revised_spiral_join_parameters_St(1,6);

    revised_spiral_join_parameters(1,1) = spiralLength;
    revised_spiral_join_parameters(1,2) = h0 - rotation_angle;
    new_spiral_XY                       = transform(transformMatrix_outOf_St,[x0 y0]);
    revised_spiral_join_parameters(1,3) = new_spiral_XY(1,1);
    revised_spiral_join_parameters(1,4) = new_spiral_XY(1,2);
    revised_spiral_join_parameters(1,5) = K0;
    revised_spiral_join_parameters(1,6) = Kf;

else
    revised_spiral_join_parameters = [];
end

end % fcn_INTERNAL_convertParametersOutOfStOrientation



%% fcn_INTERNAL_findLengthFromOffset
function sprialLength = fcn_INTERNAL_findLengthFromOffset(offset, h0, x0, y0, K0, Kf)
function_to_optimize = @(x)fcn_INTERNAL_calcUnitSpiralOffsetError(x,offset, h0, x0, y0, K0, Kf);
% options = optimset('Display','iter','PlotFcns',@optimplotfval);
% sprialLength = fminsearch(function_to_optimize,1,options);
sprialLength = fminsearch(function_to_optimize,1);

end % Ends fcn_INTERNAL_findLengthFromOffset

%% fcn_INTERNAL_findLengthFromOffset
function error = fcn_INTERNAL_calcUnitSpiralOffsetError(length_spiral, offset, h0, x0, y0, K0, Kf)

% PERFORM UNIT CIRCLE ANALYSIS
% % Initial heading, position, and curvature of the spiral
% h0 = 0;
% x0 = 0;
% y0 = 0;
% K0 = 0;
% Kf = 1; % Final curvature forced to have a radius of 1
radius = (1/Kf);

% Call the function
[x_arc,y_arc] = fcn_geometry_extractXYfromSTSpiral(length_spiral,length_spiral,h0,x0,y0,K0,Kf,-1);

% Find the angle that the spiral ends at - it is always the spiral length/2
unit_spiral_end_angles_radians = length_spiral/2;

% Find where the circle center is at
unit_tangent_vector = [cos(unit_spiral_end_angles_radians) sin(unit_spiral_end_angles_radians)];
unit_orthogonal_vector = unit_tangent_vector*[0 1; -1 0];
circle_center = radius*unit_orthogonal_vector + [x_arc(end) y_arc(end)];

actual_y_offset            = circle_center(1,2) - radius;
error = abs(offset-actual_y_offset);


end % Ends fcn_INTERNAL_findLengthFromOffset