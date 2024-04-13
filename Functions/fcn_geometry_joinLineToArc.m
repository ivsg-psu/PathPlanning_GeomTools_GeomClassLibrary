function [revised_line_parameters, revised_arc_parameters] = fcn_geometry_joinLineToArc(line_parameters, arc_parameters, flag_arc_is_first, varargin)
%% fcn_geometry_joinLineToArc
% Revises the geometric parameters of an arc and line such that they
% join. It does this by checking the offset between the two objects at the
% join location. If the offset is less than a threshold (default is 0.1
% meter), the second geometric object is shifted to force alignment.
%
% Format:
% [revised_line_parameters, revised_arc_parameters]  = fcn_geometry_joinLineToArc(line_parameters, arc_parameters, flag_arc_is_first, (threshold), (fig_num))
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
%      that this offset is removed by shifting.
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
% DEPENDENCIES:
%      
%      (none)
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_joinLineToArc
% for a full test suite.
%
% This function was written on 2024_04_12 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_04_12 - S. Brennan
% -- wrote the code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==5 && isequal(varargin{end},-1))
    flag_do_debug = 0; % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % Flag to plot the results for debugging
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
threshold = 0.1;
if (4<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        threshold = temp;
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

if flag_do_debug
    figure(debug_fig_num);
    clf;

    % Plot the inputs
    subplot(2,1,1);
    hold on;
    grid on;
    axis equal;
    xlabel('X [meters]');
    ylabel('Y [meters]')

    fcn_geometry_plotGeometry('line',line_parameters);
    fcn_geometry_plotGeometry('arc',arc_parameters);
    title('Inputs');
end

% Fix the parameters to make the line segment first, arc second, and make
% sure the line and arc point into and then out of the junction
% respectively.
[clean_line_parameters, clean_arc_parameters] = fcn_INTERNAL_fixOrientationAndOrdering(line_parameters, arc_parameters);

if flag_do_debug
    figure(debug_fig_num);

    % Plot the cleaned inputs
    subplot(2,1,2);
    hold on;
    grid on;
    axis equal;
    xlabel('X [meters]');
    ylabel('Y [meters]')

    fcn_geometry_plotGeometry('line',clean_line_parameters);
    fcn_geometry_plotGeometry('arc',clean_arc_parameters);
    title('Cleaned inputs');
end

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_unit_tangent_vector = clean_line_parameters(1,1:2);
line_base_point_xy       = clean_line_parameters(1,3:4);
line_s_start             = clean_line_parameters(1,5);
line_s_end               = clean_line_parameters(1,6);
line_start_xy            = line_base_point_xy + line_unit_tangent_vector*line_s_start;
line_end_xy              = line_base_point_xy + line_unit_tangent_vector*line_s_end;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy                = clean_arc_parameters(1,1:2);
arc_radius                   = clean_arc_parameters(1,3);
arc_start_angle_in_radians   = clean_arc_parameters(1,4);
arc_end_angle_in_radians     = clean_arc_parameters(1,5);
arc_start_xy                 = arc_center_xy + arc_radius*[cos(arc_start_angle_in_radians) sin(arc_start_angle_in_radians)];
arc_end_xy                   = arc_center_xy + arc_radius*[cos(arc_end_angle_in_radians) sin(arc_end_angle_in_radians)];
arc_is_circle                = clean_arc_parameters(1,6);
arc_is_counter_clockwise     = clean_arc_parameters(1,7);
change_in_arc_angle = arc_end_angle_in_radians-arc_start_angle_in_radians; % Find the change in angle of the arc


% With the cleaned parameters, the line vector always
% points toward the joint of the line and the arc.
to_joint_line_unit_tangent_vector = line_unit_tangent_vector;
to_joint_line_unit_ortho_vector= to_joint_line_unit_tangent_vector*[0 1; -1 0];

% Calculate the offset from the arc's circle to the line
vector_from_line_anti_join_to_arc_center = arc_center_xy - line_base_point_xy;
unit_vector_from_line_anti_join_to_arc_center = fcn_geometry_calcUnitVector(vector_from_line_anti_join_to_arc_center);
signed_distance_along_line_to_join = sum(to_joint_line_unit_tangent_vector.*vector_from_line_anti_join_to_arc_center,2);
signed_distance_ortho_line_to_arc_center = sum(to_joint_line_unit_ortho_vector.*vector_from_line_anti_join_to_arc_center,2);
crossProduct_to_find_arc_sign = cross([to_joint_line_unit_tangent_vector 0],[unit_vector_from_line_anti_join_to_arc_center 0]);
arc_direction_relative_to_line_to_join = crossProduct_to_find_arc_sign(3);

% If everything is done right, signed distance along the line to circle
% is always positive
if signed_distance_along_line_to_join<0
    error('Stop here');
end

% Depending on which side of the line the arc is located, either need
% to subtract or add on the circle radius. Note: the result is designed
% to be positive if the line doesn't quite meet the circle as a tangent
% creating a gap - in this case, a spiral is needed. The result is negative
% if the line would intersect the circle; in this case, a line shift offset
% is needed.
if signed_distance_ortho_line_to_arc_center >= 0
    offset_dist_from_line_toward_circle   = signed_distance_ortho_line_to_arc_center - arc_radius;
else
    offset_dist_from_line_toward_circle   = -1*(signed_distance_ortho_line_to_arc_center + arc_radius);
end
offset = offset_dist_from_line_toward_circle;

% Do we join the line to the arc?
if abs(offset)<threshold
    % Join the line to the arc with a direct connection, no spiral

    % Need to find where the line and arc join
    join_point_xy = line_base_point_xy + signed_distance_along_line_to_join*to_joint_line_unit_tangent_vector;


    % Calculate how much shift is needed to connect the line exactly to the
    % arc
    vector_from_arc_center_to_join = join_point_xy - arc_center_xy;
    unit_vector_from_arc_center_to_join = fcn_geometry_calcUnitVector(vector_from_arc_center_to_join);
    point_shift_xy = -1*offset_dist_from_line_toward_circle*unit_vector_from_arc_center_to_join;

    %%% Fix the line
    if 0==flag_arc_is_first
        % The line is first, do not shift it
        new_line_unit_tangent_vector = to_joint_line_unit_tangent_vector;
        new_line_base_point_xy       = line_base_point_xy;
        new_line_s_start             = 0;
        new_line_s_end               = signed_distance_along_line_to_join;
        if arc_direction_relative_to_line_to_join>=0
            % Arc veers to the left relative to the line's vector
            new_arc_is_counter_clockwise = 1;
        else
            % Arc veers to the right relative to the line's vector
            new_arc_is_counter_clockwise = 0;
        end

    else
        % The line is last, need to shift it
        new_line_unit_tangent_vector = -1*to_joint_line_unit_tangent_vector;
        new_line_base_point_xy       = join_point_xy + point_shift_xy;
        new_line_s_start             = 0;
        new_line_s_end               = signed_distance_along_line_to_join;
        if arc_direction_relative_to_line_to_join>=0
            new_arc_is_counter_clockwise = 0;
        else
            new_arc_is_counter_clockwise = 1;
        end
    end
    revised_line_parameters(1,1:2) = new_line_unit_tangent_vector;
    revised_line_parameters(1,3:4) = new_line_base_point_xy;
    revised_line_parameters(1,5)   = new_line_s_start;
    revised_line_parameters(1,6)   = new_line_s_end;

    %%% Fix the arc
    % The angle the arc starts at the join point
    arc_angle_at_join = atan2(vector_from_arc_center_to_join(1,2),vector_from_arc_center_to_join(1,1));

    if 0==flag_arc_is_first
        % The line is first, so need to shift the arc
        new_arc_center_xy = arc_center_xy - point_shift_xy;
        new_arc_start_angle_in_radians = arc_angle_at_join;
        new_arc_end_angle_in_radians   = arc_angle_at_join + change_in_arc_angle;        
    else
        % The line is last, do not shift the arc
        new_arc_center_xy = arc_center_xy;
        new_arc_end_angle_in_radians   = arc_angle_at_join;
        new_arc_start_angle_in_radians = arc_angle_at_join + change_in_arc_angle;        
    end

    revised_arc_parameters(1,1:2) = new_arc_center_xy;
    revised_arc_parameters(1,3)   = arc_radius;
    revised_arc_parameters(1,4)   = new_arc_start_angle_in_radians;
    revised_arc_parameters(1,5)   = new_arc_end_angle_in_radians;
    revised_arc_parameters(1,6)   = arc_is_circle;
    revised_arc_parameters(1,7)   = new_arc_is_counter_clockwise;

else
    warning('on','backtrace');
    warning('An error will now be thrown due to detection of a needed spiral.');
    error('The use of spirals to join lines and arcs is not coded yet')

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


    subplot(1,2,2);
    hold on;
    grid on;
    axis equal;
    title('Revised');
    xlabel('X [meters]');
    ylabel('Y [meters]')

    % Plot the outputs
    fcn_geometry_plotGeometry('segment',revised_line_parameters);
    fcn_geometry_plotGeometry('arc',revised_arc_parameters);

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
%% fcn_geometry_plotGeometry
function fcn_geometry_plotGeometry(plot_type_string, parameters)

% Get the color vector using the name
color_vector = fcn_geometry_fillColorFromNumberOrName(2,lower(plot_type_string));

% Change the plot style depending on type
switch lower(plot_type_string)
    case {'line','segment'}
        line_vector          = parameters(1,1:2);
        base_point_xy        = parameters(1,3:4);
        station_distance_min = parameters(1,5);
        station_distance_max = parameters(1,6);
        
        stations = [station_distance_min; station_distance_max];
        XY_data = stations*line_vector + ones(length(stations),1)* base_point_xy;
        plot(XY_data(:,1),XY_data(:,2),'-','LineWidth',3,'Color',color_vector);

    case 'arc'
        circleCenter         = parameters(1,1:2);
        circleRadius         = parameters(1,3);
        arcAngles            = parameters(1,4:5);
        flag_arc_is_counterclockwise = parameters(1,7);

        start_angle_in_radians = arcAngles(1);
        end_angle_in_radians   = arcAngles(2);
        degree_step = [];

        XY_data = fcn_geometry_plotArc(circleCenter, circleRadius, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step),[],[]);


        plot(XY_data(:,1),XY_data(:,2),'-','LineWidth',3,'Color',color_vector);

    otherwise
        warning('on','backtrace');
        warning('An error will now be thrown because a geometry string was not recognized.');
        error('Unknown plotting type: %s', plot_type_string);
end

% Plot green headers - calculated from vector direction
vector_direction_start = XY_data(2,1:2) - XY_data(1,1:2);
start_length = sum(vector_direction_start.^2,2).^0.5;
unit_vector_direction_start = fcn_geometry_calcUnitVector(vector_direction_start);
arrow_length = max(min(10,start_length*0.2),0.2);
offset_start = XY_data(1,1:2) + arrow_length*unit_vector_direction_start;

start_line = [XY_data(1,1:2) 0; offset_start, 0];
plot(start_line(:,1),start_line(:,2), '-','Color',[0 1 0],'Linewidth',5);

% Plot red tailers - calculated from vector direction
vector_direction_end = (XY_data(end,1:2) - XY_data(end-1,1:2));
end_length = sum(vector_direction_end.^2,2).^0.5;
unit_vector_direction_end = fcn_geometry_calcUnitVector(vector_direction_end);
arrow_length = max(min(10,end_length*0.2),0.2);
offset_end = XY_data(end,1:2) - arrow_length*unit_vector_direction_end;
end_line = [offset_end, 0; XY_data(end,1:2) 0];
plot(end_line(:,1),end_line(:,2), '-','Color',[1 0 0],'Linewidth',5);

end % Ends fcn_geometry_plotGeometry

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
if closest_pair<=2
    % Line segment is pointing away from junction, need to fix orientation,
    % base point, and station
    corrected_line_unit_tangent_vector = -line_unit_tangent_vector;
    corrected_line_base_point_xy       = line_end_xy;
    corrected_line_s_start             = 0;
    corrected_line_s_end               = line_s_end - line_s_start;
else
    % Line segment is pointint into junction, need to just fix base point and station
    corrected_line_unit_tangent_vector = line_unit_tangent_vector;
    corrected_line_base_point_xy       = line_start_xy;
    corrected_line_s_start             = 0;
    corrected_line_s_end               = line_s_end - line_s_start;    
end
clean_line_parameters(1,1:2)  = corrected_line_unit_tangent_vector;
clean_line_parameters(1,3:4)  = corrected_line_base_point_xy;
clean_line_parameters(1,5:6)  = [corrected_line_s_start corrected_line_s_end];

% Fix the arc
if 0==mod(closest_pair,2)
    % In these cases, the arc is entering the junction at its end. This is
    % not correct. Need to "flip" the arc's orientation.
    if 1==arc_is_counter_clockwise
        corrected_arc_is_counter_clockwise = 0;
    else
        corrected_arc_is_counter_clockwise = 1;
    end
    corrected_arc_start_angle_in_radians = atan2(arc_end_unit_vector(2),arc_end_unit_vector(1));
    corrected_arc_end_angle_in_radians   = corrected_arc_start_angle_in_radians - change_in_arc_angle;
else
    % In these cases, the arc is leaving the junction at its start. This is
    % correct so just pass through the variables
    corrected_arc_is_counter_clockwise = arc_is_counter_clockwise;
    corrected_arc_start_angle_in_radians = atan2(arc_start_unit_vector(2),arc_start_unit_vector(1));
    corrected_arc_end_angle_in_radians   = corrected_arc_start_angle_in_radians + change_in_arc_angle;
end
clean_arc_parameters(1,1:2)   = arc_parameters(1,1:2); % center of the arc does not change
clean_arc_parameters(1,3)     = arc_parameters(1,3);   % radius of the arc does not change
clean_arc_parameters(1,4:5)   = [corrected_arc_start_angle_in_radians corrected_arc_end_angle_in_radians];
clean_arc_parameters(1,6)     = arc_parameters(1,6);   % flag is circle
clean_arc_parameters(1,7)     = corrected_arc_is_counter_clockwise;


end % ends fcn_INTERNAL_fixOrientationAndOrdering