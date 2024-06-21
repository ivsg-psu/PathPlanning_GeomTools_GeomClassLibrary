function [revised_arc_parameters, revised_line_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcLine( arc_parameters, line_parameters, varargin)
%% fcn_geometry_alignArcLine
% Revises the geometric parameters of an arc and line such that
% they align where they join. It does this by checking the offset between
% the two objects at the join location. 
% 
% If the alignment is not feasible but the offset is less than a threshold
% (default is 0.1 meter), the line's geometric position is shifted to force
% alignment with the arc. In other words, the arc is kept stationary and
% the line is aligned to the arc.
%
% Alignment types are allowed of different types of continuity, including:
%
% C0 continuity which checks for one intersection and joins at that point; 
%
% C1 continuity which checks for a tangent point, and joins at that point;
%
% C2 continuity which creates a spiral connection between the geometries.
%
% Format:
% [revised_line_parameters, revised_arc_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters]  = ...
% fcn_geometry_alignArcLine(arc_parameters, line_parameters, (threshold), (continuity_level),  (fig_num))
%
% INPUTS:
%
%      line_parameters: the parameter set describing the line 
%      geometry. See fcn_geometry_fillEmptyDomainStructure for details,
%      specifically the structure for 'Vector regression line fit'.
%
%      arc_parameters: the parameter set describing the arc geometry. See
%      fcn_geometry_fillEmptyDomainStructure for details, specifically the
%      structure for 'Regression arc'.
%
%      (OPTIONAL INPUTS)
%
%      threshold: the offset, in meters, between the arc and the line such
%      that this offset is removed by shifting. If the offset is larger
%      than this, then the outputs are set to empty.If this is entered as a
%      2x1 or 1x2, then this specifies the threshold in St coordinated,
%      namely first in the station direction, and then in the transverse
%      direction. For example, an entry of [3 0.02] would have 3 meters
%      threshold in the station direction, but 0.02 meters threshold in the
%      transverse direction.
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
%      revised_arc_parameters: the parameter set describing the arc
%      geometry that joins the geometries. See
%      fcn_geometry_fillEmptyDomainStructure for details, specifically the
%      structure for 'Regression arc'.
%
%      revised_line_parameters: the parameter set describing the line
%      geometry that joins the geometries. See
%      fcn_geometry_fillEmptyDomainStructure for details, specifically the
%      structure for 'Vector regression line fit'.
%
%      revised_intermediate_geometry_join_type: for type C2 continuity, an
%      intermediate geometry is often inserted in the form of a spiral.
%      This output saves the geometry type as a string type.
%
%      revised_intermediate_geometry_join_parameters: the parameter set describing the
%      spiral geometry that joins the line and arc geometries if C2
%      continuity is specified. See fcn_geometry_fillEmptyDomainStructure
%      for details, specifically the structure for 'spiral'.
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_alignArcLine
% for a full test suite.
%
% This function was written on 2024_04_12 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_04_12 - Sean Brennan
% -- wrote the code
% 2024_04_19 - Sean Brennan
% -- renamed from fcn_geometry_joinLineToArc
% -- fixed bug where calculation still works if error larger than tolerance
% -- added continuity_level input
% 2024_04_20 - Sean Brennan
% -- added St conversion functions
% -- added powerful debugging plots (VERY useful - caught lots of mistakes)
% -- finished functionalizing code
% -- added C0 and C1 continuity, confirmed via script testing they work
% -- added C2 continuity and revised_spiral_join_parameters output
% -- bug fix in nargin check
% 2024_05_10 - Sean Brennan
% -- changed output list to match arc to arc alignment code
% -- removed arc is first flag
% -- functionalized code to match arc to arc
% -- changed code to force line to arc functionality only (per name of fcn)
% -- renamed function to ArcLine because LineToArc was confusing as to
%    which was first
% 2024_05_28 - S. Brennan
% -- fixed call to spiralFromCircleToCircle to use parameter vectors
% -- fixed subplot axes to match each other
% 2024_06_16 - Sean Brennan
% -- changed parameter format to new style:
%            'spiral' - 
%               [
%                x0,  % The initial x value
%                y0,  % The initial y value
%                h0,  % The initial heading
%                s_Length,  % the s-coordinate length allowed
%                K0,  % The initial curvature
%                Kf   % The final curvature
%              ] 
% 2024_06_19 - Sean Brennan
% -- changed segment references to lines 
% -- changed parameter format for line to new standard:
%             [
%              base_point_x, 
%              base_point_y, 
%              heading,
%             ]
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

% Does user want to specify threshold (this in turn turns on/off flag_perform_shift_of_line)?
threshold = 0.1;
flag_perform_shift_of_line = 1;
if (3<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        threshold = temp;
    else
        flag_perform_shift_of_line = 0;
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

%% Plot inputs?
fcn_INTERNAL_prepDebugFigure(arc_parameters, line_parameters, debug_fig_num);

%% Check to see if arc and line intersect
intersection_point1 = fcn_INTERNAL_ArcLineIntersection(arc_parameters, line_parameters, 1, debug_fig_num);
if (length(intersection_point1(:,1))>1)&&(0==continuity_level)
    warning('Multiple intersection points found between a line and arc geometry with a requested connection type of C0 continuity. Unable to resolve which intersection to use.');
end

%% Rearrange parameters so line is always the 1st input, line is 2nd
% Fix the parameters to make the line first, arc second, and make
% sure the line and arc point into and then out of the junction
% respectively. For situations where arc is actually the first input, this
% is fixed in later steps using a flag.
[clean_arc_parameters, clean_line_parameters] = fcn_INTERNAL_fixOrientationAndOrdering(arc_parameters, line_parameters, intersection_point1, debug_fig_num);


%% Get new intersection point, if arcs changed shape
intersection_point2 = fcn_INTERNAL_ArcLineIntersection(clean_arc_parameters,clean_line_parameters, 2, debug_fig_num);

%% Rotate the geometries out of XY into ST coordinates
% so that the tangent line is oriented horizontally
% and the start of the tangent line on arc1 is at the origin.
% This is to make the debugging MUCH easier, as it reduces permutations.
% Again, this is fixed in later steps.
[st_arc_parameters, st_line_parameters, St_transform_XYtoSt, flag_arc1_is_flipped] = ...
    fcn_INTERNAL_convertParametersToStOrientation(clean_arc_parameters, clean_line_parameters, continuity_level, intersection_point2, debug_fig_num);

%% Check how much shift is needed to connect line to arc
[desired_st_arc_parameters, desired_st_line_parameters, desired_st_intermediate_geometry_join_parameters, desired_intermediate_geometry_join_type] = ...
    fcn_INTERNAL_findShiftToMatchLineToArc(st_arc_parameters, st_line_parameters, continuity_level, intersection_point2, threshold, flag_perform_shift_of_line, debug_fig_num);
% Deltas are from desired to actual

%% Perform shift to join arc and line
[revised_arc_parameters_St,revised_line_parameters_St, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters_St] = ...
    fcn_INTERNAL_performShift(threshold, continuity_level, ...
    st_arc_parameters, st_line_parameters, ...
    desired_st_arc_parameters, desired_st_line_parameters, ...
    desired_st_intermediate_geometry_join_parameters, desired_intermediate_geometry_join_type, debug_fig_num);

%% Rotate results out of St back into XY
[revised_arc_parameters, revised_line_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_INTERNAL_convertParametersOutOfStOrientation(...
    revised_arc_parameters_St, revised_line_parameters_St, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters_St, St_transform_XYtoSt, flag_arc1_is_flipped, debug_fig_num);


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
    sgtitle('Arc joining with line');
    title('Original');
    xlabel('X [meters]');
    ylabel('Y [meters]')

    % Plot the inputs
    fcn_geometry_plotGeometry('line',line_parameters);
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

    % Plot the outputs
    fcn_geometry_plotGeometry('line',revised_line_parameters);
    fcn_geometry_plotGeometry('arc',revised_arc_parameters);
    fcn_geometry_plotGeometry('spiral',revised_intermediate_geometry_join_parameters);


    % % Make axis slightly larger?
    % if flag_rescale_axis
    %     temp = axis;
    %     %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
    %     axis_range_x = temp(2)-temp(1);
    %     axis_range_y = temp(4)-temp(3);
    %     percent_larger = 0.3;
    %     temp_axis2 = [temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y];
    %     good_axis = [min([temp_axis(1) temp_axis2(1)])  max([temp_axis(2) temp_axis2(2)])  min([temp_axis(3) temp_axis2(3)])  max([temp_axis(4) temp_axis2(4)]) ];
    %     axis(good_axis);
    % end

    % Force the subplots to have matching axes
    % subplot(1,2,2);
    % good_axis = axis;
    % 
    % subplot(1,2,1);
    axis(temp_axis)

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

%% fcn_INTERNAL_prepDebugFigure
function fcn_INTERNAL_prepDebugFigure(arc_parameters, line_parameters, debug_fig_num)
if ~isempty(debug_fig_num)
    figure(debug_fig_num);
    clf;

    % Plot the inputs
    subplot(3,2,1);


    fcn_geometry_plotCircle(arc_parameters(1,1:2),arc_parameters(1,3),...
        sprintf(' ''--'',''Color'',[0 0.6 0],''LineWidth'',1 '),debug_fig_num);

    fcn_geometry_plotGeometry('arc',arc_parameters);

    temp = axis;
    %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    new_axis = [temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y];
    axis([min(new_axis) max(new_axis) min(new_axis) max(new_axis)]);

    fcn_geometry_plotGeometry('line',line_parameters);

    hold on;
    grid on;
    axis equal;
    xlabel('X [meters]');
    ylabel('Y [meters]')

    title('Inputs');

end
end % Ends fcn_INTERNAL_prepDebugFigure

%% fcn_INTERNAL_ArcSegmentIntersection
function  intersection_point_arc_to_line = fcn_INTERNAL_ArcLineIntersection(arc_parameters,line_parameters, subplot_number, debug_fig_num)
firstFitType = 'arc';
firstFitType_parameters = arc_parameters;
secondFitType = 'line';
secondFitType_parameters = line_parameters;

intersection_point_arc_to_line = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, -1);


if ~isempty(debug_fig_num)
    % Plot the intersection
    figure(debug_fig_num);
    subplot(3,2,1);
    debug_axis = axis;

    subplot(3,2,subplot_number);
    plot(intersection_point_arc_to_line(:,1),intersection_point_arc_to_line(:,2),'co','MarkerSize',10,'LineWidth',2);

    axis(debug_axis);
end

end % Ends fcn_INTERNAL_ArcSegmentIntersection

%% fcn_INTERNAL_fixOrientationAndOrdering
function [clean_arc_parameters, clean_line_parameters] = fcn_INTERNAL_fixOrientationAndOrdering(arc_parameters, line_parameters, ~, debug_fig_num)
% This function takes the parameter inputs and produces parameter sets such
% that the arc is first, it is oriented so that it ends at the junction
% with the line, and line is modified so it starts at or near the
% junction. It also forces the line to start at a station value of 0.

% Get the arc fit details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy                = arc_parameters(1,1:2);
arc_radius                   = arc_parameters(1,3);
arc_start_angle_in_radians   = arc_parameters(1,4);
arc_end_angle_in_radians     = arc_parameters(1,5);
arc_start_xy                 = arc_center_xy + arc_radius*[cos(arc_start_angle_in_radians) sin(arc_start_angle_in_radians)];
arc_end_xy                   = arc_center_xy + arc_radius*[cos(arc_end_angle_in_radians) sin(arc_end_angle_in_radians)];
% arc_is_counter_clockwise     = arc_parameters(1,7);

% % Find the change in angle of the arc
% arc_start_unit_vector        = [cos(arc_start_angle_in_radians) sin(arc_start_angle_in_radians)];
% arc_end_unit_vector          = [cos(arc_end_angle_in_radians)   sin(arc_end_angle_in_radians)  ];
% if arc_is_counter_clockwise
%     cross_product_direction = 1;
% else
%     cross_product_direction = -1;
% end
% arc_change_in_angle = fcn_geometry_findAngleUsing2PointsOnCircle([0 0],1, arc_start_unit_vector, arc_end_unit_vector, cross_product_direction);



% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_base_point_xy       = line_parameters(1,1:2);
% line_unit_tangent_vector = [cos(line_parameters(1,3)) sin(line_parameters(1,3))];
line_unit_ortho_vector = [cos(line_parameters(1,3)+pi/2) sin(line_parameters(1,3)+pi/2)];




% Project from the arc to the line
vector_from_line_base_to_arc_start = arc_start_xy - line_base_point_xy;
vector_from_line_base_to_arc_end   = arc_end_xy - line_base_point_xy;

distance_arcStart_to_line = abs(sum(line_unit_ortho_vector.*vector_from_line_base_to_arc_start,2));
distance_arcEnd_to_line   = abs(sum(line_unit_ortho_vector.*vector_from_line_base_to_arc_end,2));


clean_line_parameters = line_parameters;

% Fix the arc or line depending on which combo is closest. The end of the
% arc should be the one with the least distance to the line. If this is not
% the case, then we neet to flip the orientation of the arc.
if distance_arcStart_to_line<distance_arcEnd_to_line
    % Flip the arc
    clean_arc_parameters = fcn_geometry_flipGeom('arc',  arc_parameters, -1);
else
    clean_arc_parameters = arc_parameters;
end

% % Make the arc end at the intersection, if it's there
% if ~any(isnan(intersection_point))
%     arc_endPoint = intersection_point;
% end
% 
% % Set the outputs
% 
% clean_arc_parameters(1,1:2)   = arc_parameters(1,1:2); % center of the arc does not change
% clean_arc_parameters(1,3)     = arc_parameters(1,3);   % radius of the arc does not change
% clean_arc_parameters(1,4:5)   = [corrected_arc_start_angle_in_radians corrected_arc_end_angle_in_radians];
% clean_arc_parameters(1,6)     = arc_parameters(1,6);   % flag is circle
% clean_arc_parameters(1,7)     = corrected_arc_is_counter_clockwise;
% 
% % Set the line start and end by distances only, starting at 0 and ending at
% % total distance
% clean_line_parameters(1,1:2)  = corrected_line_base_point_xy;
% clean_line_parameters(1,3)    = atan2(corrected_line_unit_tangent_vector(2),corrected_line_unit_tangent_vector(1));

if ~isempty(debug_fig_num)
    figure(debug_fig_num);
    subplot(3,2,1);
    debug_axis = axis;

    % Plot the cleaned inputs
    subplot(3,2,2);

    fcn_geometry_plotCircle(clean_arc_parameters(1,1:2),clean_arc_parameters(1,3),...
        sprintf(' ''--'',''Color'',[0 0.6 0],''LineWidth'',1 '),debug_fig_num);

    fcn_geometry_plotGeometry('arc',clean_arc_parameters);
    fcn_geometry_plotGeometry('line',clean_line_parameters);

    axis(debug_axis);
    hold on;
    grid on;
    axis equal;
    xlabel('X [meters]');
    ylabel('Y [meters]')

    title('Corrected inputs');
end


end % ends fcn_INTERNAL_fixOrientationAndOrdering

%% fcn_INTERNAL_convertParametersToStOrientation
function [st_arc_parameters, st_line_parameters, St_transform_XYtoSt, flag_arc_is_flipped] = ...
    fcn_INTERNAL_convertParametersToStOrientation(arc_parameters, line_parameters, continuity_level, intersection_point, debug_fig_num)

% Calculate needed values from parameter sets
% Get the arc fit details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy                = arc_parameters(1,1:2);
% arc_radius                   = arc_parameters(1,3);
% arc1_start_angle_in_radians   = arc1_parameters(1,4);
arc_end_angle_in_radians     = arc_parameters(1,5);
% arc_is_counter_clockwise     = arc_parameters(1,7);
% arc1_start_xy                 = arc1_center_xy + arc1_radius*[cos(arc1_start_angle_in_radians) sin(arc1_start_angle_in_radians)];
% arc_end_xy                   = arc_center_xy + arc_radius*[cos(arc_end_angle_in_radians) sin(arc_end_angle_in_radians)];

% % Find the change in angle of the arc
% % arc1_start_unit_vector        = [cos(arc1_start_angle_in_radians) sin(arc1_start_angle_in_radians)];
% % arc1_end_unit_vector          = [cos(arc1_end_angle_in_radians)   sin(arc1_end_angle_in_radians)  ];
% if arc_is_counter_clockwise
%     arc_cross_product_direction = 1;
% else
%     arc_cross_product_direction = -1;
% end
% arc1_change_in_angle = fcn_geometry_findAngleUsing2PointsOnCircle([0 0],1, arc1_start_unit_vector, arc1_end_unit_vector, cross_product_direction);


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_base_point_xy           = line_parameters(1,1:2);
line_unit_tangent_vector     = [cos(line_parameters(1,3)) sin(line_parameters(1,3))];
% line_s_start                 = line_parameters(1,5);
% line_s_end                   = line_parameters(1,6);
% line_angle = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));
% line_unit_ortho_vector = line_unit_tangent_vector*[0 1; -1 0];

% Find the closest point of the arc to the line
vector_from_line_base_to_arc_center = arc_center_xy - line_base_point_xy;
distance_along_line_to_circle_closest_point = sum(line_unit_tangent_vector.*vector_from_line_base_to_arc_center,2);
closest_point_on_line_to_circle_center = line_base_point_xy + distance_along_line_to_circle_closest_point*line_unit_tangent_vector;

% % Calculate the distance from the arc to the line
% vector_from_arc_center_to_line_joint = line_base_point_xy - arc_center_xy;
% if arc_is_counter_clockwise
%     distance_from_arc_center_to_line = sum(line_unit_ortho_vector.*vector_from_arc_center_to_line_joint,2);
% else
%     distance_from_arc_center_to_line = -sum(line_unit_ortho_vector.*vector_from_arc_center_to_line_joint,2);
% end

switch continuity_level
    case 0
        % For C0 continuity of arc to line, the closest point of the
        % desired joint to the arc is either the intersection of arc with
        % the line or where the arc ends

        % Was an intersection found? If not, use the end of the arc as the
        % final point
        if any(isnan(intersection_point))
            % If enter here, no intersection was found, use the end of the
            % arc as the actual end
            desired_angle_arc_end = arc_end_angle_in_radians;
        else
            % If enter here, an intersection was found

            %%%%
            % Find the angle where arc should end
            arc_vector_center_to_intersection = intersection_point - arc_center_xy;

            % Plot the vector (for debugging)?
            if ~isempty(debug_fig_num)
                figure(debug_fig_num);
                subplot(3,2,1);

                % Plot the projection vector
                subplot(3,2,2);
                quiver(arc_center_xy(1,1), arc_center_xy(1,2), arc_vector_center_to_intersection(1,1),  arc_vector_center_to_intersection(1,2), 0,'LineWidth',3);
            end

            desired_angle_arc_end = atan2(arc_vector_center_to_intersection(2),arc_vector_center_to_intersection(1));

        end

    case 1
        % For C1 continuity of the arc to the line, the arc and line
        % are tangent to each other. So the rotation will need to be the
        % one that produces the arc such that it is tangent with the x-axis at
        % exactly the point where the arc touches the tangent line to the
        % segement. In some cases, this will not exist and in these cases
        % the equivalent line must be found.

        % First, calculate all the tangent point. To do this,
        % project from the center of the arc orthogonal to the line
        % by the radius

        % Find the angle where arc1 should end

        % if arc_is_counter_clockwise
        %     desired_arc_end_position = arc_center_xy - line_unit_ortho_vector*arc_radius;
        % else
        %     desired_arc_end_position = arc_center_xy + line_unit_ortho_vector*arc_radius;
        % end
        desired_arc_end_position = closest_point_on_line_to_circle_center;
        arc_vector_center_to_desired_arc_end = desired_arc_end_position - arc_center_xy;

        % Plot the vector (for debugging)?
        if ~isempty(debug_fig_num)
            figure(debug_fig_num);
            subplot(3,2,1);

            % Plot the projection vector
            subplot(3,2,2);
            quiver(arc_center_xy(1,1), arc_center_xy(1,2), arc_vector_center_to_desired_arc_end(1,1),  arc_vector_center_to_desired_arc_end(1,2), 0,'LineWidth',3);
        end

        desired_angle_arc_end = atan2(arc_vector_center_to_desired_arc_end(2),arc_vector_center_to_desired_arc_end(1));

    case 2
        % For C2 continuity of an arc to a line, this will involve a
        % spiral.  For now, just use the arc to align
        desired_angle_arc_end = arc_end_angle_in_radians;

    otherwise
        error('This continuity not possible yet')
end

% Perform the rotation
desired_arc1_parameters = arc_parameters;
desired_arc1_parameters(1,5) = desired_angle_arc_end;

secondary_parameters_type_strings{1} = 'arc';
secondary_parameters{1}              = arc_parameters;
secondary_parameters_type_strings{2} = 'line';
secondary_parameters{2}              = line_parameters;

[~, st_secondary_parameters, St_transform_XYtoSt, ~, flag_arc_is_flipped] = ...
    fcn_geometry_orientGeometryXY2St('arc', desired_arc1_parameters, (secondary_parameters_type_strings), (secondary_parameters), (-1));

st_arc_parameters = st_secondary_parameters{1};
st_line_parameters = st_secondary_parameters{2};

if ~isempty(debug_fig_num)
    figure(debug_fig_num);
    subplot(3,2,1);
    debug_axis = axis;

    % Plot the rotated inputs
    subplot(3,2,3);
    hold on;
    grid on;
    axis equal;
    xlabel('S [meters]');
    ylabel('t [meters]')

    fcn_geometry_plotCircle(st_arc_parameters(1,1:2),st_arc_parameters(1,3),...
        sprintf(' ''--'',''Color'',[0 0.6 0],''LineWidth'',1 '),debug_fig_num);

    fcn_geometry_plotGeometry('arc',st_arc_parameters);
    fcn_geometry_plotGeometry('line',st_line_parameters);

    plot(st_arc_parameters(1,1),st_arc_parameters(1,2),'+','Color',[0 0.6 0]);

    axis(debug_axis);


    title('Rotated into St');
end
end % Ends fcn_INTERNAL_convertParametersToStOrientation





%% fcn_INTERNAL_findShiftToMatchSegmentToArc
function [desired_arc_parameters, desired_line_parameters, desired_intermediate_geometry_join_parameters, desired_intermediate_geometry_join_type] = ...
    fcn_INTERNAL_findShiftToMatchLineToArc(arc_parameters, line_parameters, continuity_level, intersection_point, threshold, flag_perform_shift_of_line, debug_fig_num)
% Calculates the delta amount to match the line to the arc. The delta
% values are measured FROM desired point TO actual point


% Calculate needed values from parameter sets

% Get the arc fit details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
%               [circleCenter_x.
%                circleCenter_y,
%                radius,
%                start_angle_in_radians, 
%                end_angle_in_radians,
%                flag_this_is_a_circle
%                flag_arc_is_counterclockwise
%               ] 

arc_center_xy                = arc_parameters(1,1:2);
arc_radius                   = arc_parameters(1,3);
arc_start_angle_in_radians   = arc_parameters(1,4);
arc_end_angle_in_radians     = arc_parameters(1,5);
arc_start_xy                 = arc_center_xy + arc_radius*[cos(arc_start_angle_in_radians) sin(arc_start_angle_in_radians)];
arc_end_xy                   = arc_center_xy + arc_radius*[cos(arc_end_angle_in_radians) sin(arc_end_angle_in_radians)];
arc_is_counter_clockwise     = arc_parameters(1,7);

% Find the change in angle of the arc
% arc1_start_unit_vector        = [cos(arc1_start_angle_in_radians) sin(arc1_start_angle_in_radians)];
% arc1_end_unit_vector          = [cos(arc1_end_angle_in_radians)   sin(arc1_end_angle_in_radians)  ];
if 1==arc_is_counter_clockwise
    arc_is_counter_clockwise = 1;
else
    arc_is_counter_clockwise = -1;
end
% arc1_change_in_angle = fcn_geometry_findAngleUsing2PointsOnCircle([0 0],1, arc1_start_unit_vector, arc1_end_unit_vector, arc1_is_counter_clockwise);



% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_base_point_xy           = line_parameters(1,1:2);
line_unit_tangent_vector     = [cos(line_parameters(1,3)) sin(line_parameters(1,3))];
% line_angle = atan2(line_unit_tangent_vector(2),line_unit_tangent_vector(1));
% line_unit_ortho_vector = line_unit_tangent_vector*[0 1; -1 0];

% Find the closest point of the arc to the line
vector_from_line_base_to_arc_center = arc_center_xy - line_base_point_xy;
distance_along_line_to_circle_closest_point = sum(line_unit_tangent_vector.*vector_from_line_base_to_arc_center,2);
closest_point_on_line_to_circle_center = line_base_point_xy + distance_along_line_to_circle_closest_point*line_unit_tangent_vector;

distance_from_arc_center_to_line = sum((closest_point_on_line_to_circle_center -arc_center_xy).^2,2).^0.5; 
% % Calculate the distance from the arc to the line
% vector_from_arc_center_to_line_joint = line_base_point_xy - arc_center_xy;
% if arc_is_counter_clockwise
%     distance_from_arc_center_to_line = -sum(line_unit_ortho_vector.*vector_from_arc_center_to_line_joint,2);
% else
%     distance_from_arc_center_to_line = sum(line_unit_ortho_vector.*vector_from_arc_center_to_line_joint,2);
% end


% Calculate the distance between the circles and the join point
space_from_arc_to_line = distance_from_arc_center_to_line - arc_radius;


% Check that the threshold is given in St format or not. If in St format,
% just use the t portion.
if length(threshold)>1
    t_threshold = threshold(2);
else
    t_threshold = threshold;
end


% With the cleaned parameters, the line vector always
% points in the direction of the join point.
% to_joint_unit_tangent_vector = [1 0];
% to_joint_unit_ortho_vector   = [0 1];

% For arc connections to lines, the only allowable intermediate geometry is
% a spiral
desired_intermediate_geometry_join_type       = 'spiral'; % Intermediate geometry will be a spiral
desired_intermediate_geometry_join_parameters = nan(1,6);


switch continuity_level
    case 0
        % For C0 continuity of arc to line, the closest point of the
        % desired joint to the arc is simply the end of arc

        if ~any(isnan(intersection_point))
            % Arc will end at the intersection, which is at the origin
            vector_from_arc_center_to_line_joint = [0 0] - arc_center_xy;
            angle_of_intersection = atan2(vector_from_arc_center_to_line_joint(2),vector_from_arc_center_to_line_joint(1));

            desired_arc_parameters        = arc_parameters;
            desired_arc_parameters(1,5)   = angle_of_intersection; % Update where the spiral ends

            desired_line_parameters        = line_parameters;
            desired_line_parameters(1,1:2) = [0 0]; % Move line_base_point_xy

        else
            if flag_perform_shift_of_line==1
                if abs(space_from_arc_to_line)<=t_threshold
                    % Yes, line can be moved enough to be tangent. So put
                    % line's start at the end of the arc. Do not change the
                    % arc.
                    desired_arc_parameters        = arc_parameters;

                    desired_line_parameters        = line_parameters;
                    desired_line_parameters(1,1:2) = [0 0]; % Move line_base_point_xy

                else
                    % Not possible to shift enough to allow connection
                    desired_arc_parameters        = nan(size(arc_parameters));
                    desired_line_parameters    = nan(1,3);
                end
            else
                % No shift allowed by user entry, so not possible
                desired_arc_parameters            = nan(size(arc_parameters));
                desired_line_parameters        = nan(1,3);
            end
        end

    case 1
        % For C1 continuity of arc to line, the closest point of the
        % desired joint to the arc is always going to be the point where
        % the arc touches the tangent line of the line. Since the
        % line is the only one that can be moved, the connecting point
        % for the arc is simply the location where the arc is tangent, which by
        % construction is [0 0]. The connection point for the line is at
        % this same point


        if flag_perform_shift_of_line==1
            if abs(space_from_arc_to_line)<=t_threshold
                % Yes, arc2 can be moved enough to be tangent. So put
                % arc2's center in correct place
                desired_arc_parameters        = arc_parameters;
                desired_arc_parameters(1,5)   = -pi/2;

                desired_line_parameters        = line_parameters;
                desired_line_parameters(1,1:2) = [0 0]; % Move line_base_point_xy

            else
                % Not possible to shift enough to allow connection
                desired_arc_parameters        = nan(size(arc_parameters));
                desired_line_parameters    = nan(1,3);
            end
        else
            % No shift allowed by user entry, so not possible
            desired_arc_parameters            = nan(size(arc_parameters));
            desired_line_parameters        = nan(1,3);
        end


    case 2
        
        % For C2 continuity of a line to an arc, the arc must connect to the
        % line via a spiral. Calculation of the spiral is difficult and
        % requires numerical iteration, and the inputs require the calculation
        % of the offset of the outer part of the circle from the line. This
        % offset in the direction of the circle center MUST be strictly
        % positive or the spiral cannot work.
        

        flag_spiral_was_calculated = 0;
        if space_from_arc_to_line>0 
            % spiral_join_parameters = [spiralLength,h0,x0,y0,K0,Kf];

            circle1_parameters = [0 arc_radius  arc_radius];
            circle2_parameters = [[0 -space_from_arc_to_line], inf];

            [desired_intermediate_geometry_join_parameters, ~] = ...
                fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, [], -1);

            % [desired_intermediate_geometry_join_parameters, ~] = ...
            %     fcn_geometry_spiralFromCircleToCircle(arc_radius, inf, [0 -space_between_arc_and_line], [], -1);
            flag_spiral_was_calculated = 1;
        else
            % Not enough space between arc and line for a spiral, not
            % without modifying it
            if flag_perform_shift_of_line==1
                if abs(space_from_arc_to_line)<=threshold
                    unit_direction_vector_to_shift_line = [0 -1];

                    revised_line_parameters = line_parameters;
                    revised_line_parameters(1,1:2) = line_parameters(1,1:2) + unit_direction_vector_to_shift_line*(abs(space_from_arc_to_line)+0.001);

                    % Call function again with revised parameters that
                    % should work
                    [desired_arc_parameters, desired_line_parameters, desired_intermediate_geometry_join_parameters] = ...
                        fcn_INTERNAL_findShiftToMatchLineToArc(arc_parameters, revised_line_parameters, continuity_level, intersection_point, [], 0, debug_fig_num);
                    flag_spiral_was_calculated = 1;
                else
                    % Not possible
                    desired_arc_parameters = nan(size(arc_parameters));
                    desired_line_parameters = nan(1,3);
                end
            else
                % Not allowed by user entry
                desired_arc_parameters = nan(size(arc_parameters));
                desired_line_parameters = nan(1,3);
            end
        end % Ends if statement to check if spiral is possible

        % Do we need to fill in the other parameters?
        if 1==flag_spiral_was_calculated
            desired_intermediate_geometry_join_type = 'spiral';

            % Check results?
            if 1==0

                % Call the function fcn_geometry_extractXYfromSTSpiral to predict the
                % spiral and calculate the offsets, plotting the results in
                % figure 1234
                figure(1234);
                clf;
                hold on;
                grid on;
                axis equal;

                fcn_geometry_plotGeometry('arc',arc_parameters);
                fcn_geometry_plotGeometry('line',line_parameters);
                fcn_geometry_plotGeometry(desired_intermediate_geometry_join_type,desired_intermediate_geometry_join_parameters);

            end % Ends plotting

            % Find the angle and position that the spiral starts and ends at
            % spiralLength = desired_intermediate_geometry_join_parameters(1,4);
            h0           = desired_intermediate_geometry_join_parameters(1,3);
            % K0           = desired_intermediate_geometry_join_parameters(1,5);
            % Kf           = desired_intermediate_geometry_join_parameters(1,6);
            % analytical_end_angle   = h0 + (Kf-K0)*spiralLength/2 + K0*spiralLength;


            flag_spiral_is_bad = 0;
            % Make sure the spiral's start position is within the angle range
            % of arc1. If not, return Nan values for everything
            arc_angle_where_spiral_starts = h0 - pi/2;
            spiral_arc_join_xy = arc_center_xy + arc_radius*[cos(arc_angle_where_spiral_starts) sin(arc_angle_where_spiral_starts)];
            intersection_is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(arc_start_xy, spiral_arc_join_xy, arc_end_xy,-1);
            if arc_is_counter_clockwise ~= intersection_is_counterClockwise
                % Spiral does not connect with the arc
                desired_arc_parameters = nan(size(arc_parameters));
                desired_line_parameters = nan(1,3);
                flag_spiral_is_bad = 1;
            end

            % Make sure the spiral's end position is within the distance range of the segement. If not, return Nan values for everything
            % Find the x-value of the spiral result
            XY_spiral = fcn_geometry_plotGeometry(desired_intermediate_geometry_join_type,desired_intermediate_geometry_join_parameters,[],[],-1);
            spiral_end_XY = XY_spiral(end,:);
            if spiral_end_XY(1,1)<0 
                % Spiral does not connect with the line.
                desired_arc_parameters = nan(size(arc_parameters));
                desired_line_parameters = nan(1,3);
                flag_spiral_is_bad = 1;
            end

            if 0==flag_spiral_is_bad
                % If made it here, then the spiral is good. Fill outputs and then
                % return results
                desired_arc_parameters = arc_parameters;
                desired_arc_parameters(1,5)   = arc_angle_where_spiral_starts; % Update where the spiral ends

                desired_line_parameters        = line_parameters;
                desired_line_parameters(1,1:2) = spiral_end_XY; % line_base_point_xy
            else
                % Bad spiral - set all to nan
                desired_intermediate_geometry_join_parameters = nan(size(desired_intermediate_geometry_join_parameters));
            end
        end


        
    otherwise
        error('This continuity not possible yet')
end


if ~isempty(debug_fig_num)
    % Plot the bounding box
    figure(debug_fig_num);
    subplot(3,2,1);
    debug_axis = axis;

    subplot(3,2,4);
    fcn_geometry_plotCircle(desired_arc_parameters(1,1:2),desired_arc_parameters(1,3),...
        sprintf(' ''--'',''Color'',[0 0.6 0],''LineWidth'',1 '),debug_fig_num);

    fcn_geometry_plotGeometry('arc',desired_arc_parameters);
    fcn_geometry_plotGeometry('line',desired_line_parameters);
    fcn_geometry_plotGeometry(desired_intermediate_geometry_join_type,desired_intermediate_geometry_join_parameters);

    axis(debug_axis);

    hold on;
    grid on;
    axis equal;
    xlabel('S [meters]');
    ylabel('t [meters]')

    title('Desired St geometries');
end
end % Ends fcn_INTERNAL_findShiftToMatchSegmentToArc


%% fcn_INTERNAL_performShift
function [revised_arc_parameters_St,revised_line_parameters_St, revised_intermediate_geometry_type, revised_intermediate_geometry_parameters_St] = ...
    fcn_INTERNAL_performShift(St_thresholds, continuity_level, ...
    st_arc_parameters, st_line_parameters, ...
    desired_st_arc_parameters, desired_st_line_parameters, ...
    desired_intermediate_geometry_join_parameters, desired_intermediate_geometry_join_type, ...
    debug_fig_num)

% Check that the threshold is given in St format or not. If in St format,
% just use the t portion.
if length(St_thresholds)>1
    t_threshold = St_thresholds(2);
else
    t_threshold = St_thresholds;
end

% Find out how much line is shifting by looking at how much the base point of
% line is moving in the line's t-direction
st_line_unit_tangent_vector     = [cos(st_line_parameters(1,3)) sin(st_line_parameters(1,3))];
st_line_unit_ortho_vector       = st_line_unit_tangent_vector*[0 1; -1 0];


vector_st_from_desired_line_base_xy_to_original_base_xy  = st_line_parameters(1,1:2)-desired_st_line_parameters(1,1:2);
t_distance_from_desired_line_base_xy_to_original_base_xy = sum(vector_st_from_desired_line_base_xy_to_original_base_xy.*st_line_unit_ortho_vector,2);

% Check to see if shift is even possible
t_shift_distance = abs(t_distance_from_desired_line_base_xy_to_original_base_xy);

if (t_shift_distance>=t_threshold) && (continuity_level~=2)
    % Not possible to shift
    revised_arc_parameters_St                   = nan(size(st_arc_parameters));
    revised_line_parameters_St                  = nan(size(st_line_parameters));
    revised_intermediate_geometry_parameters_St = nan(size(desired_intermediate_geometry_join_parameters));
    revised_intermediate_geometry_type          = desired_intermediate_geometry_join_type; 
else
    revised_arc_parameters_St                   = desired_st_arc_parameters;
    revised_line_parameters_St                  = desired_st_line_parameters;
    revised_intermediate_geometry_parameters_St = desired_intermediate_geometry_join_parameters;
    
    switch continuity_level
        case 0
            % For C0 continuity of arc1 to arc2, the closest point of the
            % desired joint to the arc is either the intersection of arc1
            % arc2, or where arc1 ends
            revised_intermediate_geometry_type = '';
        case 1
            % For C1 continuity of arc1 to arc2, the closest point of the desired
            % joint to arc1 and arc2 is always going to be the point where each arc
            % touches the tangent line connecting them. Since arc2 is the only one
            % that can be moved, the connecting point is simply the location where
            % arc1 is tangent, which by construction is [0 0].
            revised_intermediate_geometry_type = 'line';
        case 2
            revised_intermediate_geometry_type = 'spiral';
        otherwise
            error('This continuity not possible yet')
    end
end
if ~isempty(debug_fig_num)
    % Plot the results
    figure(debug_fig_num);
    subplot(3,2,1);
    debug_axis = axis;

    subplot(3,2,5);

    fcn_geometry_plotGeometry('arc',revised_arc_parameters_St);
    fcn_geometry_plotGeometry('line',revised_line_parameters_St);
    fcn_geometry_plotGeometry(desired_intermediate_geometry_join_type,revised_intermediate_geometry_parameters_St);

    title('St outputs after shift');
    axis(debug_axis);
end

end % Ends fcn_INTERNAL_performShift


%% fcn_INTERNAL_convertParametersOutOfStOrientation
function [revised_arc_parameters, revised_line_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_INTERNAL_convertParametersOutOfStOrientation(...
    revised_arc_parameters_St, revised_line_parameters_St, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters_St, St_transform_XYtoSt, flag_arc_is_flipped, debug_fig_num)

% Call the function to convert from ST back to XY
st_parameters_type_strings{1} = 'arc';
st_parameters_type_strings{2} = 'line';
st_parameters_type_strings{3} = revised_intermediate_geometry_join_type;
st_parameters{1} = revised_arc_parameters_St;
st_parameters{2} = revised_line_parameters_St;
st_parameters{3} = revised_intermediate_geometry_join_parameters_St;

[XY_parameters] = ...
fcn_geometry_orientGeometrySt2XY(st_parameters_type_strings, st_parameters, St_transform_XYtoSt, flag_arc_is_flipped, (-1));

revised_arc_parameters = XY_parameters{1};
revised_line_parameters = XY_parameters{2};
revised_intermediate_geometry_join_parameters = XY_parameters{3};

if ~isempty(debug_fig_num)
    % Plot the results
    figure(debug_fig_num);
    subplot(3,2,1);
    debug_axis = axis;

    subplot(3,2,6);

    fcn_geometry_plotGeometry('arc',revised_arc_parameters);
    fcn_geometry_plotGeometry('line',revised_line_parameters);
    fcn_geometry_plotGeometry(revised_intermediate_geometry_join_type,revised_intermediate_geometry_join_parameters);
    
    title('Final outputs');
    axis(debug_axis);
end
end % Ends fcn_INTERNAL_convertParametersOutOfStOrientation
