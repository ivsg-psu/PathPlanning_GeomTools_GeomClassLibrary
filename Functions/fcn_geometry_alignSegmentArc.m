function [revised_segment_parameters, revised_arc_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignSegmentArc(segment_parameters, arc_parameters, varargin)
%% fcn_geometry_alignSegmentArc
% Revises the geometric parameters of an arc and line segment such that
% they align where they join. It does this by checking the offset between
% the two objects at the join location. 
% 
% If the alignment is not feasible but the offset is less than a threshold
% (default is 0.1 meter), the arc's geometric position is shifted to force
% alignment with the line. In other words, the line is kept stationary and
% the arc is aligned to the line segment.
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
% [revised_segment_parameters, revised_arc_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters]  = ...
% fcn_geometry_alignSegmentArc(segment_parameters, arc_parameters, (threshold), (continuity_level),  (fig_num))
%
% INPUTS:
%
%      segment_parameters: the parameter set describing the line segment
%      geometry. See fcn_geometry_fillEmptyDomainStructure for details,
%      specifically the structure for 'Vector regression segment fit'.
%
%      arc_parameters: the parameter set describing the arc geometry. See
%      fcn_geometry_fillEmptyDomainStructure for details, specifically the
%      structure for 'Regression arc'.
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
%      revised_segment_parameters: the parameter set describing the line
%      segment geometry that joins the geometries. See
%      fcn_geometry_fillEmptyDomainStructure for details, specifically the
%      structure for 'Vector regression segment fit'.
%
%      revised_arc_parameters: the parameter set describing the arc
%      segment geometry that joins the geometries. See
%      fcn_geometry_fillEmptyDomainStructure for details, specifically the
%      structure for 'Regression arc'.
%
%      revised_intermediate_geometry_join_type: for type C2 continuity, an
%      intermediate geometry is often inserted in the form of a spiral.
%      This output saves the geometry type as a string type.
%
%      revised_intermediate_geometry_join_parameters: the parameter set describing the
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
% See the script: script_test_fcn_geometry_alignSegmentArc
% for a full test suite.
%
% This function was written on 2024_05_13 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_05_13 - Sean Brennan
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
    debug_fig_num = 7564;
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
if (3<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        threshold = temp;
    else 
        threshold = [];
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
fcn_INTERNAL_prepDebugFigure(arc_parameters, segment_parameters, debug_fig_num);

%% Check to see if arc and segment intersect
intersection_point1 = fcn_INTERNAL_ArcSegmentIntersection(arc_parameters, segment_parameters, 1, debug_fig_num);
if (length(intersection_point1(:,1))>1)&&(0==continuity_level)
    warning('Multiple intersection points found between a segment and arc geometry with a requested connection type of C0 continuity. Unable to resolve which intersection to use.');
end

%% Flip the arc's ordering
arc_parameters_flipped = fcn_geometry_flipGeom('arc',  arc_parameters, -1);
segment_parameters_flipped = fcn_geometry_flipGeom('segment',  segment_parameters, -1);


if ~isempty(debug_fig_num)
    % Plot the results
    figure(debug_fig_num);
    subplot(3,2,1);
    debug_axis = axis;

    subplot(3,2,2);

    fcn_geometry_plotGeometry('arc',arc_parameters_flipped);
    fcn_geometry_plotGeometry('segment',segment_parameters_flipped);
        
    title('Flipped ordering');
    axis(debug_axis);
end

%% Call alignArcSegment to get inverse solution

[revised_inverse_arc_parameters, revised_inverse_segment_parameters, revised_inverse_intermediate_geometry_join_type, revised_inverse_intermediate_geometry_join_parameters] = ...
    fcn_geometry_alignArcSegment(...
    arc_parameters_flipped, segment_parameters, (threshold), (continuity_level), (2345));

if ~isempty(debug_fig_num)
    % Plot the results
    figure(debug_fig_num);
    subplot(3,2,1);
    debug_axis = axis;

    subplot(3,2,3);

    fcn_geometry_plotGeometry('arc',revised_inverse_arc_parameters);
    fcn_geometry_plotGeometry('segment',revised_inverse_segment_parameters);
    fcn_geometry_plotGeometry(revised_inverse_intermediate_geometry_join_type,revised_inverse_intermediate_geometry_join_parameters);
    
    title('Inverse ArcSegment solution');
    axis(debug_axis);
end

%% Flip the direction of the geometries

% Flip the arc ordering
revised_offset_arc_parameters = fcn_geometry_flipGeom('arc',  revised_inverse_arc_parameters, -1);

% Flip the segment ordering
revised_offset_segment_parameters = fcn_geometry_flipGeom('segment',  revised_inverse_segment_parameters, -1);

% Flip the intermediate geometry
revised_offset_intermediate_geometry_join_parameters = fcn_geometry_flipGeom(revised_inverse_intermediate_geometry_join_type,  revised_inverse_intermediate_geometry_join_parameters, -1);

if ~isempty(debug_fig_num)
    % Plot the results
    figure(debug_fig_num);
    subplot(3,2,1);
    debug_axis = axis;

    subplot(3,2,4);

    fcn_geometry_plotGeometry('arc',revised_offset_arc_parameters);
    fcn_geometry_plotGeometry('segment',revised_offset_segment_parameters);
    fcn_geometry_plotGeometry(revised_inverse_intermediate_geometry_join_type,revised_offset_intermediate_geometry_join_parameters);
    
    title('Solution offset but flipped');
    axis(debug_axis);
end


%% Find transform that moves the line segment to desired solution
% Renormalize the segment result
revised_offset_segment_parameters = fcn_INTERNAL_renormalizeSegment(revised_offset_segment_parameters);

% Find how much the arc moved
segment_old_base_point = segment_parameters(1,3:4);
segment_new_base_point = revised_offset_segment_parameters(1,3:4);
translation_vector = segment_new_base_point - segment_old_base_point;

% create a transformation that moves one to the other
transformMatrix_translation_into_origin = se2(0,'theta',translation_vector);
St_transform_XYtoSt = transformMatrix_translation_into_origin;

% Apply the inverse of the transform to the results

[revised_arc_parameters, revised_segment_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_INTERNAL_convertParametersOutOfStOrientation(...
    revised_offset_arc_parameters, revised_offset_segment_parameters, revised_inverse_intermediate_geometry_join_type, revised_offset_intermediate_geometry_join_parameters, St_transform_XYtoSt, 0, debug_fig_num);


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
    fcn_geometry_plotGeometry('segment',segment_parameters);
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
    fcn_geometry_plotGeometry('segment',revised_segment_parameters);
    fcn_geometry_plotGeometry('arc',revised_arc_parameters);
    fcn_geometry_plotGeometry('spiral',revised_intermediate_geometry_join_parameters);


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
    fcn_geometry_plotGeometry('segment',line_parameters);

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

end
end % Ends fcn_INTERNAL_prepDebugFigure

%% fcn_INTERNAL_ArcSegmentIntersection
function  intersection_point_arc_to_segment = fcn_INTERNAL_ArcSegmentIntersection(arc_parameters,segment_parameters, subplot_number, debug_fig_num)
firstFitType = 'arc';
firstFitType_parameters = arc_parameters;
secondFitType = 'segment';
secondFitType_parameters = segment_parameters;

intersection_point_arc_to_segment = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, -1);


if ~isempty(debug_fig_num)
    % Plot the intersection
    figure(debug_fig_num);
    subplot(3,2,1);
    debug_axis = axis;

    subplot(3,2,subplot_number);
    plot(intersection_point_arc_to_segment(:,1),intersection_point_arc_to_segment(:,2),'co','MarkerSize',10,'LineWidth',2);

    axis(debug_axis);
end

end % Ends fcn_INTERNAL_ArcSegmentIntersection


%% fcn_INTERNAL_findShiftToMatchSegmentToArc
function [desired_arc_parameters, desired_segment_parameters, desired_intermediate_geometry_join_parameters, desired_intermediate_geometry_join_type] = ...
    fcn_INTERNAL_findShiftToMatchSegmentToArc(arc_parameters, segment_parameters, continuity_level, intersection_point, threshold, flag_perform_shift_of_segment, debug_fig_num)
% Calculates the delta amount to match the segment to the arc. The delta
% values are measured FROM desired point TO actual point

if length(threshold)==1
    transverse_threshold = threshold;
else
    transverse_threshold = threshold(2);
end

% Calculate needed values from parameter sets
% Calculate needed values from parameter sets
% Get the arc fit details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
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
segment_unit_tangent_vector     = fcn_geometry_calcUnitVector(segment_parameters(1,1:2));
segment_base_point_xy           = segment_parameters(1,3:4);
segment_s_start                 = segment_parameters(1,5);
segment_s_end                   = segment_parameters(1,6);
% segment_angle = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_unit_ortho_vector = segment_unit_tangent_vector*[0 1; -1 0];
segment_length                  = segment_s_end - segment_s_start;

% Calculate the distance from the arc to the segment
vector_from_arc_center_to_segment_joint = segment_base_point_xy - arc_center_xy;
if arc_is_counter_clockwise
    distance_from_arc_center_to_segment = -sum(segment_unit_ortho_vector.*vector_from_arc_center_to_segment_joint,2);
else
    distance_from_arc_center_to_segment = sum(segment_unit_ortho_vector.*vector_from_arc_center_to_segment_joint,2);
end


% Calculate the distance between the circles and the join point
space_between_arc_and_segment = distance_from_arc_center_to_segment - arc_radius;


% With the cleaned parameters, the line vector always
% points in the direction of the join point.
% to_joint_unit_tangent_vector = [1 0];
% to_joint_unit_ortho_vector   = [0 1];

switch continuity_level
    case 0
        % For C0 continuity of arc to segment, the closest point of the
        % desired joint to the arc is simply the end of arc

        desired_intermediate_geometry_join_type       = ''; % Intermediate geometry will be a line segement
        desired_intermediate_geometry_join_parameters = nan(1,6);

        if ~any(isnan(intersection_point))
            % Arc will end at the intersection, which is at the origin
            vector_from_arc_center_to_segment_joint = [0 0] - arc_center_xy;
            angle_of_intersection = atan2(vector_from_arc_center_to_segment_joint(2),vector_from_arc_center_to_segment_joint(1));

            desired_arc_parameters        = arc_parameters;
            desired_arc_parameters(1,5)   = angle_of_intersection; % Update where the spiral ends

            desired_segment_parameters        = segment_parameters;
            desired_segment_parameters(1,3:4) = [0 0]; % Move segment_base_point_xy

        else
            if flag_perform_shift_of_segment==1
                if abs(space_between_arc_and_segment)<=transverse_threshold
                    % Yes, segment can be moved enough to be tangent. So put
                    % segment's start at the end of the arc. Do not change the
                    % arc.
                    desired_arc_parameters        = arc_parameters;

                    desired_segment_parameters        = segment_parameters;
                    desired_segment_parameters(1,3:4) = [0 0]; % Move segment_base_point_xy

                else
                    % Not possible to shift enough to allow connection
                    desired_arc_parameters        = nan(size(arc_parameters));
                    desired_segment_parameters    = nan(1,6);
                end
            else
                % No shift allowed by user entry, so not possible
                desired_arc_parameters            = nan(size(arc_parameters));
                desired_segment_parameters        = nan(1,6);
            end
        end

    case 1
        % For C1 continuity of arc to segment, the closest point of the
        % desired joint to the arc is always going to be the point where
        % the arc touches the tangent line of the segment. Since the
        % segment is the only one that can be moved, the connecting point
        % for the arc is simply the location where the arc is tangent, which by
        % construction is [0 0]. The connection point for the segment is at
        % this same point

        desired_intermediate_geometry_join_type       = 'segment'; % Intermediate geometry will be a line segement
        desired_intermediate_geometry_join_parameters = nan(1,6);

        if flag_perform_shift_of_segment==1
            if abs(space_between_arc_and_segment)<=transverse_threshold
                % Yes, arc2 can be moved enough to be tangent. So put
                % arc2's center in correct place
                desired_arc_parameters        = arc_parameters;
                desired_arc_parameters(1,5)   = -pi/2;

                desired_segment_parameters        = segment_parameters;
                desired_segment_parameters(1,3:4) = [0 0]; % Move segment_base_point_xy

            else
                % Not possible to shift enough to allow connection
                desired_arc_parameters        = nan(size(arc_parameters));
                desired_segment_parameters    = nan(1,6);
            end
        else
            % No shift allowed by user entry, so not possible
            desired_arc_parameters            = nan(size(arc_parameters));
            desired_segment_parameters        = nan(1,6);
        end


    case 2
        
        % For C2 continuity of a line to an arc, the arc must connect to the
        % line via a spiral. Calculation of the spiral is difficult and
        % requires numerical iteration, and the inputs require the calculation
        % of the offset of the outer part of the circle from the line. This
        % offset in the direction of the circle center MUST be strictly
        % positive or the spiral cannot work.
        
        desired_intermediate_geometry_join_type       = 'spiral'; % Intermediate geometry will be a line segement
        desired_intermediate_geometry_join_parameters = nan(1,6);

        flag_spiral_was_calculated = 0;
        if space_between_arc_and_segment>0
            % spiral_join_parameters = [spiralLength,h0,x0,y0,K0,Kf];

            [desired_intermediate_geometry_join_parameters, ~] = ...
                fcn_geometry_spiralFromCircleToCircle(arc_radius, inf, [0 -space_between_arc_and_segment], [], -1);
            flag_spiral_was_calculated = 1;
        else
            % Not enough space between arc and segment for a spiral, not
            % without modifying it
            if flag_perform_shift_of_segment==1
                if abs(space_between_arc_and_segment)<=transverse_threshold
                    unit_direction_vector_to_shift_segment = [0 -1];

                    revised_segment_parameters = segment_parameters;
                    revised_segment_parameters(1,3:4) = segment_parameters(1,3:4) + unit_direction_vector_to_shift_segment*(abs(space_between_arc_and_segment)+0.001);

                    % Call function again with revised parameters that
                    % should work
                    [desired_arc_parameters, desired_segment_parameters, desired_intermediate_geometry_join_parameters] = ...
                        fcn_INTERNAL_findShiftToMatchSegmentToArc(arc_parameters, revised_segment_parameters, continuity_level, intersection_point, [], 0, debug_fig_num);
                    flag_spiral_was_calculated = 1;
                else
                    % Not possible
                    desired_arc_parameters = nan(size(arc_parameters));
                    desired_segment_parameters = nan(1,6);
                end
            else
                % Not allowed by user entry
                desired_arc_parameters = nan(size(arc_parameters));
                desired_segment_parameters = nan(1,6);
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
                fcn_geometry_plotGeometry('segment',segment_parameters);
                fcn_geometry_plotGeometry(desired_intermediate_geometry_join_type,desired_intermediate_geometry_join_parameters);

            end % Ends plotting

            % Find the angle and position that the spiral starts and ends at
            % spiralLength = desired_intermediate_geometry_join_parameters(1,1);
            h0           = desired_intermediate_geometry_join_parameters(1,2);
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
                desired_segment_parameters = nan(1,6);
                flag_spiral_is_bad = 1;
            end

            % Make sure the spiral's end position is within the distance range of the segement. If not, return Nan values for everything
            % Find the x-value of the spiral result
            XY_spiral = fcn_geometry_plotGeometry(desired_intermediate_geometry_join_type,desired_intermediate_geometry_join_parameters,[],[],-1);
            spiral_end_XY = XY_spiral(end,:);
            if spiral_end_XY(1,1)<0 || spiral_end_XY(1,1)>segment_length
                % Spiral does not connect with the segment.
                desired_arc_parameters = nan(size(arc_parameters));
                desired_segment_parameters = nan(1,6);
                flag_spiral_is_bad = 1;
            end

            if 0==flag_spiral_is_bad
                % If made it here, then the spiral is good. Fill outputs and then
                % return results
                desired_arc_parameters = arc_parameters;
                desired_arc_parameters(1,5)   = arc_angle_where_spiral_starts; % Update where the spiral ends

                desired_segment_parameters        = segment_parameters;
                desired_segment_parameters(1,3:4) = spiral_end_XY; % segment_base_point_xy
                desired_segment_parameters(1,5)   = 0; % segment_s_start
                desired_segment_parameters(1,6)   = segment_length - spiral_end_XY(1,1); % segment_s_end
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
    fcn_geometry_plotGeometry('segment',desired_segment_parameters);
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



%% fcn_INTERNAL_convertParametersOutOfStOrientation
function [revised_arc_parameters, revised_segment_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_INTERNAL_convertParametersOutOfStOrientation(...
    revised_arc_parameters_St, revised_segment_parameters_St, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters_St, St_transform_XYtoSt, flag_arc_is_flipped, debug_fig_num)

% Call the function to convert from ST back to XY
st_parameters_type_strings{1} = 'segment';
st_parameters_type_strings{2} = 'arc';
st_parameters_type_strings{3} = revised_intermediate_geometry_join_type;
st_parameters{1} = revised_segment_parameters_St;
st_parameters{2} = revised_arc_parameters_St;
st_parameters{3} = revised_intermediate_geometry_join_parameters_St;

[XY_parameters] = ...
fcn_geometry_orientGeometrySt2XY(st_parameters_type_strings, st_parameters, St_transform_XYtoSt, flag_arc_is_flipped, (-1));

revised_segment_parameters = XY_parameters{1};
revised_arc_parameters     = XY_parameters{2};
revised_intermediate_geometry_join_parameters = XY_parameters{3};

if ~isempty(debug_fig_num)
    % Plot the results
    figure(debug_fig_num);
    subplot(3,2,1);
    debug_axis = axis;

    subplot(3,2,6);

    fcn_geometry_plotGeometry('arc',revised_arc_parameters);
    fcn_geometry_plotGeometry('segment',revised_segment_parameters);
    fcn_geometry_plotGeometry(revised_intermediate_geometry_join_type,revised_intermediate_geometry_join_parameters);
    
    title('Final outputs');
    axis(debug_axis);
end
end % Ends fcn_INTERNAL_convertParametersOutOfStOrientation

%% fcn_INTERNAL_renormalizeSegment
function new_segment_parameters = fcn_INTERNAL_renormalizeSegment(old_segment_parameters)
segment_unit_tangent_vector = old_segment_parameters(1,1:2);
segment_base_point_xy       = old_segment_parameters(1,3:4);
segment_s_start             = old_segment_parameters(1,5);
segment_s_end               = old_segment_parameters(1,6);
% Make sure the line segment is well-formed, e.g. the station at the end is
% larger than the station at the start. If not, need to correct
if segment_s_end<segment_s_start
    % Flip the order
    segment_s_start         = segment_parameters(1,6);
    segment_s_end           = segment_parameters(1,5);

    % Flip the vector
    segment_unit_tangent_vector = -segment_unit_tangent_vector;

end
segment_start_xy            = segment_base_point_xy + segment_unit_tangent_vector*segment_s_start;
segment_end_xy              = segment_base_point_xy + segment_unit_tangent_vector*segment_s_end;
segment_length = sum((segment_end_xy - segment_start_xy).^2,2).^0.5;

new_segment_parameters(1,1:2) = segment_unit_tangent_vector;
new_segment_parameters(1,3:4) = segment_start_xy;
new_segment_parameters(1,5)   = 0;
new_segment_parameters(1,6)   = segment_length;
end % Ends fcn_INTERNAL_renormalizeSegment