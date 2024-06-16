function [revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = fcn_geometry_alignArcArc(arc1_parameters, arc2_parameters, varargin)
%% fcn_geometry_alignArcArc
% Revises the geometric parameters of an arc connecting to another arc such
% that they align where they join. It does this by checking the offset
% between the two objects at the join location. If the offset is less than
% a threshold (default is 0.1 meter), the second geometric object is
% shifted to force alignment or, for increased continuity, an intermediate
% geometric object such as a spiral is inserted between the two arcs.
%
% Format:
% [revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters]  = ...
% fcn_geometry_alignArcArc(arc1_parameters, arc2_parameters, (threshold), (continuity_level),  (fig_num))
%
% INPUTS:
%
%      arc1_parameters: the parameter set describing the 1st arc geometry. See
%      fcn_geometry_fillEmptyDomainStructure for details, specifically the
%      structure for 'Regression arc'.
%
%      arc2_parameters: the parameter set describing the 2nd arc geometry. See
%      fcn_geometry_fillEmptyDomainStructure for details, specifically the
%      structure for 'Regression arc'.
%
%      (OPTIONAL INPUTS)
%
%      threshold: the offset, in meters, between the arc2 and arc2 such
%      that this offset is removed by shifting the geometry of arc2. If the
%      offset is larger than the threshold, then the outputs are set to
%      empty. If threshold is entered as a 2x1 or 1x2, then this specifies
%      the threshold first in the transverse direction, and then in the
%      station direction. For example, an entry of [0.02 3] would have 0.02
%      meters threshold in the transverse direction, but 3 meters threshold
%      in the station direction. To force a joining of geometries without
%      using shifts, leave threshold empty, e.g. threshold = [];
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
%      revised_arc1_parameters: the parameter set describing the 1st arc
%      segment geometry that joins the geometries. See
%      fcn_geometry_fillEmptyDomainStructure for details, specifically the
%      structure for 'Regression arc'.
%
%      revised_arc2_parameters: the parameter set describing the 2nd arc
%      segment geometry that joins the geometries. See
%      fcn_geometry_fillEmptyDomainStructure for details, specifically the
%      structure for 'Regression arc'.
%
%      revised_intermediate_geometry_join_type: for type C1 and C2
%      continuity, an intermediate geometry is often inserted in the form
%      of a line segment or spiral, respectively. This output saves the
%      geometry type as a string type.
%
%      revised_intermediate_geometry_join_parameters: the parameter set describing the
%      spiral segment geometry that joins the arc1 and arc2 geometries if C2
%      continuity is specified. See fcn_geometry_fillEmptyDomainStructure
%      for details, specifically the structure for 'spiral'.
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_alignArcArc
% for a full test suite.
%
% This function was written on 2024_04_21 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_04_21 - S. Brennan
% -- wrote the code
% 2024_05_09 - S. Brennan
% -- added intermediate type outputs to allow line segments or spirals
% 2024_05_28 - S. Brennan
% -- fixed call to spiralFromCircleToCircle to use parameter vectors
% -- fixed bug related to angle checks in arcs when oriented opposite
% 2024_06_16 - Sean Brennan
% -- changed parameter format to new style:
%            'spiral' - 
%
%               [
%                x0,  % The initial x value
%                y0,  % The initial y value
%                h0,  % The initial heading
%                s_Length,  % the s-coordinate length allowed
%                K0,  % The initial curvature
%                Kf   % The final curvature
%              ] 

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
    debug_fig_num = 999978;
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
        narginchk(2,5);

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
flag_perform_shift_of_arc2 = 1;
if (3<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        threshold = temp;
    else
        flag_perform_shift_of_arc2 = 0;
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
fcn_INTERNAL_prepDebugFigure(arc1_parameters, arc2_parameters, debug_fig_num);

%% Check to see if arc1 and arc2 intersect
intersection_point1 = fcn_INTERNAL_ArcArcIntersection(arc1_parameters, arc2_parameters, 1, debug_fig_num);

%% Rearrange parameters so line is always the 1st input, arc is 2nd
% Fix the parameters to make the ordering is correct
[clean_arc1_parameters, clean_arc2_parameters] = fcn_INTERNAL_fixOrientationAndOrdering(arc1_parameters, arc2_parameters, intersection_point1, debug_fig_num);

%% Get new intersection point, if arcs changed shape
intersection_point2 = fcn_INTERNAL_ArcArcIntersection(clean_arc1_parameters,clean_arc2_parameters, 2, debug_fig_num);

%% Rotate the geometries out of XY into ST coordinates
% so that the tangent line is oriented horizontally
% and the start of the tangent line on arc1 is at the origin.
% This is to make the debugging MUCH easier, as it reduces permutations.
% Again, this is fixed in later steps.
[st_arc1_parameters, st_arc2_parameters, St_transform_XYtoSt, flag_arc1_is_flipped] = ...
    fcn_INTERNAL_convertParametersToStOrientation(clean_arc1_parameters, clean_arc2_parameters, continuity_level, intersection_point2, debug_fig_num); 

%% Check how much shift is needed to connect arc1 to arc2
[desired_st_arc1_parameters, desired_st_arc2_parameters, desired_st_intermediate_geometry_join_parameters, desired_intermediate_geometry_join_type] = ...
    fcn_INTERNAL_findShiftToMatchArcToArc(st_arc1_parameters, st_arc2_parameters, continuity_level, intersection_point2, threshold, flag_perform_shift_of_arc2, debug_fig_num);
% Deltas are from desired to actual

%% Perform shift to join arcs 1 and 2
[revised_arc1_parameters_St,revised_arc2_parameters_St, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters_St] = ...
    fcn_INTERNAL_performShift(threshold, continuity_level, ...
    st_arc1_parameters, st_arc2_parameters, ...
    desired_st_arc1_parameters, desired_st_arc2_parameters, ...
    desired_st_intermediate_geometry_join_parameters, desired_intermediate_geometry_join_type, debug_fig_num);

%% Rotate results out of St back into XY
[revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_INTERNAL_convertParametersOutOfStOrientation(...
    revised_arc1_parameters_St, revised_arc2_parameters_St, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters_St, St_transform_XYtoSt, flag_arc1_is_flipped, debug_fig_num);

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
    sgtitle('Arc joining with arc');
    title('Original');
    xlabel('X [meters]');
    ylabel('Y [meters]')


    % Plot the inputs
    fcn_geometry_plotCircle(arc1_parameters(1,1:2),arc1_parameters(1,3),...
        sprintf(' ''--'',''Color'',[0 0.6 0],''LineWidth'',1 '),fig_num);
    fcn_geometry_plotCircle(arc2_parameters(1,1:2),arc2_parameters(1,3),...
        sprintf(' ''--'',''Color'',[0.6 0 0],''LineWidth'',1 '),fig_num);
    
    fcn_geometry_plotGeometry('arc',arc1_parameters);
    fcn_geometry_plotGeometry('arc',arc2_parameters);

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
    fcn_geometry_plotCircle(revised_arc1_parameters(1,1:2),revised_arc1_parameters(1,3),...
        sprintf(' ''--'',''Color'',[0 0.6 0],''LineWidth'',1 '),fig_num);
    fcn_geometry_plotCircle(revised_arc2_parameters(1,1:2),revised_arc2_parameters(1,3),...
        sprintf(' ''--'',''Color'',[0.6 0 0],''LineWidth'',1 '),fig_num);

    fcn_geometry_plotGeometry('arc',revised_arc1_parameters);
    fcn_geometry_plotGeometry('arc',revised_arc2_parameters);
     fcn_geometry_plotGeometry(revised_intermediate_geometry_join_type,revised_intermediate_geometry_join_parameters);


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
function [clean_arc1_parameters, clean_arc2_parameters] = fcn_INTERNAL_fixOrientationAndOrdering(arc1_parameters, arc2_parameters, intersection_point, debug_fig_num)
% This function takes the parameter inputs and produces parameter sets such
% that the arc1 is first, it is oriented so that it ends at the junction
% with arc2, and the arc2 starts at the junction. 

% Get the arc fit details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_center_xy                = arc1_parameters(1,1:2);
arc1_radius                   = arc1_parameters(1,3);
arc1_start_angle_in_radians   = arc1_parameters(1,4);
arc1_end_angle_in_radians     = arc1_parameters(1,5);
arc1_start_xy                 = arc1_center_xy + arc1_radius*[cos(arc1_start_angle_in_radians) sin(arc1_start_angle_in_radians)];
arc1_end_xy                   = arc1_center_xy + arc1_radius*[cos(arc1_end_angle_in_radians) sin(arc1_end_angle_in_radians)];
arc1_is_counter_clockwise     = arc1_parameters(1,7);

% Find the change in angle of the arc
arc1_start_unit_vector        = [cos(arc1_start_angle_in_radians) sin(arc1_start_angle_in_radians)];
arc1_end_unit_vector          = [cos(arc1_end_angle_in_radians)   sin(arc1_end_angle_in_radians)  ];
if arc1_is_counter_clockwise
    cross_product_direction = 1;
else
    cross_product_direction = -1;
end
arc1_change_in_angle = fcn_geometry_findAngleUsing2PointsOnCircle([0 0],1, arc1_start_unit_vector, arc1_end_unit_vector, cross_product_direction);


% Get the arc fit details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_center_xy                = arc2_parameters(1,1:2);
arc2_radius                   = arc2_parameters(1,3);
arc2_start_angle_in_radians   = arc2_parameters(1,4);
arc2_end_angle_in_radians     = arc2_parameters(1,5);
arc2_start_xy                 = arc2_center_xy + arc2_radius*[cos(arc2_start_angle_in_radians) sin(arc2_start_angle_in_radians)];
arc2_end_xy                   = arc2_center_xy + arc2_radius*[cos(arc2_end_angle_in_radians) sin(arc2_end_angle_in_radians)];
arc2_is_counter_clockwise     = arc2_parameters(1,7);

% Find the change in angle of the arc
arc2_start_unit_vector        = [cos(arc2_start_angle_in_radians) sin(arc2_start_angle_in_radians)];
arc2_end_unit_vector          = [cos(arc2_end_angle_in_radians)   sin(arc2_end_angle_in_radians)  ];
if arc2_is_counter_clockwise
    cross_product_direction = 1;
else
    cross_product_direction = -1;
end
arc2_change_in_angle = fcn_geometry_findAngleUsing2PointsOnCircle([0 0],1, arc2_start_unit_vector, arc2_end_unit_vector, cross_product_direction);


% Find arc1 and arc2 join points, e.g. where they meet. This can
% happen at either end
if ~any(isnan(intersection_point))
    arc1_endPoint = intersection_point;
else
    arc1_endPoint = arc1_end_xy;
end

distances_to_check = sum((...
    [arc1_start_xy; arc1_start_xy; arc1_endPoint; arc1_endPoint] - [arc2_start_xy; arc2_end_xy; arc2_start_xy; arc2_end_xy]).^2,2).^0.5;
[~,closest_pair] = min(distances_to_check);

% Fix the arc or line depending on which combo is closest
switch closest_pair
    case 1 % arc1 start, arc2 start
        % The arc1 is entering the junction at its start. This is
        % not correct. Need to "flip" the arc1's orientation.
        if 1==arc1_is_counter_clockwise
            corrected_arc1_is_counter_clockwise = 0;
        else
            corrected_arc1_is_counter_clockwise = 1;
        end
        corrected_arc1_start_angle_in_radians = atan2(arc1_end_unit_vector(2),arc1_end_unit_vector(1));
        corrected_arc1_end_angle_in_radians   = corrected_arc1_start_angle_in_radians - arc1_change_in_angle;

        % The arc2 is leaving the junction at its start. This is
        % correct so just pass through the variables
        corrected_arc2_is_counter_clockwise   = arc2_is_counter_clockwise;
        corrected_arc2_start_angle_in_radians = atan2(arc2_start_unit_vector(2),arc2_start_unit_vector(1));
        corrected_arc2_end_angle_in_radians   = corrected_arc2_start_angle_in_radians + arc2_change_in_angle;

    case 2 % arc1 start, arc2 end
        % The arc1 is entering the junction at its start. This is
        % not correct. Need to "flip" the arc1's orientation.
        if 1==arc1_is_counter_clockwise
            corrected_arc1_is_counter_clockwise = 0;
        else
            corrected_arc1_is_counter_clockwise = 1;
        end
        corrected_arc1_start_angle_in_radians = atan2(arc1_end_unit_vector(2),arc1_end_unit_vector(1));
        corrected_arc1_end_angle_in_radians   = corrected_arc1_start_angle_in_radians - arc1_change_in_angle;

        % The arc2 is entering the junction at its end. This is
        % not correct. Need to "flip" the arc2's orientation.
        if 1==arc2_is_counter_clockwise
            corrected_arc2_is_counter_clockwise = 0;
        else
            corrected_arc2_is_counter_clockwise = 1;
        end
        corrected_arc2_start_angle_in_radians = atan2(arc2_end_unit_vector(2),arc2_end_unit_vector(1));
        corrected_arc2_end_angle_in_radians   = corrected_arc2_start_angle_in_radians - arc2_change_in_angle;

    case 3 % arc1 end, arc2 start
        % The arc1 is entering the junction at its end. This is
        % correct so just pass through the variables
        corrected_arc1_is_counter_clockwise   = arc1_is_counter_clockwise;
        corrected_arc1_start_angle_in_radians = atan2(arc1_start_unit_vector(2),arc1_start_unit_vector(1));
        corrected_arc1_end_angle_in_radians   = corrected_arc1_start_angle_in_radians + arc1_change_in_angle;

        % The arc2 is leaving the junction at its start. This is
        % correct so just pass through the variables
        corrected_arc2_is_counter_clockwise = arc2_is_counter_clockwise;
        corrected_arc2_start_angle_in_radians = atan2(arc2_start_unit_vector(2),arc2_start_unit_vector(1));
        corrected_arc2_end_angle_in_radians   = corrected_arc2_start_angle_in_radians + arc2_change_in_angle;

    case 4 % arc1 end, arc2 end
        % The arc1 is entering the junction at its end. This is
        % correct so just pass through the variables
        corrected_arc1_is_counter_clockwise   = arc1_is_counter_clockwise;
        corrected_arc1_start_angle_in_radians = atan2(arc1_start_unit_vector(2),arc1_start_unit_vector(1));
        corrected_arc1_end_angle_in_radians   = corrected_arc1_start_angle_in_radians + arc1_change_in_angle;

        % The arc2 is entering the junction at its end. This is
        % not correct. Need to "flip" the arc2's orientation.
        if 1==arc2_is_counter_clockwise
            corrected_arc2_is_counter_clockwise = 0;
        else
            corrected_arc2_is_counter_clockwise = 1;
        end
        corrected_arc2_start_angle_in_radians = atan2(arc2_end_unit_vector(2),arc2_end_unit_vector(1));
        corrected_arc2_end_angle_in_radians   = corrected_arc2_start_angle_in_radians - arc2_change_in_angle;

    otherwise
        error('Impossible case encountered - must stop!');
end

% Set the outputs

clean_arc1_parameters(1,1:2)   = arc1_parameters(1,1:2); % center of the arc does not change
clean_arc1_parameters(1,3)     = arc1_parameters(1,3);   % radius of the arc does not change
clean_arc1_parameters(1,4:5)   = [corrected_arc1_start_angle_in_radians corrected_arc1_end_angle_in_radians];
clean_arc1_parameters(1,6)     = arc1_parameters(1,6);   % flag is circle
clean_arc1_parameters(1,7)     = corrected_arc1_is_counter_clockwise;


clean_arc2_parameters(1,1:2)   = arc2_parameters(1,1:2); % center of the arc does not change
clean_arc2_parameters(1,3)     = arc2_parameters(1,3);   % radius of the arc does not change
clean_arc2_parameters(1,4:5)   = [corrected_arc2_start_angle_in_radians corrected_arc2_end_angle_in_radians];
clean_arc2_parameters(1,6)     = arc2_parameters(1,6);   % flag is circle
clean_arc2_parameters(1,7)     = corrected_arc2_is_counter_clockwise;

if ~isempty(debug_fig_num)
    figure(debug_fig_num);
    subplot(3,2,1);
    debug_axis = axis;

    % Plot the cleaned inputs
    subplot(3,2,2);

    fcn_geometry_plotCircle(clean_arc1_parameters(1,1:2),clean_arc1_parameters(1,3),...
        sprintf(' ''--'',''Color'',[0 0.6 0],''LineWidth'',1 '),debug_fig_num);
    fcn_geometry_plotCircle(clean_arc2_parameters(1,1:2),clean_arc2_parameters(1,3),...
        sprintf(' ''--'',''Color'',[0.6 0 0],''LineWidth'',1 '),debug_fig_num);

    fcn_geometry_plotGeometry('arc',clean_arc1_parameters);
    fcn_geometry_plotGeometry('arc',clean_arc2_parameters);

    axis(debug_axis);
    hold on;
    grid on;
    axis equal;
    xlabel('X [meters]');
    ylabel('Y [meters]')

    title('Corrected inputs');
end


end % ends fcn_INTERNAL_fixOrientationAndOrdering


%% fcn_INTERNAL_ArcArcIntersection
function  intersection_point_both_arcs = fcn_INTERNAL_ArcArcIntersection(arc1_parameters,arc2_parameters, subplot_number, debug_fig_num)

% Calculate needed values from parameter sets
% Get the arc fit details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_center_xy                = arc1_parameters(1,1:2);
arc1_radius                   = arc1_parameters(1,3);
arc1_start_angle_in_radians   = arc1_parameters(1,4);
arc1_end_angle_in_radians     = arc1_parameters(1,5);
arc1_start_xy                 = arc1_center_xy + arc1_radius*[cos(arc1_start_angle_in_radians) sin(arc1_start_angle_in_radians)];
arc1_end_xy                   = arc1_center_xy + arc1_radius*[cos(arc1_end_angle_in_radians) sin(arc1_end_angle_in_radians)];
arc1_is_counter_clockwise     = arc1_parameters(1,7);

% Find the change in angle of the arc
% arc1_start_unit_vector        = [cos(arc1_start_angle_in_radians) sin(arc1_start_angle_in_radians)];
% arc1_end_unit_vector          = [cos(arc1_end_angle_in_radians)   sin(arc1_end_angle_in_radians)  ];
% if arc1_is_counter_clockwise
%     cross_product_direction = 1;
% else
%     cross_product_direction = -1;
% end
% arc1_change_in_angle = fcn_geometry_findAngleUsing2PointsOnCircle([0 0],1, arc1_start_unit_vector, arc1_end_unit_vector, cross_product_direction);


% Get the arc fit details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_center_xy                = arc2_parameters(1,1:2);
arc2_radius                   = arc2_parameters(1,3);
arc2_start_angle_in_radians   = arc2_parameters(1,4);
arc2_end_angle_in_radians     = arc2_parameters(1,5);
arc2_start_xy                 = arc2_center_xy + arc2_radius*[cos(arc2_start_angle_in_radians) sin(arc2_start_angle_in_radians)];
arc2_end_xy                   = arc2_center_xy + arc2_radius*[cos(arc2_end_angle_in_radians) sin(arc2_end_angle_in_radians)];
arc2_is_counter_clockwise     = arc2_parameters(1,7);

% Find the change in angle of the arc
% arc2_start_unit_vector        = [cos(arc2_start_angle_in_radians) sin(arc2_start_angle_in_radians)];
% arc2_end_unit_vector          = [cos(arc2_end_angle_in_radians)   sin(arc2_end_angle_in_radians)  ];
% if arc2_is_counter_clockwise
%     cross_product_direction = 1;
% else
%     cross_product_direction = -1;
% end
% arc2_change_in_angle = fcn_geometry_findAngleUsing2PointsOnCircle([0 0],1, arc2_start_unit_vector, arc2_end_unit_vector, cross_product_direction);


% With the cleaned parameters, the line vector always
% points toward the joint of the line and the arc.
% to_joint_line_unit_tangent_vector = line_unit_tangent_vector;
% to_joint_line_unit_ortho_vector= to_joint_line_unit_tangent_vector*[0 1; -1 0];

% Use MATLAB's circcirc algorithm to find intersections between two circles
[xout,yout] = circcirc(arc1_center_xy(1,1),arc1_center_xy(1,2),arc1_radius,arc2_center_xy(1,1),arc2_center_xy(1,2),arc2_radius);

% Check results of above
if ~isempty(debug_fig_num)
    % Plot the intersection
    figure(debug_fig_num);
    subplot(3,2,subplot_number);

    % Plot the intersections
    intersections = [xout', yout'];
    plot(intersections(:,1),intersections(:,2),'k.','MarkerSize',20);
end

if ~any(isnan(xout))
    % intersection points were found! To be an intersection, the point must
    % be on both arc1 and arc2

    % Which point(s) to keep?
    circle_intersection_points = [xout', yout'];

    % Are the intersections within the arc range that we were given? To
    % check this, we use the three points on each arc - the start, the
    % intersection, and the end to calculate the arc direction. We then
    % check to see if it is the same as the given direction for that arc -
    % if it is, the point is on the arc. The way we search is to initialize
    % the potential arc intersection points to the circle intersections,
    % and remove any of the arc intersection points that are not on both of
    % the arcs.
    potential_arc_intersection_points = circle_intersection_points;
    for ith_row = 1:length(circle_intersection_points(:,1))
        
        % Check arc1
        intersection_is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(arc1_start_xy, circle_intersection_points(ith_row,:), arc1_end_xy,-1);
        if arc1_is_counter_clockwise ~= intersection_is_counterClockwise
            potential_arc_intersection_points(ith_row,:) = [nan nan];
        end

        % Check arc2
        intersection_is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(arc2_start_xy, circle_intersection_points(ith_row,:), arc2_end_xy,-1);
        if arc2_is_counter_clockwise ~= intersection_is_counterClockwise
            potential_arc_intersection_points(ith_row,:) = [nan nan];
        end

    end

    if isequal(potential_arc_intersection_points(1,:),potential_arc_intersection_points(2,:))
        intersection_point_both_arcs = potential_arc_intersection_points(1,:);
    else
        % Find which point is closest to arc1's start point
        if ~any(isnan(potential_arc_intersection_points(1,:)))
            arc_angle_point1  = fcn_geometry_arcAngleFrom3Points(arc1_start_xy, potential_arc_intersection_points(1,:), arc1_end_xy,(-1));
        else
            arc_angle_point1 = nan;
        end
        if ~any(isnan(potential_arc_intersection_points(2,:)))
            arc_angle_point2  = fcn_geometry_arcAngleFrom3Points(arc1_start_xy, potential_arc_intersection_points(2,:), arc1_end_xy,(-1));
        else
            arc_angle_point2 = nan;
        end

        if arc_angle_point1<arc_angle_point2
            intersection_point_both_arcs = potential_arc_intersection_points(1,:);
        else
            intersection_point_both_arcs = potential_arc_intersection_points(2,:);
        end
    end

else
    intersection_point_both_arcs = [nan nan];
end

if ~isempty(debug_fig_num)
    % Plot the intersection
    figure(debug_fig_num);
    subplot(3,2,1);
    debug_axis = axis;

    subplot(3,2,subplot_number);
    plot(intersection_point_both_arcs(:,1),intersection_point_both_arcs(:,2),'co','MarkerSize',10,'LineWidth',2);
    
    axis(debug_axis);
end

end % Ends fcn_INTERNAL_ArcArcIntersection

%% fcn_INTERNAL_convertParametersToStOrientation
function [st_arc1_parameters, st_arc2_parameters, St_transform_XYtoSt, flag_arc1_is_flipped] = ...
    fcn_INTERNAL_convertParametersToStOrientation(arc1_parameters, arc2_parameters, continuity_level, intersection_point, debug_fig_num)

% Calculate needed values from parameter sets
% Get the arc fit details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_center_xy                = arc1_parameters(1,1:2);
arc1_radius                   = arc1_parameters(1,3);
% arc1_start_angle_in_radians   = arc1_parameters(1,4);
arc1_end_angle_in_radians     = arc1_parameters(1,5);
% arc1_start_xy                 = arc1_center_xy + arc1_radius*[cos(arc1_start_angle_in_radians) sin(arc1_start_angle_in_radians)];
arc1_end_xy                   = arc1_center_xy + arc1_radius*[cos(arc1_end_angle_in_radians) sin(arc1_end_angle_in_radians)];
arc1_is_counter_clockwise     = arc1_parameters(1,7);

% Find the change in angle of the arc
% arc1_start_unit_vector        = [cos(arc1_start_angle_in_radians) sin(arc1_start_angle_in_radians)];
% arc1_end_unit_vector          = [cos(arc1_end_angle_in_radians)   sin(arc1_end_angle_in_radians)  ];
if arc1_is_counter_clockwise
    arc1_cross_product_direction = 1;
else
    arc1_cross_product_direction = -1;
end
% arc1_change_in_angle = fcn_geometry_findAngleUsing2PointsOnCircle([0 0],1, arc1_start_unit_vector, arc1_end_unit_vector, cross_product_direction);


% Get the arc fit details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_center_xy                = arc2_parameters(1,1:2);
arc2_radius                   = arc2_parameters(1,3);
arc2_start_angle_in_radians   = arc2_parameters(1,4);
% arc2_end_angle_in_radians     = arc2_parameters(1,5);
arc2_start_xy                 = arc2_center_xy + arc2_radius*[cos(arc2_start_angle_in_radians) sin(arc2_start_angle_in_radians)];
% arc2_end_xy                   = arc2_center_xy + arc2_radius*[cos(arc2_end_angle_in_radians) sin(arc2_end_angle_in_radians)];
arc2_is_counter_clockwise     = arc2_parameters(1,7);

% Find the change in angle of the arc
% arc2_start_unit_vector        = [cos(arc2_start_angle_in_radians) sin(arc2_start_angle_in_radians)];
% arc2_end_unit_vector          = [cos(arc2_end_angle_in_radians)   sin(arc2_end_angle_in_radians)  ];
if arc2_is_counter_clockwise
    arc2_cross_product_direction = 1;
else
    arc2_cross_product_direction = -1;
end
% arc2_change_in_angle = fcn_geometry_findAngleUsing2PointsOnCircle([0 0],1, arc2_start_unit_vector, arc2_end_unit_vector, cross_product_direction);


switch continuity_level
    case 0
        % For C0 continuity of arc1 to arc2, the closest point of the
        % desired joint to the arc is either the intersection of arc1 with
        % arc2, or where arc1 ends

        % Was an intersection found? If not, use the end of the arc1 as the
        % final point
         if any(isnan(intersection_point))
             % If enter here, no intersection was found
            angle_arc1_end = arc1_end_angle_in_radians;
         else
             % If enter here, an intersection was found

             %%%%
             % Find the angle where arc1 should end
             arc1_vector_center_to_intersection = intersection_point - arc1_center_xy;

             % Plot the vector (for debugging)?
             if ~isempty(debug_fig_num)
                 figure(debug_fig_num);
                 subplot(3,2,1);

                 % Plot the projection vector
                 subplot(3,2,2);
                 quiver(arc1_center_xy(1,1), arc1_center_xy(1,2), arc1_vector_center_to_intersection(1,1),  arc1_vector_center_to_intersection(1,2), 0,'LineWidth',3);
             end

             angle_arc1_end = atan2(arc1_vector_center_to_intersection(2),arc1_vector_center_to_intersection(1));

         end      

        % Perform the rotation
        desired_arc1_parameters = arc1_parameters;
        desired_arc1_parameters(1,5) = angle_arc1_end;

        secondary_parameters_type_strings{1} = 'arc';
        secondary_parameters{1}              = arc1_parameters;
        secondary_parameters_type_strings{2} = 'arc';
        secondary_parameters{2}              = arc2_parameters;

        [~, st_secondary_parameters, St_transform_XYtoSt, ~, flag_arc1_is_flipped] = ...
            fcn_geometry_orientGeometryXY2St('arc', desired_arc1_parameters, (secondary_parameters_type_strings), (secondary_parameters), (-1));

        st_arc1_parameters = st_secondary_parameters{1};
        st_arc2_parameters = st_secondary_parameters{2};

    case 1
        % For C1 continuity of arc1 to arc2, arc1 and arc2 are tangent to
        % each other. So the rotation will need to be the one that produces
        % arc1 such that it is tangent with the x-axis at exactly the point
        % where arc1 touches the tangent line to arc2. In some cases, this
        % will not exist and in these cases the equivalent line must be
        % found.

        % First, calculate all the tangent points on circles defined by
        % centers and radii for arc1 and arc2, where tangents can be either
        % inner tangents or outer tangents. The variable 
        % flag_inside_or_outside forces one or the other to be used. The
        % key factor of this function versus the one that follows is that
        % one can specify the voting points where the tangent is given that
        % is closest to the given votes. FORMAT:
        %
        % [...
        %     points_tangent_start, ...
        %     points_tangent_end] ...
        %     = ...
        %     fcn_geometry_findTangentPointsTwoCircles(...
        %     centers_start,...
        %     centers_end,...
        %     radii_start,...
        %     radii_end,...
        %     (flag_inside_or_out),...
        %     (voting_points_start,voting_points_end),...
        %     (fig_num))


        % Is the tangent on the inside or outside. For tangents to be on the
        % outside, they have to have the same clockwise or counterclockwise
        % orientation
        if (arc1_is_counter_clockwise==arc2_is_counter_clockwise)
            flag_inside_or_outside = 1;  % The tangent will be on the outside
        else
            flag_inside_or_outside = -1; % The tangent will be on the inside
        end

        [points_tangent_start1, points_tangent_end1] = ...
            fcn_geometry_findTangentPointsTwoCircles(...
            arc1_center_xy,...
            arc2_center_xy,...
            arc1_radius,...
            arc2_radius,...
            flag_inside_or_outside,...
            arc1_end_xy, ...
            arc2_start_xy,...
            -1);

        %%%%% OR

        % finds tangent points from one set of circles to another, returning only
        % the one set of points that matches the given cross products.
        % NOTE: the function below will NOT give the right answer if there are 2
        % possible tangents
        [...
            points_tangent_start2, ...
            points_tangent_end2] ...
            = ...
            fcn_geometry_findTangentPointTwoCircles(...
            arc1_center_xy,...
            arc2_center_xy,...
            arc1_radius,...
            arc2_radius,...
            arc1_cross_product_direction,...
            arc2_cross_product_direction,...
            -1);

        % Was a good tangent found? If it is good, there will be no NaN
        % values in the tangent points.
        if ~any(isnan([points_tangent_start1 points_tangent_start2 points_tangent_end1 points_tangent_end2]))
            % If enter here, the tangent points are good! Need to determine
            % which one to keep, since there are almost always multiple
            % ones

            % Make sure start and end points both agree, when calculated two different
            % ways. 
            rounded_start_values = round([points_tangent_start1 points_tangent_start2],4);
            assert(isequal(rounded_start_values(1,1:2),rounded_start_values(1,3:4)));
            rounded_end_values = round([points_tangent_end1 points_tangent_end2],4);
            assert(isequal(rounded_end_values(1,1:2),rounded_end_values(1,3:4)));

            %%%%
            % Find the angle where arc1 should end
            arc1_vector_center_to_arc_end = points_tangent_start1 - arc1_center_xy;

            % Plot the vector (for debugging)?
            if ~isempty(debug_fig_num)
                figure(debug_fig_num);
                subplot(3,2,1);

                % Plot the projection vector
                subplot(3,2,2);
                quiver(arc1_center_xy(1,1), arc1_center_xy(1,2), arc1_vector_center_to_arc_end(1,1),  arc1_vector_center_to_arc_end(1,2), 0,'LineWidth',3);
            end

            angle_arc1_end = atan2(arc1_vector_center_to_arc_end(2),arc1_vector_center_to_arc_end(1));

            % %%%%%%%%%%%%%%%%%
            % % Generate the base line
            % tangent_vector = points_tangent_end1 - points_tangent_start1;
            % 
            % tangent_vector_length = sum((tangent_vector).^2,2).^0.5;
            % 
            % 
            % circle_to_circle_unit_tangent_vector = fcn_geometry_calcUnitVector(tangent_vector);
            % circle_to_circle_unit_orthogo_vector = circle_to_circle_unit_tangent_vector*[0 1; -1 0];
            % 
            % line_angle = -atan2(circle_to_circle_unit_tangent_vector(2),circle_to_circle_unit_tangent_vector(1));
            % 
            % desired_circle1_center =  [0 0] + arc1_is_counter_clockwise*arc1_radius*circle_to_circle_unit_orthogo_vector;
            % desired_circle2_center =  [tangent_vector_length 0] + arc2_is_counter_clockwise*arc2_radius*[0 1];


        else
            % The code will enter here if the circle tangents are not
            % defined. This happens with inner tangents occur when arcs are
            % in opposite directions, for example counter-clockwise for
            % arc1 and clockwise for arc2 and if, in this case, the arcs
            % overlap. Inner tangents will not exist when circles intersect
            % each other AND the internal tangent is requested. If this is
            % the case, the alignment point on arc1 is on the
            % center-to-center line between the circles for arc1 and arc2.
            %
            % or
            %
            % circle tangents are also not defined with outer tangents when
            % one circle is inside the other circle, and the orientations
            % are the same. If this is the case, the alignment vector is
            % the one pointing "out" of first circle's direction, tangent
            % to the first circle, centered at the origin.

            circle1_to_circle2_projection_vector = arc2_center_xy - arc1_center_xy;

            % Plot the vector (for debugging)?
            if ~isempty(debug_fig_num)
                figure(debug_fig_num);
                subplot(3,2,1);

                % Plot the projection vector
                subplot(3,2,2);
                quiver(arc1_center_xy(1,1), arc1_center_xy(1,2), circle1_to_circle2_projection_vector(1,1),  circle1_to_circle2_projection_vector(1,2), 0,'LineWidth',3);
            end

            if arc1_is_counter_clockwise~=arc2_is_counter_clockwise
                % Should only enter here if an internal tangent is requested for
                % overlapping circles

                % The point of contact is defined as the first circle's radius
                % projected toward the second circle.
                angle_arc1_end = atan2(circle1_to_circle2_projection_vector(2),circle1_to_circle2_projection_vector(1));

            else
                % Should only enter here when one circle is inside another
                % The point of contact is defined as the first circle's radius
                % projected toward the second circle. 
                if arc1_radius>=arc2_radius
                    % If the arc1 radius is larger than arc2, then angle is
                    % toward projection
                    angle_arc1_end = atan2(circle1_to_circle2_projection_vector(2),circle1_to_circle2_projection_vector(1));
                else
                    % If the arc1 radius is larger than arc2, then angle is
                    % away from projection
                    angle_arc1_end = atan2(-1*circle1_to_circle2_projection_vector(2),-1*circle1_to_circle2_projection_vector(1));
                end
            end
        end

        % Perform the rotation
        desired_arc1_parameters = arc1_parameters;
        desired_arc1_parameters(1,5) = angle_arc1_end;

        secondary_parameters_type_strings{1} = 'arc';
        secondary_parameters{1}              = arc1_parameters;
        secondary_parameters_type_strings{2} = 'arc';
        secondary_parameters{2}              = arc2_parameters;

        [~, st_secondary_parameters, St_transform_XYtoSt, ~, flag_arc1_is_flipped] = ...
            fcn_geometry_orientGeometryXY2St('arc', desired_arc1_parameters, (secondary_parameters_type_strings), (secondary_parameters), (-1));

        st_arc1_parameters = st_secondary_parameters{1};
        st_arc2_parameters = st_secondary_parameters{2};

    case 2
        % For C2 continuity of an arc to an arc, the spiral must change
        % curvature. Calculation of the spiral is difficult and requires
        % numerical iteration, and the inputs require the calculation
        % of the offset of the spirals relative to each other - e.g. the
        % "space" between the arcs. There are many cases where spirals cannot
        % exist between arcs, and specifically if the space_between_circles is
        % negative, the spiral cannot exist. To determine if a spiral is
        % possible between arcs, we call the function to determine if the
        % spiral is possible between the circles containing the arcs, and if
        % not, how much space is needed to allow the spiral.
        secondary_parameters_type_strings{1} = 'arc';
        secondary_parameters{1}              = arc2_parameters;
        [st_primary_parameters, st_secondary_parameters, St_transform_XYtoSt, ~, flag_arc1_is_flipped] = ...
            fcn_geometry_orientGeometryXY2St('arc', arc1_parameters, (secondary_parameters_type_strings), (secondary_parameters), (-1));
        st_arc1_parameters = st_primary_parameters;
        st_arc2_parameters = st_secondary_parameters{1};

    otherwise
        error('This continuity not possible yet')
end

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

    fcn_geometry_plotCircle(st_arc1_parameters(1,1:2),st_arc1_parameters(1,3),...
        sprintf(' ''--'',''Color'',[0 0.6 0],''LineWidth'',1 '),debug_fig_num);
    fcn_geometry_plotCircle(st_arc2_parameters(1,1:2),st_arc2_parameters(1,3),...
        sprintf(' ''--'',''Color'',[0.6 0 0],''LineWidth'',1 '),debug_fig_num);

    fcn_geometry_plotGeometry('arc',st_arc1_parameters);
    fcn_geometry_plotGeometry('arc',st_arc2_parameters);

    plot(st_arc1_parameters(1,1),st_arc1_parameters(1,2),'+','Color',[0 0.6 0]);
    plot(st_arc2_parameters(1,1),st_arc2_parameters(1,2),'+','Color',[0.6 0 0]);

    axis(debug_axis);


    title('Rotated into St');
end 
end % Ends fcn_INTERNAL_convertParametersToStOrientation

%% fcn_INTERNAL_findShiftToMatchArcToArc
function [desired_arc1_parameters, desired_arc2_parameters, desired_intermediate_geometry_join_parameters, desired_intermediate_geometry_join_type] = ...
    fcn_INTERNAL_findShiftToMatchArcToArc(arc1_parameters, arc2_parameters, continuity_level, intersection_point, threshold, flag_perform_shift_of_arc2, debug_fig_num)
% Calculates the delta amount to match the arc to the arc. The delta
% values are measured FROM desired point TO actual point

if length(threshold)==1 || isempty(threshold)
    transverse_threshold = threshold;
else
    transverse_threshold = threshold(2);
end

desired_intermediate_geometry_join_parameters = nan(1,6); % Initialize output to be a "blank" spiral
desired_intermediate_geometry_join_type       = 'spiral';

% Calculate needed values from parameter sets
% Calculate needed values from parameter sets
% Get the arc fit details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_center_xy                = arc1_parameters(1,1:2);
arc1_radius                   = arc1_parameters(1,3);
arc1_start_angle_in_radians   = arc1_parameters(1,4);
% arc1_end_angle_in_radians     = arc1_parameters(1,5);
arc1_start_xy                 = arc1_center_xy + arc1_radius*[cos(arc1_start_angle_in_radians) sin(arc1_start_angle_in_radians)];
% arc1_end_xy                   = arc1_center_xy + arc1_radius*[cos(arc1_end_angle_in_radians) sin(arc1_end_angle_in_radians)];
arc1_is_counter_clockwise     = arc1_parameters(1,7);

% Find the change in angle of the arc
% arc1_start_unit_vector        = [cos(arc1_start_angle_in_radians) sin(arc1_start_angle_in_radians)];
% arc1_end_unit_vector          = [cos(arc1_end_angle_in_radians)   sin(arc1_end_angle_in_radians)  ];
if 1==arc1_is_counter_clockwise
    arc1_is_counter_clockwise = 1;
else
    arc1_is_counter_clockwise = -1;
end
% arc1_change_in_angle = fcn_geometry_findAngleUsing2PointsOnCircle([0 0],1, arc1_start_unit_vector, arc1_end_unit_vector, arc1_is_counter_clockwise);


% Get the arc fit details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_center_xy                = arc2_parameters(1,1:2);
arc2_radius                   = arc2_parameters(1,3);
% arc2_start_angle_in_radians   = arc2_parameters(1,4);
arc2_end_angle_in_radians     = arc2_parameters(1,5);
% arc2_start_xy                 = arc2_center_xy + arc2_radius*[cos(arc2_start_angle_in_radians) sin(arc2_start_angle_in_radians)];
arc2_end_xy                   = arc2_center_xy + arc2_radius*[cos(arc2_end_angle_in_radians) sin(arc2_end_angle_in_radians)];
flag_arc2_is_counter_clockwise     = arc2_parameters(1,7);

% Find the change in angle of the arc
% arc2_start_unit_vector        = [cos(arc2_start_angle_in_radians) sin(arc2_start_angle_in_radians)];
% arc2_end_unit_vector          = [cos(arc2_end_angle_in_radians)   sin(arc2_end_angle_in_radians)  ];
if 1==flag_arc2_is_counter_clockwise
    arc2_sign_counter_clockwise = 1;
else
    arc2_sign_counter_clockwise = -1;
end
% arc2_change_in_angle = fcn_geometry_findAngleUsing2PointsOnCircle([0 0],1, arc2_start_unit_vector, arc2_end_unit_vector, arc2_is_counter_clockwise);

% Calculate the distance between the circles and the join point
space_between_circles = fcn_geometry_gapCircleToCircle(arc1_radius, arc2_radius, arc2_center_xy, arc2_sign_counter_clockwise,-1);


% With the cleaned parameters, the line vector always
% points in the direction of the join point.
% to_joint_unit_tangent_vector = [1 0];
% to_joint_unit_ortho_vector   = [0 1];

switch continuity_level
    case 0
        % For C0 continuity of arc1 to arc2, the closest point of the
        % desired joint to the arc is simply the end of arc1

        desired_intermediate_geometry_join_type       = ''; % Intermediate geometry will be a line segement
        desired_intermediate_geometry_join_parameters = nan(1,6);

        if ~any(isnan(intersection_point))
            % Arc1 will end at the intersection, which is at the origin
            vector_from_arc1_center_to_intersection_joint = [0 0] - arc1_center_xy;
            angle_of_intersection_arc1 = atan2(vector_from_arc1_center_to_intersection_joint(2),vector_from_arc1_center_to_intersection_joint(1));

            vector_from_arc2_center_to_intersection_joint = [0 0] - arc2_center_xy;
            angle_of_intersection_arc2 = atan2(vector_from_arc2_center_to_intersection_joint(2),vector_from_arc2_center_to_intersection_joint(1));

            desired_arc1_parameters        = arc1_parameters;
            desired_arc1_parameters(1,5)   = angle_of_intersection_arc1; % Update where the spiral ends

            desired_arc2_parameters        = arc2_parameters;
            desired_arc2_parameters(1,4)   = angle_of_intersection_arc2; % Update where the spiral starts

        else

            if flag_perform_shift_of_arc2==1
                if abs(space_between_circles)<transverse_threshold
                    % Yes, arc2 can be moved enough to be tangent. So put
                    % arc2's center in correct place
                    desired_arc1_parameters        = arc1_parameters;
                    desired_arc1_parameters(1,5)   = -pi/2;

                    desired_arc2_parameters        = arc2_parameters;
                    if (1==arc2_sign_counter_clockwise)
                        desired_arc2_parameters(1,1:2) = [0 arc2_radius]; % Move arc2's center to align correctly
                        desired_arc2_parameters(1,4)   = -pi/2;
                    else
                        desired_arc2_parameters(1,1:2) = [0 -arc2_radius]; % Move arc2's center to align correctly
                        desired_arc2_parameters(1,4)   = pi/2;
                    end


                else
                    % Not possible to shift enough to allow connection
                    desired_arc1_parameters        = nan(size(arc1_parameters));
                    desired_arc2_parameters        = nan(size(arc1_parameters));
                end
            else
                % No shift allowed by user entry, so not possible
                desired_arc1_parameters            = nan(size(arc1_parameters));
                desired_arc2_parameters            = nan(size(arc1_parameters));
            end
        end

    case 1
        % For C1 continuity of arc1 to arc2, the closest point of the desired
        % joint to arc1 and arc2 is always going to be the point where each arc
        % touches the tangent line connecting them. Since arc2 is the only one
        % that can be moved, the connecting point for arc1 is simply the location where
        % arc1 is tangent, which by construction is [0 0]. The connection
        % point for arc2 is where the arc2's circle touches the x-axis,
        % which by construction is the radius distance below the center of
        % the circle

        if (space_between_circles>0) && (1==flag_arc2_is_counter_clockwise)
            % In this case, one circle is inside the other and an external
            % tangent is requested - and this is not possible. So need to
            % check to see if arc2 can be moved to create an external
            % tangent.
            if flag_perform_shift_of_arc2==1
                if abs(space_between_circles)<transverse_threshold
                    % Yes, arc2 can be moved enough to be tangent. So put
                    % arc2's center in correct place
                    desired_arc1_parameters        = arc1_parameters;
                    desired_arc1_parameters(1,5)   = -pi/2;

                    desired_arc2_parameters        = arc2_parameters;
                    desired_arc2_parameters(1,1:2) = [0 arc2_radius]; % Move arc2's center to align correctly
                    desired_arc2_parameters(1,4)   = -pi/2;

                    segment_parameters(1,1:2) = [1 0]; % line_unit_tangent_vector
                    segment_parameters(1,3:4) = [0 0]; % line_base_point_xy
                    segment_parameters(1,5)   = 0;     % line_s_start
                    segment_parameters(1,6)   = abs(arc2_center_xy(1)); % line_s_end
                    desired_intermediate_geometry_join_type       = 'segment'; % Intermediate geometry will be a line segement
                    desired_intermediate_geometry_join_parameters = segment_parameters; 
                   
                else
                    % Not possible to shift enough to allow connection
                    desired_arc1_parameters        = nan(size(arc1_parameters));
                    desired_arc2_parameters        = nan(size(arc1_parameters));
                    desired_intermediate_geometry_join_type       = 'segment'; % Intermediate geometry will be a line segement
                    desired_intermediate_geometry_join_parameters = nan(1,6);  % "empty" line 
                end
            else
                % No shift allowed by user entry, so not possible
                desired_arc1_parameters            = nan(size(arc1_parameters));
                desired_arc2_parameters            = nan(size(arc1_parameters));
                desired_intermediate_geometry_join_type       = 'segment'; % Intermediate geometry will be a line segement
                desired_intermediate_geometry_join_parameters = nan(1,6);  % "empty" line
            end
        elseif (space_between_circles<=0) && (1==flag_arc2_is_counter_clockwise)
            % In this case, one circle is OUTSIDE the other and an external
            % tangent is requested - and this is always possible.
            
            % Put arcs start/end angles in correct places
            desired_arc1_parameters        = arc1_parameters;
            desired_arc1_parameters(1,5)   = -pi/2; % End of arc 1 should be at -90 degrees

            desired_arc2_parameters        = arc2_parameters;
            desired_arc2_parameters(1,4)   = -pi/2; % Start of arc 2 should be at -90 degrees

            segment_parameters(1,1:2) = [1 0]; % line_unit_tangent_vector
            segment_parameters(1,3:4) = [0 0]; % line_base_point_xy
            segment_parameters(1,5)   = 0;     % line_s_start
            segment_parameters(1,6)   = abs(arc2_center_xy(1)); % line_s_end
            desired_intermediate_geometry_join_type       = 'segment'; % Intermediate geometry will be a line segement
            desired_intermediate_geometry_join_parameters = segment_parameters;

            
        elseif (space_between_circles>0)  && (1~=flag_arc2_is_counter_clockwise)
            % The circles do not overlap, and an outside tangent is
            % requested. This is always possible.
            
            % Put arcs start/end angles in correct places
            desired_arc1_parameters        = arc1_parameters;
            desired_arc1_parameters(1,5)   = -pi/2; % End of arc 1 should be at -90 degrees

            desired_arc2_parameters        = arc2_parameters;
            desired_arc2_parameters(1,4)   = pi/2; % Start of arc 2 should be at 90 degrees

            segment_parameters(1,1:2) = [1 0]; % line_unit_tangent_vector
            segment_parameters(1,3:4) = [0 0]; % line_base_point_xy
            segment_parameters(1,5)   = 0;     % line_s_start
            segment_parameters(1,6)   = abs(arc2_center_xy(1)); % line_s_end
            desired_intermediate_geometry_join_type       = 'segment'; % Intermediate geometry will be a line segement
            desired_intermediate_geometry_join_parameters = segment_parameters;
            

        elseif (space_between_circles<=0) && (1~=flag_arc2_is_counter_clockwise)
            % The circles overlap, and an outside tangent is requested.
            % This isn't possible unless arc2 can be moved
            if flag_perform_shift_of_arc2==1
                if abs(space_between_circles)<transverse_threshold
                    % Yes, arc2 can be moved enough to be tangent. So put
                    % arc2's arc start and center in correct place.
                    desired_arc1_parameters        = arc1_parameters;
                    desired_arc1_parameters(1,5)   = -pi/2;

                    desired_arc2_parameters        = arc2_parameters;
                    desired_arc2_parameters(1,1:2) = [0 -arc2_radius]; % Move arc2's center to align correctly
                    desired_arc2_parameters(1,4)   = pi/2; % Start of arc 2 should be at 90 degrees
            
                    segment_parameters(1,1:2) = [1 0]; % line_unit_tangent_vector
                    segment_parameters(1,3:4) = [0 0]; % line_base_point_xy
                    segment_parameters(1,5)   = 0;     % line_s_start
                    segment_parameters(1,6)   = 0; % line_s_end
                    desired_intermediate_geometry_join_type       = 'segment'; % Intermediate geometry will be a line segement
                    desired_intermediate_geometry_join_parameters = segment_parameters;

                else
                    % Not possible to shift enough to allow connection
                    desired_arc1_parameters        = nan(size(arc1_parameters));
                    desired_arc2_parameters        = nan(size(arc1_parameters));
                    desired_intermediate_geometry_join_type       = 'segment'; % Intermediate geometry will be a line segement
                    desired_intermediate_geometry_join_parameters = nan(1,6);  % "empty" line
                end
            else
                % No shift allowed by user entry, so not possible
                desired_arc1_parameters            = nan(size(arc1_parameters));
                desired_arc2_parameters            = nan(size(arc1_parameters));
                desired_intermediate_geometry_join_type       = 'segment'; % Intermediate geometry will be a line segement
                desired_intermediate_geometry_join_parameters = nan(1,6);  % "empty" line
            end

        else
            error('Unknown error encountered - it should not be possible to enter this case!');            
        end





    case 2
        % For C2 continuity of an arc to an arc, the spiral must change
        % curvature. Calculation of the spiral is difficult and requires
        % numerical iteration, and the inputs require the calculation
        % of the offset of the spirals relative to each other - e.g. the
        % "space" between the arcs. There are many cases where spirals cannot
        % exist between arcs, and specifically if the space_between_circles is
        % negative, the spiral cannot exist. To determine if a spiral is
        % possible between arcs, we call the function to determine if the
        % spiral is possible between the circles containing the arcs, and if
        % not, how much space is needed to allow the spiral.


        if space_between_circles>0
            % spiral_join_parameters = [spiralLength,h0,x0,y0,K0,Kf];

            circle1_parameters = [0 arc1_radius  arc1_radius];
            circle2_parameters = [arc2_center_xy arc2_radius];

            [desired_intermediate_geometry_join_parameters, ~] = ...
                fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, arc2_sign_counter_clockwise, -1);

            % 'spiral' -
            % 
            % [
            %     x0,  % The initial x value
            %     y0,  % The initial y value
            %     h0,  % The initial heading
            %     s_Length,  % the s-coordinate length allowed
            %     K0,  % The initial curvature
            %     Kf   % The final curvature
            %     ]

            h0           = desired_intermediate_geometry_join_parameters(1,3);
            spiralLength = desired_intermediate_geometry_join_parameters(1,4);
            % x0           = desired_spiral_join_parameters(1,3);
            % y0           = desired_spiral_join_parameters(1,4);
            K0           = desired_intermediate_geometry_join_parameters(1,5);
            Kf           = desired_intermediate_geometry_join_parameters(1,6);

            % Check results?
            if 1==0

                % Plot the results in
                % figure 1234
                figure(1234);
                clf;
                hold on;
                grid on;
                axis equal;

                fcn_geometry_plotGeometry('arc',arc1_parameters);
                fcn_geometry_plotGeometry('arc',arc2_parameters);
                fcn_geometry_plotGeometry('spiral',desired_intermediate_geometry_join_parameters);
            end % Ends plotting

            % Find the angle and position that the spiral ends at
            analytical_end_angle   = h0 + (Kf-K0)*spiralLength/2 + K0*spiralLength;

            % The checks should confirm that the
            % arc's start and end is after the start of arc1 (not within) and before the
            % end of arc2 (not within)
            flag_spiral_is_bad = 0;

            % Precalculate some quantities
            arc1_angle_where_spiral_starts = h0 - pi/2;
            spiral_arc1_join_xy = arc1_center_xy + arc1_radius*[cos(arc1_angle_where_spiral_starts) sin(arc1_angle_where_spiral_starts)];
            if 1==flag_arc2_is_counter_clockwise
                arc2_angle_where_spiral_ends = analytical_end_angle - pi/2;
            else
                arc2_angle_where_spiral_ends = pi/2+analytical_end_angle;
            end
            spiral_arc2_join_xy = arc2_center_xy + arc2_radius*[cos(arc2_angle_where_spiral_ends) sin(arc2_angle_where_spiral_ends)];


            % Make sure the spiral's start position is within the angle range
            % of arc1's start to the spiral end. If not, return Nan values for everything

            intersection_is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(arc1_start_xy, spiral_arc1_join_xy, spiral_arc2_join_xy,-1);
            if arc1_is_counter_clockwise ~= intersection_is_counterClockwise
                % St_shift = [nan nan];
                desired_arc1_parameters = nan(size(arc1_parameters));
                desired_arc2_parameters = nan(size(arc1_parameters));
                flag_spiral_is_bad = 1;
            end

            % Make sure the spiral's end position is within the angle range
            % of arc2's end. If not, return Nan values for everything

            intersection_is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(spiral_arc1_join_xy, spiral_arc2_join_xy, arc2_end_xy,-1);
            if arc2_sign_counter_clockwise ~= intersection_is_counterClockwise
                % St_shift = [nan nan];
                desired_arc1_parameters = nan(size(arc1_parameters));
                desired_arc2_parameters = nan(size(arc1_parameters));
                flag_spiral_is_bad = 1;
            end

            if 0==flag_spiral_is_bad
                % If made it here, then the spiral is good. Fill outputs and then
                % return results
                desired_arc1_parameters = arc1_parameters;
                desired_arc1_parameters(1,4)   = arc1_start_angle_in_radians;
                desired_arc1_parameters(1,5)   = arc1_angle_where_spiral_starts;

                desired_arc2_parameters        = arc2_parameters;
                desired_arc2_parameters(1,4)   = arc2_angle_where_spiral_ends;
                desired_arc2_parameters(1,5)   = arc2_end_angle_in_radians;
            end
        else
            % Not enough space between circles for a spiral, not without
            % modifying it
            if flag_perform_shift_of_arc2==1
                if abs(space_between_circles)<transverse_threshold
                    direction_vector_to_shift_arc2 = arc2_center_xy - arc1_center_xy;
                    unit_direction_vector_to_shift_arc2 = fcn_geometry_calcUnitVector(direction_vector_to_shift_arc2);

                    if 1==flag_arc2_is_counter_clockwise
                        % Small circle must be completely inside the large circle
                        sign_of_vector = -1;
                    else
                        % Small circle must be completely outside the large circle
                        sign_of_vector = 1;
                    end

                    revised_arc2_parameters = arc2_parameters;
                    revised_arc2_parameters(1,1:2) = arc2_parameters(1,1:2) + unit_direction_vector_to_shift_arc2*sign_of_vector*(abs(space_between_circles)+0.001);

                    if ~isempty(debug_fig_num)
                        plotting_fig_num = 454545;
                        figure(plotting_fig_num);
                        clf;

                        subplot(1,2,1);
                        hold on;
                        axis equal
                        grid on;
                        fcn_geometry_plotCircle(arc1_parameters(1,1:2),arc1_parameters(1,3),...
                            sprintf(' ''--'',''Color'',[0 0.6 0],''LineWidth'',1 '),plotting_fig_num);
                        fcn_geometry_plotCircle(arc2_parameters(1,1:2),arc2_parameters(1,3),...
                            sprintf(' ''--'',''Color'',[0.6 0 0],''LineWidth'',1 '),plotting_fig_num);
                        fcn_geometry_plotGeometry('arc',arc1_parameters);
                        fcn_geometry_plotGeometry('arc',arc2_parameters);

                        subplot(1,2,2);
                        hold on;
                        axis equal
                        grid on;
                        fcn_geometry_plotCircle(arc1_parameters(1,1:2),arc1_parameters(1,3),...
                            sprintf(' ''--'',''Color'',[0 0.6 0],''LineWidth'',1 '),plotting_fig_num);
                        fcn_geometry_plotCircle(revised_arc2_parameters(1,1:2),revised_arc2_parameters(1,3),...
                            sprintf(' ''--'',''Color'',[0.6 0 0],''LineWidth'',1 '),plotting_fig_num);
                        fcn_geometry_plotGeometry('arc',arc1_parameters);
                        fcn_geometry_plotGeometry('arc',revised_arc2_parameters);
                    end

                    [desired_arc1_parameters, desired_arc2_parameters, desired_intermediate_geometry_join_parameters] = ...
                        fcn_INTERNAL_findShiftToMatchArcToArc(arc1_parameters, revised_arc2_parameters, continuity_level, intersection_point, [], 0, debug_fig_num);
                else
                    % Not possible
                    desired_arc1_parameters = nan(size(arc1_parameters));
                    desired_arc2_parameters = nan(size(arc1_parameters));
                end
            else
                % Not allowed by user entry
                desired_arc1_parameters = nan(size(arc1_parameters));
                desired_arc2_parameters = nan(size(arc1_parameters));
            end
        end % Ends if statement to check if spiral is possible
    otherwise
        error('This continuity not possible yet')
end

% % Depending on which side of the line the arc is located, either need to
% % subtract or add on the circle radius. Note: the result is designed to be
% % positive if the line doesn't quite meet the circle as a tangent creating
% % a gap - in this case, a spiral is needed. The result is negative if the
% % line would intersect the circle; in this case, a line shift offset is
% % needed.
% unit_St_vector_to_arc2_center = to_joint_unit_ortho_vector*arc2_is_counter_clockwise;
% if 0 == continuity_level
%     difference_vector_from_line_end_to_arc_start = arc2_start_xy - desired_closest_arc1_point_to_joint;
%     delta_transverse = sum(to_joint_unit_ortho_vector.*difference_vector_from_line_end_to_arc_start,2);
%     delta_station    = sum(to_joint_unit_tangent_vector.*difference_vector_from_line_end_to_arc_start,2);
%     St_shift = [delta_station, delta_transverse];
% elseif 1 == continuity_level
%     desired_arc2_center_xy = desired_closest_arc1_point_to_joint + arc2_radius*unit_St_vector_to_arc2_center;
%     difference_vector_from_desired_to_actual_circle_center = arc2_center_xy - desired_arc2_center_xy;
%     delta_transverse = sum(to_joint_unit_ortho_vector.*difference_vector_from_desired_to_actual_circle_center,2);
%     delta_station    = sum(to_joint_unit_tangent_vector.*difference_vector_from_desired_to_actual_circle_center,2);
%     % if signed_distance_ortho_line_to_arc_center < 0
%     %     delta_transverse   = -1*delta_transverse;
%     % end
%     St_shift = [delta_station, delta_transverse];
% elseif 2 == continuity_level
%     % Do nothing
% else
%     error('This continuity not possible yet')
% end


if ~isempty(debug_fig_num)
    % Plot the bounding box
    figure(debug_fig_num);
    subplot(3,2,1);
    debug_axis = axis;

    subplot(3,2,4);
    fcn_geometry_plotCircle(desired_arc1_parameters(1,1:2),desired_arc1_parameters(1,3),...
        sprintf(' ''--'',''Color'',[0 0.6 0],''LineWidth'',1 '),debug_fig_num);
    fcn_geometry_plotCircle(desired_arc2_parameters(1,1:2),desired_arc2_parameters(1,3),...
        sprintf(' ''--'',''Color'',[0.6 0 0],''LineWidth'',1 '),debug_fig_num);

    fcn_geometry_plotGeometry('arc',desired_arc1_parameters);
    fcn_geometry_plotGeometry('arc',desired_arc2_parameters);
    fcn_geometry_plotGeometry(desired_intermediate_geometry_join_type,desired_intermediate_geometry_join_parameters);
    
    axis(debug_axis);

    hold on;
    grid on;
    axis equal;
    xlabel('S [meters]');
    ylabel('t [meters]')

    title('Desired St geometries');
end
end % Ends fcn_INTERNAL_findShiftToMatchArcToArc

%% fcn_INTERNAL_performShift
function [revised_arc1_parameters_St,revised_arc2_parameters_St, revised_intermediate_geometry_type, revised_intermediate_geometry_parameters_St] = ...
    fcn_INTERNAL_performShift(threshold, continuity_level, ...
    st_arc1_parameters, st_arc2_parameters, ...
    desired_st_arc1_parameters, desired_st_arc2_parameters, ...
    desired_intermediate_geometry_join_parameters, desired_intermediate_geometry_join_type, ...
    debug_fig_num)

% Find out how much arc2 is shifting by looking at how much the center of
% arc2 is moving
St_shift = st_arc2_parameters(1,1:2)-desired_st_arc2_parameters(1,1:2);

% Check to see if shift is even possible
flag_shift_is_possible = 0;

% The threshold can be a [1 x 2] representing S and t tolerances or [1 x 1]
% representing total distances. Check each.
if length(threshold)==1
    shift_distance = abs(sum(St_shift.^2,2).^0.5);
    if shift_distance<=threshold
        flag_shift_is_possible = 1;
    end
else
    if abs(St_shift(1))<=threshold(1) && abs(St_shift(2))<=threshold(2)
        flag_shift_is_possible = 1;
    end

end
if 0==flag_shift_is_possible    
    % Not possible to shift
    revised_arc1_parameters_St         = nan(size(st_arc1_parameters));
    revised_arc2_parameters_St         = nan(size(st_arc2_parameters));
    revised_intermediate_geometry_parameters_St = nan(size(desired_intermediate_geometry_join_parameters));
    revised_intermediate_geometry_type = desired_intermediate_geometry_join_type; 
else
    revised_arc1_parameters_St        = desired_st_arc1_parameters;
    revised_arc2_parameters_St        = desired_st_arc2_parameters;
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
            revised_intermediate_geometry_type = 'segment';
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

    fcn_geometry_plotGeometry('arc',revised_arc1_parameters_St);
    fcn_geometry_plotGeometry('arc',revised_arc2_parameters_St);
    fcn_geometry_plotGeometry(desired_intermediate_geometry_join_type,revised_intermediate_geometry_parameters_St);

    title('St outputs after shift');
    axis(debug_axis);
end

end % Ends fcn_INTERNAL_performShift

%% fcn_INTERNAL_convertParametersOutOfStOrientation
function [revised_arc1_parameters, revised_arc2_parameters, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters] = ...
    fcn_INTERNAL_convertParametersOutOfStOrientation(...
    revised_arc1_parameters_St, revised_arc2_parameters_St, revised_intermediate_geometry_join_type, revised_intermediate_geometry_join_parameters_St, St_transform_XYtoSt, flag_arc1_is_flipped, debug_fig_num)

% Call the function to convert from ST back to XY
st_parameters_type_strings{1} = 'arc';
st_parameters_type_strings{2} = 'arc';
st_parameters_type_strings{3} = revised_intermediate_geometry_join_type;
st_parameters{1} = revised_arc1_parameters_St;
st_parameters{2} = revised_arc2_parameters_St;
st_parameters{3} = revised_intermediate_geometry_join_parameters_St;

[XY_parameters] = ...
fcn_geometry_orientGeometrySt2XY(st_parameters_type_strings, st_parameters, St_transform_XYtoSt, flag_arc1_is_flipped, (-1));

revised_arc1_parameters = XY_parameters{1};
revised_arc2_parameters = XY_parameters{2};
revised_intermediate_geometry_join_parameters = XY_parameters{3};

if ~isempty(debug_fig_num)
    % Plot the results
    figure(debug_fig_num);
    subplot(3,2,1);
    debug_axis = axis;

    subplot(3,2,6);

    fcn_geometry_plotGeometry('arc',revised_arc1_parameters);
    fcn_geometry_plotGeometry('arc',revised_arc2_parameters);
    fcn_geometry_plotGeometry(revised_intermediate_geometry_join_type,revised_intermediate_geometry_join_parameters);
    
    title('Final outputs');
    axis(debug_axis);
end
end % Ends fcn_INTERNAL_convertParametersOutOfStOrientation

%% fcn_INTERNAL_prepDebugFigure
function fcn_INTERNAL_prepDebugFigure(arc1_parameters, arc2_parameters, debug_fig_num)
if ~isempty(debug_fig_num)
    figure(debug_fig_num);
    clf;

    % Plot the inputs
    subplot(3,2,1);


    fcn_geometry_plotCircle(arc1_parameters(1,1:2),arc1_parameters(1,3),...
        sprintf(' ''--'',''Color'',[0 0.6 0],''LineWidth'',1 '),debug_fig_num);
    fcn_geometry_plotCircle(arc2_parameters(1,1:2),arc2_parameters(1,3),...
        sprintf(' ''--'',''Color'',[0.6 0 0],''LineWidth'',1 '),debug_fig_num);


    fcn_geometry_plotGeometry('arc',arc1_parameters);
    fcn_geometry_plotGeometry('arc',arc2_parameters);

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



