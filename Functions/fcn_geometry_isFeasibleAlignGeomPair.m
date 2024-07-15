function [flag_is_feasible, feasibility_distance, closest_feasible_input2_parameters] = ...
    fcn_geometry_isFeasibleAlignGeomPair(input1_type, input1_parameters, input2_type, input2_parameters, continuity_level, varargin)
%% fcn_geometry_isFeasibleAlignGeomPair
% Given two geometries and a given type of continuity, checks if it is
% possible to achieve user-requested levels of continuity between the two.
%
% The following matches are supported with either C0, C1, and C2
% continuity, and user-defined tolerances:
%
% * connections from line segments to arcs using
%   fcn_geometry_alignSegmentArc 
%
% * connections from arcs to line segments using
%   fcn_geometry_alignArcSegment
%
% * connections from line arcs to arcs, using
%   fcn_geometry_alignArcArc 
%
% Format:
% [flag_is_feasible, feasibility_distance] = ...
%    fcn_geometry_isFeasibleAlignGeomPair(...
%    input1_type, input1_parameters, ...
%    input2_type, input2_parameters, ...
%    continuity, (threshold), (fig_num))
%
% INPUTS:
%      input1_type: a string denoting the first geometry type, consisting
%      of: 'circle','arc','line','segment'
%
%      input1_parameters: a vector of the parameters consistent with the
%      input type string. See fcn_geometry_fillEmptyDomainStructure for
%      details
%
%      input2_type: a string denoting the first geometry type, consisting
%      of: 'circle','arc','line','segment'
%
%      input2_parameters: a vector of the parameters consistent with the
%      input type string. See fcn_geometry_fillEmptyDomainStructure for
%      details
%
%      continuity_level: the level of continuity being sought between adjacent
%      geometries. The levels of continuity include:
%
%      C0 continuous: the ith+1 segment is forced to intersect the ith segment.
%      For two segments to have C0 continuity in agreement, they must
%      intersect at at least 1 point.
%
%      C1 continous: the ith+1 segment is C0 continous with the ith segment,
%      and the tangent vector of the ith+1 segment is forced to match the
%      tangent vector of the ith segment at the intersection point. In other
%      words, the derivative dt/ds in St-coordinates are matched. Or, in XY
%      coordinates, the radius of ith fit at point of contact is matched to the
%      radius of the ith+1 fit at the same contact point. NOTE: a function can
%      be C1 continuous in St coordinates but NOT C1 continuous in XY
%      coordinates. For example, an arc joining a vertical line has an
%      undefined slope at the vertical line, but still be tangent to the
%      circle.
%
%      C2 continous: the ith+1 segment is C0 and C1 continous with the ith
%      segment, and the curvature (1/R) of the ith+1 segment is forced to match
%      the curvature of the ith segment at the intersection point. In other
%      words, the 2nd derivative dt/ds in St-coordinates are matched such that
%      there are no discontinuities in the curvature. For vehicles, steering
%      angle is kinematically (but not dynamically) associated with the
%      radius of the local path. Enforcing C2 curvature thus forces the curve
%      fit to produce a path such that the steering angle of a vehicle would
%      not have to instantaneously change from one setting to another in order
%      for a vehicle to remain on a path.
%
%      (OPTIONAL INPUTS)
%
%      threshold: the offset, in meters, between the arc and the line such
%      that this offset is removed by shifting. If the offset is larger
%      than this, then the outputs are set to empty. If this is entered as
%      a 2x1 or 1x2, then this specifies the threshold in St coordinates,
%      e.g. first in the station direction, and then in the transverse
%      direction. For example, an entry of [3 0.2] would have 3 meters
%      threshold in the station direction, but 0.2 meters threshold in the
%      transverse direction. Default is 0.
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      revised_fitSequence_parameters: the parameters for each of the N
%      fits such that they are aligned.
%
%      closest_feasible_input2_parameters: the closest parameters that give
%      feasible continuity.
%
% DEPENDENCIES:
%
%      fcn_geometry_fitArcRegressionFromHoughFit
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_isFeasibleAlignGeomPair
% for a full test suite.
%
% This function was written on 2024_04_19 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_04_19 - S Brennan
% -- wrote the code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==7 && isequal(varargin{end},-1))
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
        narginchk(5,7);

    end
end


% Does user want to specify threshold?
threshold = 0;
if (6<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        threshold = temp;
    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if  (0==flag_max_speed) && (7<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;

        % Does user want to specify animation_figure_handles?
        % flag_plot_subfigs = 0;

        if length(temp)>1
            fig_num           = temp(1);
            % h_plotPoints      = temp(2);
            % h_plotPercentage  = temp(3);
            % h_plotFitShape    = temp(4);
            % flag_plot_subfigs = 1;
        end
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

input1_type = fcn_INTERNAL_covertComplexShapeNamesToSimpleNames(input1_type);
input2_type = fcn_INTERNAL_covertComplexShapeNamesToSimpleNames(input2_type);

switch input1_type
    case 'segment'
        [flag_is_feasible, feasibility_distance, closest_feasible_input2_parameters] = fcn_INTERNAL_confirmSegmentToX(input1_parameters, input2_type, input2_parameters, continuity_level, threshold);
    case 'arc'
        [flag_is_feasible, feasibility_distance, closest_feasible_input2_parameters] = fcn_INTERNAL_confirmArcToX(input1_parameters, input2_type, input2_parameters, continuity_level, threshold);
    otherwise
        warning('on','backtrace');
        warning('An error will be thrown at this point due to missing code.');
        error('Alignments are not yet supported for curves from fit type: %s',current_fit_type);
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

    % if flag_plot_subfigs
    %
    %     figure(get(h_plotFitShape.Parent.Parent, 'Number'));
    %
    %     % Plot the results in the subplot
    %     flag_rescale_axis = 0;
    %
    %     % Match subplot 4 axis with that from subplot 1
    %     subplot(2,2,1);
    %     original_axis = axis;
    %
    % else
    % Plot the results in the given figure number
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end
    % end

    % if flag_plot_subfigs == 1
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     % Final fit
    %     figure(fig_num);
    %     subplot(2,2,4);
    %     hold on;
    %     grid on;
    %     axis equal;
    %     xlabel('X [meters]');
    %     ylabel('Y [meters]');
    %
    %     % Plot the fit results
    %     for ith_domain = 1:length(fitSequence_points)
    %         % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain,fitSequence_bestFitType{ith_domain},[],-1);
    %         current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain,[],[],-1);
    %         current_fitSequence_points = fitSequence_points{ith_domain};
    %         current_fitSequence_shape  = fitSequence_shapes{ith_domain};
    %         plot(current_fitSequence_points(:,1),current_fitSequence_points(:,2),'.','Color',current_color*0.8,'MarkerSize',10);
    %         plot(current_fitSequence_shape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',1,'EdgeAlpha',0);
    %     end
    %
    %     % Plot the domain fits
    %     fcn_geometry_plotFitSequences(fitSequence_bestFitType, fitSequence_parameters,(fig_num));
    %
    %     axis(original_axis);
    %
    %
    % else


    % Plot the input geometries 
    % subplot(1,2,1);
    hold on;
    grid on;
    axis equal;
    xlabel('X [meters]');
    ylabel('Y [meters]');

    fcn_geometry_plotFitSequences({input1_type, input2_type},{input1_parameters,input2_parameters},(fig_num));

    title(sprintf('Is continuity %.0d feasible?: %.0d, dist: %.2f',continuity_level,flag_is_feasible,feasibility_distance));


    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end

    % good_axis = axis;
    %
    %  % Plot the aligned geometries on the left
    % subplot(1,2,2);
    % hold on;
    % grid on;
    % axis equal;
    % xlabel('X [meters]');
    % ylabel('Y [meters]');
    %
    % fcn_geometry_plotFitSequences(revised_fitSequence_types, revised_fitSequence_parameters,(fig_num));
    %
    % % Match axis
    % axis(good_axis);

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
%% fcn_INTERNAL_covertComplexShapeNamesToSimpleNames
function simple_name_string = fcn_INTERNAL_covertComplexShapeNamesToSimpleNames(complex_name_string)

switch lower(complex_name_string)
    case {'arc','line','segment','spiral','','none','circle'}
        simple_name_string = complex_name_string;
    case {'regression arc'}
        simple_name_string = 'arc';
    case {'vector regression segment fit'}
        simple_name_string = 'segment';
    otherwise
        warning('on','backtrace');
        warning('An error will be thrown due to unrecognized fitting name type, inside fcn_INTERNAL_covertComplexShapeNamesToSimpleNames.');
        error('Unrecognized fit string: %s', complex_name_string);
end

end % Ends fcn_INTERNAL_covertComplexShapeNamesToSimpleNames


%% fcn_INTERNAL_confirmArcToX
function [flag_is_feasible, feasibility_distance, closest_feasible_input2_parameters] = fcn_INTERNAL_confirmArcToX(arc_parameters, X_fitType, X_parameters, continuity_level, threshold)

switch X_fitType
    case 'segment'
        % Arc to Segment - same as Segment to Arc
        [flag_is_feasible, feasibility_distance, closest_feasible_input2_parameters] = fcn_INTERNAL_confirmSegmentToArc(X_parameters, arc_parameters, continuity_level, threshold);

    case 'arc'
        % Arc to Arc
        [flag_is_feasible, feasibility_distance, closest_feasible_input2_parameters] = fcn_INTERNAL_confirmArcToArc(arc_parameters, X_parameters, continuity_level, threshold);

    otherwise
        warning('on','backtrace');
        warning('An error will be thrown at this point due to missing code.');
        error('Alignments are not yet supported for curves from fit type: %s',current_fit_type);
end
end % Ends fcn_INTERNAL_confirmArcToX


%% fcn_INTERNAL_confirmSegmentToX
function [flag_is_feasible, feasibility_distance, closest_feasible_input2_parameters] = fcn_INTERNAL_confirmSegmentToX(segment_parameters, X_fitType, X_parameters, continuity_level, threshold)

switch X_fitType
    case 'segment'
        % Segment to Segment
        warning('on','backtrace');
        warning('An error will be thrown at this point due to missing code.');
        error('Alignments from segment to segment are not yet supported.');

    case 'arc'
        % Segment to Arc
        [flag_is_feasible, feasibility_distance, closest_feasible_input2_parameters] = fcn_INTERNAL_confirmSegmentToArc(segment_parameters, X_parameters, continuity_level, threshold);

    otherwise
        warning('on','backtrace');
        warning('An error will be thrown at this point due to missing code.');
        error('Alignments are not yet supported for curves from fit type: %s',current_fit_type);
end

end % Ends fcn_INTERNAL_confirmSegmentToX

%% fcn_INTERNAL_confirmSegmentToArc
function  [flag_is_feasible, feasibility_distance, closest_feasible_input2_parameters] = fcn_INTERNAL_confirmSegmentToArc(segment_parameters, arc_parameters, continuity_level, threshold)


% Calculate needed values from parameter sets
% Get the arc fit details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy                = arc_parameters(1,1:2);
arc_radius                   = arc_parameters(1,3);
% arc_start_angle_in_radians   = arc_parameters(1,4);
% arc_end_angle_in_radians     = arc_parameters(1,5);
% arc_start_xy                 = arc_center_xy + arc_radius*[cos(arc_start_angle_in_radians) sin(arc_start_angle_in_radians)];
% arc_end_xy                   = arc_center_xy + arc_radius*[cos(arc_end_angle_in_radians) sin(arc_end_angle_in_radians)];
% arc_is_counter_clockwise     = arc_parameters(1,7);

% Find the change in angle of the arc
% arc1_start_unit_vector        = [cos(arc1_start_angle_in_radians) sin(arc1_start_angle_in_radians)];
% arc1_end_unit_vector          = [cos(arc1_end_angle_in_radians)   sin(arc1_end_angle_in_radians)  ];
% if 1==arc_is_counter_clockwise
%     arc_is_counter_clockwise = 1;
% else
%     arc_is_counter_clockwise = -1;
% end
% arc1_change_in_angle = fcn_geometry_findAngleUsing2PointsOnCircle([0 0],1, arc1_start_unit_vector, arc1_end_unit_vector, arc1_is_counter_clockwise);


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_base_point_xy           = segment_parameters(1,1:2);
segment_unit_tangent_vector     = [cos(segment_parameters(1,3)) sin(segment_parameters(1,3))];
% segment_s_start                 = segment_parameters(1,5);
% segment_s_end                   = segment_parameters(1,6);
% segment_angle = atan2(segment_unit_tangent_vector(2),segment_unit_tangent_vector(1));
segment_unit_ortho_vector = segment_unit_tangent_vector*[0 1; -1 0];
% segment_length                  = segment_s_end - segment_s_start;

% Calculate the distance from the arc's circle to the segment
vector_from_arc_center_to_segment_base_point = segment_base_point_xy - arc_center_xy;
distance_from_arc_center_to_segment = abs(sum(segment_unit_ortho_vector.*vector_from_arc_center_to_segment_base_point,2));

% Calculate the distance between the circles and the join point
space_between_arc_circle_and_segment = distance_from_arc_center_to_segment - arc_radius;

if length(threshold)==1 || isempty(threshold)
    transverse_threshold = threshold;
else
    transverse_threshold = threshold(2);
end

switch continuity_level
    case 0
        if space_between_arc_circle_and_segment<=0
            firstFitType = 'arc';
            firstFitType_parameters = arc_parameters;
            secondFitType = 'segment';
            secondFitType_parameters = segment_parameters;

            intersection_point_arc_to_segment = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, -1);
            if ~any(isnan(intersection_point_arc_to_segment),'all')
                flag_is_feasible = 1;
            else
                flag_is_feasible = 0;
            end
        else
            flag_is_feasible = 0;
        end
        feasibility_distance = space_between_arc_circle_and_segment;
        % BUG - SNB - 2024_05_26 - need to calculate feasibility if there is a threshold
        % given and there is NOT any intersections. If there is no
        % intersection, then need to check if the end of the line segment
        % if projected by the s-distance outward against the arc projected
        % s-distance along arc is within the transverse distance. This
        % isn't a trivial calculation - but it should be done.

    case 1
        % C1 connectivity is feasible if the arc could be tangent to the
        % segment, if allowed to move as much as the threshold
        if abs(space_between_arc_circle_and_segment)<=transverse_threshold
            flag_is_feasible = 1;
        else
            flag_is_feasible = 0;
        end
        feasibility_distance = abs(space_between_arc_circle_and_segment)-transverse_threshold;

    case 2
        % C2 connectivity is only feasible if the space between the circle
        % and line segment is positive
        open_space_distance = space_between_arc_circle_and_segment + transverse_threshold;
        if open_space_distance>=0
            flag_is_feasible = 1;
        else
            flag_is_feasible = 0;
        end
        feasibility_distance = -open_space_distance;

    otherwise
        warning('on','backtrace');
        warning('An error will be thrown at this point due to bad continuity number.');
        error('Alignments are not yet supported for continuities between segments and arcs with continuity type: %.0d',continuity_level);
end
closest_feasible_input2_parameters = [];
end % Ends fcn_INTERNAL_confirmSegmentToArc

%% fcn_INTERNAL_confirmArcToArc
function  [flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_INTERNAL_confirmArcToArc(arc1_parameters, arc2_parameters, continuity_level, threshold)

% Fill defaults
flag_is_feasible = 0;
feasibility_distance = [];
closest_feasible_arc2_parameters = [];


% Calculate needed values from parameter sets
% Calculate needed values from parameter sets
% Get the arc fit details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_center_xy                = arc1_parameters(1,1:2);
arc1_radius                   = arc1_parameters(1,3);
% arc1_start_angle_in_radians   = arc1_parameters(1,4);
% arc1_end_angle_in_radians     = arc1_parameters(1,5);
% arc1_start_xy                 = arc1_center_xy + arc1_radius*[cos(arc1_start_angle_in_radians) sin(arc1_start_angle_in_radians)];
% arc1_end_xy                   = arc1_center_xy + arc1_radius*[cos(arc1_end_angle_in_radians) sin(arc1_end_angle_in_radians)];
% arc1_is_counter_clockwise     = arc1_parameters(1,7);

% Find the change in angle of the arc
% arc1_start_unit_vector        = [cos(arc1_start_angle_in_radians) sin(arc1_start_angle_in_radians)];
% arc1_end_unit_vector          = [cos(arc1_end_angle_in_radians)   sin(arc1_end_angle_in_radians)  ];
% if 1==arc1_is_counter_clockwise
%     arc1_is_counter_clockwise = 1;
% else
%     arc1_is_counter_clockwise = -1;
% end
% % arc1_change_in_angle = fcn_geometry_findAngleUsing2PointsOnCircle([0 0],1, arc1_start_unit_vector, arc1_end_unit_vector, arc1_is_counter_clockwise);


% Get the arc fit details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_center_xy                = arc2_parameters(1,1:2);
arc2_radius                   = arc2_parameters(1,3);
% arc2_start_angle_in_radians   = arc2_parameters(1,4);
% arc2_end_angle_in_radians     = arc2_parameters(1,5);
% arc2_start_xy                 = arc2_center_xy + arc2_radius*[cos(arc2_start_angle_in_radians) sin(arc2_start_angle_in_radians)];
% arc2_end_xy                   = arc2_center_xy + arc2_radius*[cos(arc2_end_angle_in_radians) sin(arc2_end_angle_in_radians)];
arc2_is_counter_clockwise     = arc2_parameters(1,7);

% Find the change in angle of the arc
% arc2_start_unit_vector        = [cos(arc2_start_angle_in_radians) sin(arc2_start_angle_in_radians)];
% arc2_end_unit_vector          = [cos(arc2_end_angle_in_radians)   sin(arc2_end_angle_in_radians)  ];
if 1==arc2_is_counter_clockwise
    arc2_is_counter_clockwise = 1;
else
    arc2_is_counter_clockwise = -1;
end
% arc2_change_in_angle = fcn_geometry_findAngleUsing2PointsOnCircle([0 0],1, arc2_start_unit_vector, arc2_end_unit_vector, arc2_is_counter_clockwise);

% Calculate the distance between the circles and the join point
shift_xy = arc1_center_xy - [0 arc1_radius]; % The function below assumes that arc1 is centered at [0 r1], so we need to shift everything to this location
space_between_circles = fcn_geometry_gapCircleToCircle(arc1_radius, arc2_radius, arc2_center_xy-shift_xy, arc2_is_counter_clockwise,-1);

if length(threshold)==1 || isempty(threshold)
    transverse_threshold = threshold;
else
    transverse_threshold = threshold(2);
end

switch continuity_level
    case 0
        % C0 continuity requires an intersection
        feasibility_distance = space_between_circles;
        if space_between_circles<=0
            firstFitType = 'arc';
            firstFitType_parameters = arc1_parameters;
            secondFitType = 'arc';
            secondFitType_parameters = arc2_parameters;

            intersection_point_arc_to_segment = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, -1);
            if ~any(isnan(intersection_point_arc_to_segment),'all')
                flag_is_feasible = 1;
            else
                flag_is_feasible = 0;
                feasibility_distance = -space_between_circles; % THIS IS WRONG - need to calculate the result.
            end
        else
            flag_is_feasible = 0;
        end

        % BUG - SNB - 2024_05_26 - need to calculate feasibility if there is a threshold
        % given and there is NOT any intersections. If there is no
        % intersection, then need to check if the end of the line segment
        % if projected by the s-distance outward against the arc projected
        % s-distance along arc is within the transverse distance. This
        % isn't a trivial calculation - but it should be done.

    case 1
        % For C1 continuity of arc1 to arc2, the closest point of the desired
        % joint to arc1 and arc2 is always going to be the point where each arc
        % touches the tangent line connecting them. Since arc2 is the only one
        % that can be moved, the connecting point for arc1 is simply the location where
        % arc1 is tangent, which by construction is [0 0]. The connection
        % point for arc2 is where the arc2's circle touches the x-axis,
        % which by construction is the radius distance below the center of
        % the circle

        if (space_between_circles>0) && (1==arc2_is_counter_clockwise)
            % In this case, one circle is inside the other and an external
            % tangent is requested - and this is not possible. So need to
            % check to see if arc2 can be moved to create an external
            % tangent.

            if abs(space_between_circles)<transverse_threshold
                % Yes, arc2 can be moved enough to be tangent. 
                flag_is_feasible = 1;
            else
                flag_is_feasible = 0;
            end
            feasibility_distance = abs(space_between_circles)-transverse_threshold;

        elseif (space_between_circles<=0) && (1==arc2_is_counter_clockwise)
            % In this case, one circle is OUTSIDE the other and an external
            % tangent is requested - and this is always possible.
            
            flag_is_feasible = 1;
            feasibility_distance = space_between_circles-transverse_threshold;
            
        elseif (space_between_circles>0)  && (1~=arc2_is_counter_clockwise)
            % The circles do not overlap, and an outside tangent is
            % requested. This is always possible.
            
            flag_is_feasible = 1;
            feasibility_distance = -space_between_circles-transverse_threshold;
            

        elseif (space_between_circles<=0) && (1~=arc2_is_counter_clockwise)
            % The circles overlap, and an outside tangent is requested.
            % This isn't possible unless arc2 can be moved

            if abs(space_between_circles)<transverse_threshold
                % Yes, arc2 can be moved enough to be tangent.
                flag_is_feasible = 1;
            else
                % Not possible to shift enough to allow connection
                flag_is_feasible = 0;
            end
            feasibility_distance = -transverse_threshold+space_between_circles;

        else
            error('Unknown error encountered - it should not be possible to enter this case!');            
        end

        % This part is not yet done
        closest_feasible_arc2_parameters = nan;
        warning('This section is not completely finished.');
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

        in_boundary_margin = [];
        [flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleArcToArc(arc1_parameters, arc2_parameters,(threshold), (in_boundary_margin), (-1));

    otherwise
        warning('on','backtrace');
        warning('An error will be thrown at this point due to bad continuity number.');
        error('Alignments are not yet supported for continuities between segments and arcs with continuity type: %.0d',continuity_level);
end
end % Ends fcn_INTERNAL_confirmArcToArc

