function [flag_is_feasible, feasibility_distance, closest_feasible_arc_parameters] = ...
    fcn_geometry_isC2FeasibleLineToArc(segment_parameters, arc_parameters, varargin)
%% fcn_geometry_isC2FeasibleLineToArc
% Given a fixed line geometry and variable arc geometry, including a
% possible threshold to move the arc geometry, checks if the line and arc
% can be joined with C2 continuity. The method is to calculate the closest
% C2 feasible arc parameter set given the line, then use this to calculate
% the distance in parameter space to this closest parameter set, and then -
% if this distance is less than either zero (default) or a positive
% tolerance, calculates flag_is_feasible.
%
% Format:
%  [flag_is_feasible, feasibility_distance, closest_feasible_arc_parameters] = fcn_geometry_isC2FeasibleArcToArc(...
%    segment_parameters, arc_parameters, ...
%    (threshold), (in_boundary_margin), (fig_num))
%
% INPUTS:
%
%      segment_parameters: a vector of segment parameters consistent with
%            line or segment geometries, e.g.
%  
%            [base_point_x base_point_y angle_of_vector]. 
% 
%            See fcn_geometry_fillEmptyDomainStructure for details
%
%      arc_parameters: a vector of arc parameters consistent with arc or
%            circle geometries, e.g. 
%  
%            [circle_center_x circle_center_y radius]. 
% 
%            See fcn_geometry_fillEmptyDomainStructure for details
%
%      (OPTIONAL INPUTS)
%
%      threshold: the offset, in meters, that either the radius or circle
%      center can be moved to create C2 feasibility. If the
%      needed feasible offset is larger than this value, then the
%      flag_is_feasible is set to 0. If threshold is entered as a 2x1 or
%      1x2, then this specifies the threshold in St coordinates, e.g. first
%      in the station direction, and then in the transverse direction. For
%      example, an entry of [3 0.2] would have 3 meters threshold in the
%      station direction, but 0.2 meters threshold in the transverse
%      direction. Note that only the transverse distances are used to
%      determine feasibility in shifting. The default threshold is 0.
%
%      in_boundary_margin: the distance, in meters, that the corrected
%      parameters must be moved within the boundary for feasibility. This
%      margin is added because feasible parameters that lie exactly on the
%      boundary of feasibility may still yield solutions that are
%      technically not feasible spiral connections. Default is 0.
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      flag_is_feasible: a value of 1 if C2 continuity between arc1 and
%      arc2 is feasible for the given tolerance, 0 otherwise. 
%
%      feasibility_distance: the distance between arc1 and arc2 that gives
%      feasibile C2 continuity.
%
%      closest_feasible_arc_parameters: the closest parameters that give
%      feasible continuity between arc2 and arc1.
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_isC2FeasibleLineToArc
% for a full test suite.
%
% This function was written on 2024_06_24 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_06_24 - S Brennan
% -- wrote the code

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
        narginchk(2,5);

    end
end


% Does user want to specify threshold?
threshold = 0;
if (3<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        threshold = temp;
    end
end

% Does user want to specify in_boundary_margin?
in_boundary_margin = 0;
if (4<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        in_boundary_margin = temp;
    end
end


% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if  (0==flag_max_speed) && (5<= nargin)
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


%% Solve for the closest boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Determine if the inputs are segments or lines, arcs or circles
% Segments have 4 parameters, lines have 3
if length(segment_parameters(1,:))==4
    flag_segmentInput_is_segment = 1;
else
    flag_segmentInput_is_segment = 0;
end

% Arcs have 7 parameters, circles have 3
if length(arc_parameters(1,:))==7
    flag_arc2_is_arc = 1;
else
    flag_arc2_is_arc = 0;
end

%%%% Do not need this section - it is a carry-over from ArcArc checking
%
% if 1==flag_segmentInput_is_segment && 1==flag_arc2_is_arc
%     flag_check_arcs = 1;
%     if segment_parameters(1,7) == arc_parameters(1,7)
%         flag_oriented_same_direction = 1;
%     else
%         flag_oriented_same_direction = 0;
%     end
% else
%     flag_check_arcs = 0;
% end

%% Extract needed values from parameter sets
% Get the details from arc1 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_base_xy                = segment_parameters(1,1:2);
segment_angle                  = segment_parameters(1,3);

% Get the details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy                = arc_parameters(1,1:2);
r2                            = arc_parameters(1,3);

% Find distance, d12 from segment to the circle by using dot product of
% orthogonal projection
segment_unit_tangent_vector = [cos(segment_angle) sin(segment_angle)];
segment_unit_orthogo_vector = segment_unit_tangent_vector*[0 1; -1 0];
vector_from_segment_base_to_circle_center = arc_center_xy - segment_base_xy;

d12 = real(sum(vector_from_segment_base_to_circle_center.*segment_unit_orthogo_vector,2));

query_point = [r2 d12];

% Method: 
% An circle is C2 feasible connected to a fixed line if it does not
% intersect the line.
%
% This requirement produces an inequality requirement for the radius of the
% arc (r2) and the orthogonal distance from the line 
% and arc2 (d12).
%
% Specifically:
% We require (r2 < d12). 
%
% If this is plotted with d12 on the y-axis and r2 on the x-axis, then it
% creates a 45-degree line with slope of 1 and the origin as the
% y-intercept. 
%
% To determine the distance of a query point (r2,d12) to this line, a
% vector can be created to the query point from the origin.
% The dot product of this vector with the projection normal from the
% boundary, in the [-1 1] direction, is the distance from the boundary to
% the query point.

URHERE

cases_test_vectors = ...
    query_point - [0, 0];

if flag_check_arcs
    if flag_oriented_same_direction
        cases_unit_projection_vectors = [...
            1 1;
            -1 1;
            nan nan]/(2^0.5);
    else
        cases_unit_projection_vectors = [...
            nan nan;
            nan nan;
            1 -1]/(2^0.5);
    end
else
    cases_unit_projection_vectors = [...
        1 1;
        -1 1;
        1 -1]/(2^0.5);
end

distances_by_case = sum(cases_unit_projection_vectors.*cases_test_vectors,2);
[closest_distance,case_number] = min(distances_by_case);

if closest_distance<=0
    % In this case, the join is fully feasible without correction
    flag_is_feasible = 1;
    feasibility_distance = closest_distance;
    closest_feasible_arc_parameters = arc_parameters;
    corrected_parameters_r2_d12 = [r2 d12];
else
    % Some correction is necessary for the join to be feasible

    % Find closest parameters
    corrected_parameters_r2_d12 = query_point - cases_unit_projection_vectors(case_number,:)*(closest_distance + in_boundary_margin);
    closest_feasible_arc_parameters = arc_parameters;

    % The radius is the first value
    closest_feasible_arc_parameters(1,3) = corrected_parameters_r2_d12(1,1);

    % The new center can be calculated by the center-to-center vector
    updated_d12 = corrected_parameters_r2_d12(1,2);
    vector_center1_to_center2 = arc_center_xy - segment_base_xy;
    unit_vector_center1_to_center2 = fcn_geometry_calcUnitVector(vector_center1_to_center2);
    updated_arc2_center_xy = segment_base_xy + unit_vector_center1_to_center2*updated_d12;
    closest_feasible_arc_parameters(1,1:2) = updated_arc2_center_xy;
    feasibility_distance = closest_distance;

    if closest_distance <= threshold
        flag_is_feasible = 1;
    else
        flag_is_feasible = 0;
    end
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

    % Plot the results in the given figure number
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end
    % end

    % Plot the input geometries 
    subplot(1,2,1);
    hold on;
    grid on;
    axis equal;
    xlabel('X [meters]');
    ylabel('Y [meters]');

    % Plot the inputs
    segment_length = [];
    format1_string = sprintf(' ''-'',''Color'',[0.6 0 0],''LineWidth'',7 ');
    
    if 1==flag_segmentInput_is_segment
        fcn_geometry_plotGeometry('arc', segment_parameters,segment_length, format1_string, (fig_num));
    else
        fcn_geometry_plotGeometry('circle', segment_parameters,segment_length, format1_string, (fig_num));
    end

    format2_string = sprintf(' ''-'',''Color'',[1 0 0],''LineWidth'',4 ');
    if 1==flag_arc2_is_arc
        fcn_geometry_plotGeometry('arc', arc_parameters,segment_length, format2_string, (fig_num));
    else
        fcn_geometry_plotGeometry('circle', arc_parameters,segment_length, format2_string, (fig_num));
    end

    % Plot the corrected solution?
    if 0==flag_is_feasible
        format_string = sprintf(' ''-'',''Color'',[0 1 0],''LineWidth'',4 ');
        if flag_check_arcs
            fcn_geometry_plotGeometry('arc', closest_feasible_arc_parameters,segment_length, format_string, (fig_num));
        else
            fcn_geometry_plotGeometry('circle', closest_feasible_arc_parameters,segment_length, format_string, (fig_num));
        end
    end

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end

    % Plot the aligned geometries on the left
    subplot(1,2,2);
    hold on;
    grid on;
    axis equal;
    xlabel('r2/r1 [meters]');
    ylabel('d12/r1 [meters]');
    
    % Plot the cases
    plot([0 1],[1 0],'b-','LineWidth',3);
    plot([1 2],[0 1],'m-','LineWidth',3);
    plot([0 1],[1 2],'c-','LineWidth',3);
    plot(query_point(1,1)/segment_angle,query_point(1,2)/segment_angle,'r.','MarkerSize',30);
    plot(corrected_parameters_r2_d12(1,1)/segment_angle,corrected_parameters_r2_d12(1,2)/segment_angle,'g.','MarkerSize',30);
    
    legend('Case1','Case2','Case3','Query','Fixed')

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end

    % Add a title
    if 0==flag_is_feasible
        sgtitle(sprintf('NOT C2 feasible, dist: %.2f',feasibility_distance));
    else
        sgtitle(sprintf('C2 feasible, dist: %.2f',feasibility_distance));
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
