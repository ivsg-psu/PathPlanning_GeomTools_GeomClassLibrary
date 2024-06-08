function [flag_is_feasible, feasibility_distance, closest_feasible_arc2_parameters] = fcn_geometry_isC2FeasibleArcToArc(arc1_parameters, arc2_parameters, varargin)
%% fcn_geometry_isC2FeasibleArcToArc
% Given two arc geometries, checks if the arcs can be joined with C2
% continuity. If not, it finds the closest C2 feasible arc2 given arc1.
%
% Format:
% [flag_is_feasible, feasibility_distance] = ...
%    arc1_parameters, arc2_parameters, ...
%    (threshold), (in_boundary_margin), (fig_num))
%
% INPUTS:
%
%      arc1_parameters: a vector of arc1 parameters consistent with arc or
%      circle type geometries. See fcn_geometry_fillEmptyDomainStructure
%      for details
%
%      arc2_parameters: a vector of arc2 parameters consistent with arc or
%      circle type geometries. See fcn_geometry_fillEmptyDomainStructure
%      for details
%
%      (OPTIONAL INPUTS)
%
%      threshold: the offset, in meters, between the arcs such that this
%      offset can be removed by shifting. If the needed feasible offset is
%      larger than this value, then the flag_is_feasible is set to 0. If
%      threshold is entered as a 2x1 or 1x2, then this specifies the
%      threshold in St coordinates, e.g. first in the station direction,
%      and then in the transverse direction. For example, an entry of [3
%      0.2] would have 3 meters threshold in the station direction, but 0.2
%      meters threshold in the transverse direction. Note that only the
%      transverse distances are used to determine feasibility. The default
%      threshold is 0
%
%      in_boundary_margin: the distance, in meters, that the offset must be
%      moved within the boundary for feasibility. This margin is added
%      because feasible solutions that exact exactly on the boundary
%      (within the tolerance) may still yield solutions that are
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
%      closest_feasible_arc2_parameters: the closest parameters that give
%      feasible continuity between arc2 and arc1.
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_isC2FeasibleArcToArc
% for a full test suite.
%
% This function was written on 2024_05_26 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_05_26 - S Brennan
% -- wrote the code
% 2024_06_08 - S Brennan
% -- added in_boundary_margin input option

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

% Extract needed values from parameter sets
% Get the details from arc1 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc1_center_xy                = arc1_parameters(1,1:2);
r1                            = arc1_parameters(1,3);

% Get the details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc2_center_xy                = arc2_parameters(1,1:2);
r2                            = arc2_parameters(1,3);

d12 = real(sum((arc1_center_xy-arc2_center_xy).^2,2).^0.5);

query_point = [r2 d12];

% Method: 
% An arc2 is C2 feasible connected to arc1 if it is one of 3 cases:
% 1) both circles are oriented the same direction, and the circle of arc2
% is completely encircled by the circle for arc1
% 2) both circles are oriented the same direction, and the circle of arc2
% completely encircles the circle for arc1
% 3) both circles are oriented in opposite direction, and the circle of arc2
% exists entirely outside the circle for arc1
%
% Each of these cases produces an inequality requirement for the radius of
% arc2 (r2) and the distance from the center of the circles for arc1
% and arc2 (d12), given a radius of arc1 (r1).
%
% Specifically:
%
% CASE 1: This only occurs if (r2 < r1), both CW or both CCW
%
% For case 1, when r2 is smaller than r1, the circle2 is inside circle1
% only if:
%  
% r2 < (r1 - d12) 
%
% this can be rewritten with d12 as the dependent variable (y-axis):
%
%  d12 < r1 - r2
% 
% If this is plotted with d12 on the y-axis and r2 on the x-axis, then it
% creates a line with slope of -1 and r1 as the y-intercept. The line
% intercepts the x-axis at the location where r2=r1 and d12 is zero, as
% expected. 
%
% To determine the distance of a query point (r2,d12) to this line, a
% vector can be created to the query point from either (0, r1) or (r1,0).
% The dot product of this vector with the projection normal from the
% boundary, in the [1 1] direction, is the distance to the boundary.
%
% CASE 2: This only occurs if (r2 > r1), both CW or both CCW
%
% For case 2, when r2 is larger than r1, the circle2 is encircling circle1
% only if:
%  
% r2 > (r1 + d12) 
%
% this can be rewritten with d12 as the dependent variable (y-axis):
%
%  d12 < r2 - r1
% 
% If this is plotted with d12 on the y-axis and r2 on the x-axis, then it
% creates a line with slope of 1 and -r1 as the y-intercept. The line
% intercepts the x-axis at the location where r2=r1 and d12 is zero, as
% expected. The feasible range exists only for d12 greater than zero, e.g.
% within a right triangle whose corner starts at (r1, 0).
%
% To determine the distance of a query point (r2,d12) to this line, a
% vector can be created to the query point from either (0, -r1) or (r1,0).
% The dot product of this vector with the projection normal from the
% boundary, in the [-1 1] direction, is the distance to the boundary.
%
% CASE 3: one arc is CW the other is CCW
%
% For case 3, the only feasible solution is when the circle for arc2 is
% completely outside the cirlce for arc1. This occurs only if:
%  
% d12 > (r1 + r2) 
%
% If this is plotted with d12 on the y-axis and r2 on the x-axis, then it
% creates a line with slope of 1 and r1 as the intercept. The line 
% intercepts the x-axis at the location where r2=-r1 and d12 is zero which
% is not physically possible (negative radius). The feasible range exists
% only for d12 greater than zero, e.g. within a right triangle whose corner
% starts at (0, r1).
%
% To determine the distance of a query point (r2,d12) to this line, a
% vector can be created to the query point from either (0, r1) or (-r1,0).
% The dot product of this vector with the projection normal from the
% boundary, in the [1 -1] direction, is the distance to the boundary.
%
% These inequalities therefore give the method to find both feasibility and
% the closest location: project the test point (r2, d12) into each of the 3
% directions. Only 1 of them can be negative, and if this is the case, then
% the C2 connectivity is feaible. If none are negative, then take the
% minimum distance solution, project this distance into the feasible space
% to find the feasible solution. As well, compare this projection distance
% to the threshold to see if the projection would be allowed to determine
% feasibility.

cases_test_vectors = [...
    query_point - [r1, 0];
    query_point - [r1, 0];
    query_point - [0, r1];
    ];

cases_projection_vectors = [...
    1 1;
    -1 1;
    1 -1];

cases_unit_projection_vectors = fcn_geometry_calcUnitVector(cases_projection_vectors);
distances_by_case = sum(cases_unit_projection_vectors.*cases_test_vectors,2);
[closest_distance,case_number] = min(distances_by_case);

if closest_distance<=0
    % In this case, the join is fully feasible without correction
    flag_is_feasible = 1;
    feasibility_distance = closest_distance;
    closest_feasible_arc2_parameters = arc2_parameters;
    corrected_parameters_r2_d12 = [r2 d12];
else
    % Some correction is necessary for the join to be feasible

    % Find closest parameters
    corrected_parameters_r2_d12 = query_point - cases_unit_projection_vectors(case_number,:)*(closest_distance + in_boundary_margin);
    closest_feasible_arc2_parameters = arc2_parameters;

    % The radius is the first value
    closest_feasible_arc2_parameters(1,3) = corrected_parameters_r2_d12(1,1);

    % The new center can be calculated by the center-to-center vector
    updated_d12 = corrected_parameters_r2_d12(1,2);
    vector_center1_to_center2 = arc2_center_xy - arc1_center_xy;
    unit_vector_center1_to_center2 = fcn_geometry_calcUnitVector(vector_center1_to_center2);
    updated_arc2_center_xy = arc1_center_xy + unit_vector_center1_to_center2*updated_d12;
    closest_feasible_arc2_parameters(1,1:2) = updated_arc2_center_xy;
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

    segment_length = [];
    format_string = sprintf(' ''-'',''Color'',[0.6 0 0],''LineWidth'',7 ');
    fcn_geometry_plotGeometry('circle', arc1_parameters,segment_length, format_string, (fig_num));
    format_string = sprintf(' ''-'',''Color'',[1 0 0],''LineWidth'',4 ');
    fcn_geometry_plotGeometry('circle', arc2_parameters,segment_length, format_string, (fig_num));

    if 0==flag_is_feasible
        title(sprintf('NOT C2 feasible, dist: %.2f',feasibility_distance));
        format_string = sprintf(' ''-'',''Color'',[0 1 0],''LineWidth'',4 ');
        fcn_geometry_plotGeometry('circle', closest_feasible_arc2_parameters,segment_length, format_string, (fig_num));
    else
        title(sprintf('C2 feasible, dist: %.2f',feasibility_distance));
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
    xlabel('r2 [meters]');
    ylabel('d12 [meters]');
    
    % Plot the cases
    plot([0 r1],[r1 0],'b-','LineWidth',3);
    plot([r1 2*r1],[0 r1],'m-','LineWidth',3);
    plot([0 r1],[r1 2*r1],'c-','LineWidth',3);
    plot(query_point(1,1),query_point(1,2),'r.','MarkerSize',30);
    plot(corrected_parameters_r2_d12(1,1),corrected_parameters_r2_d12(1,2),'g.','MarkerSize',30);
    
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
