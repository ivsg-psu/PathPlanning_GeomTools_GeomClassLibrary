function intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, varargin) 
%% fcn_geometry_intersectGeom
% This function finds the intersection point(s) of two geometries, such as
% line-arc, line-circle, line segment-arc, arc-arc, arc-circle, arc-line,
% etc.
% 
% For segments, arcs, and spirals that may have multiple intersection
% points, the function always returns the point that is nearest to the
% start of the first segment, arc, or spiral, with "nearest" meaning in
% station distance, not geometric position. In other words, returns the
% intersection that is encountered first when traversing the "from"
% geometry starting at its "start" position.
%
% For circles and lines that are the "starting" geometry, the intersection
% point is returned that is the first point encountered in the 2nd geometry
% if the 2nd geometry is a segment, arc, or spiral.
%
% For geometries where both the first and second geometries have no start
% or end, for example circles intersecting with circles or lines, both
% intersection points are returned if there are more than 1 intersections.
%
% FORMAT: 
%
% intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num) 
% 
% INPUTS:
%
%      firstFitType: The string input of the first geometry. Can be one of the
%      following strings:
%   
%           'arc', 'line', 'circle', 'line segment' or 'segment'
%           (Pending: 'spiral')
%
%      firstFitType_parameters: The parameters of the first geometry as a
%      [1xN] vector.  See fcn_geometry_fillEmptyDomainStructure for details
%      on how to fill the vector for different geometry types.
%
%      secondFitType: The string input of the second geometry that uses the
%      same string inputs as firstFitType.
%
%      secondFitType_parameters: The parameters of the second geometry,
%      using the same format as firstFitType_parameters.
%
%
% (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
% intersection_points: This is a nx2 intersection points matrix of the
% given geometries
% 
%
% DEPENDENCIES:
%
%   fcn_geometry_arcDirectionFrom3Points
%   fcn_geometry_arcAngleFrom3Points
%   fcn_geometry_plotGeometry
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_intersectGeom
% for a full test suite.
%
% This function was written on 2024_05_02 by Aneesh Batchu with minor
% revisions by Sean Brennan
% Questions or comments? abb6486@psu.edu or sbrennan@psu.edu

% Revision History
% 2024_05_02 - Aneesh Batchu
% -- wrote the code 
% 2024_05_06 - Aneesh Batchu 
% -- Fixed BUGS in Arc-Arc, line-arc, arc-line, line segment-arc
% intersection cases. Added a few conditional statements to remove NaNs
% from the potential_intersection_points to compute the true "intersection
% points"
% 2024_05_10 - Sean Brennan
% -- added some more comments
% -- added 'segment' alongside 'line segment' option
% -- added TO-DO list
% -- added backtrace on to each error call, to allow fast debugging later
% -- changed plot style to better see difference between 1st and 2nd input
% -- functionalized the arc inputs codes (need to do the rest)
% -- added arc to segment case 
% -- reordered functions to make cross-dependencies easier, for example
% arc-to-segment calls circle-to-line, then arc-to-line, before doing
% arc-to-segment check. That way, if bugs are found/fixed in
% circle-to-line, they automatically get fixed in arc-to-segment.
% 2024_05_13 - Aneesh Batchu
% -- Fixed a bug in the arc-segment case where the arc direction (when
% arc_parameters(1,7) is undefined (=0)). The function
% "fcn_INTERNAL_findPointsInArc" calculates the arc angle using arc_start,
% arc_intersection, and arc_end, and outputs the intersection point that is
% closest to the arc.

% TO-DO:
% 2025_05_10 - added by Sean Brennan
% -- add arc to circle calculations, and many others that are missing. See
% test script which reveals many missing cases
% -- main code is hard to read because it is just a long (LONG) combination
% of cases. Need to functionalize the code. I started fixing this on arc, but
% organization is needed for all the other cases.

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
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(4,5);

    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (5<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

%% Main Code starts from here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch lower(firstFitType)
    case 'circle'
        intersection_points = fcn_INTERNAL_intersectCircleGeomtries(secondFitType, firstFitType_parameters, secondFitType_parameters); 
    case 'arc'
        intersection_points = fcn_INTERNAL_intersectArcGeomtries(secondFitType, firstFitType_parameters, secondFitType_parameters); 
    case 'line'
        intersection_points = fcn_INTERNAL_intersectLineGeometries(secondFitType, firstFitType_parameters, secondFitType_parameters);
    case {'segment'}
        intersection_points = fcn_INTERNAL_intersectSegmentGeomtries(secondFitType, firstFitType_parameters, secondFitType_parameters);          
    case 'spiral'
        warning('on','backtrace');
        warning('An error will be thrown at this point due to missing code.');
        error('fcn_geometry_intersectionGeom case is not yet ready for any spiral start case');
    otherwise
        warning('on','backtrace');
        warning('An error will be thrown at this point due to missing code.');
        error('fcn_geometry_intersectionGeom is not yet ready for curves from fit type: %s',firstFitType);

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

    title('Testing Intersection of Geometries');
    xlabel('X [meters]');
    ylabel('Y [meters]')

    % Plot the inputs
    fcn_geometry_plotGeometry(lower(firstFitType),firstFitType_parameters,[],sprintf(' ''LineWidth'',4 '));
    fcn_geometry_plotGeometry(lower(secondFitType),secondFitType_parameters,[],sprintf(' ''LineWidth'',2 '));
    plot(intersection_points(:,1),intersection_points(:,2),'k.','MarkerSize',20);

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end


end

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
%% fcn_INTERNAL_intersectCircleGeomtries
function  intersection_points = fcn_INTERNAL_intersectCircleGeomtries(secondFitType, first_circle_parameters, secondFitType_parameters)
switch lower(secondFitType)
    case 'circle'
        intersection_points = fcn_INTERNAL_intersectCircleCircle(first_circle_parameters, secondFitType_parameters);
    case 'arc'
        intersection_points = fcn_INTERNAL_intersectCircleArc(first_circle_parameters, secondFitType_parameters);
    case 'line'
        intersection_points = fcn_INTERNAL_intersectLineCircle(secondFitType_parameters, first_circle_parameters);
    case 'segment'        
        intersection_points = fcn_INTERNAL_intersectCircleSegment(first_circle_parameters, secondFitType_parameters); 
    case 'spiral'
        warning('on','backtrace');
        warning('An error will be thrown at this point due to missing code.');
        error('fcn_geometry_intersectionGeom case is not yet ready for circle-spiral case');
    otherwise
        warning('on','backtrace');
        warning('An error will be thrown at this point due to missing code.');
        error('fcn_geometry_intersectionGeom is not yet ready for curves from fit type: %s',firstFitType);
end
end % Ends fcn_INTERNAL_intersectCircleGeomtries

%% fcn_INTERNAL_intersectCircleCircle
function intersection_points = fcn_INTERNAL_intersectCircleCircle(firstFitType_parameters, secondFitType_parameters)
% Calculate needed values from parameter sets

% Get the circle1 fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle1_center_xy                = firstFitType_parameters(1,1:2);
circle1_radius                   = firstFitType_parameters(1,3);

% Get the circle2 fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle2_center_xy                = secondFitType_parameters(1,1:2);
circle2_radius                   = secondFitType_parameters(1,3);

% Are the circles the same to 8 decimal places?
if isequal(round(circle1_center_xy,8),round(circle2_center_xy,8)) && isequal(round(circle1_radius,8),round(circle2_radius,8))
    intersection_points = [inf inf];
else
    % Use MATLAB's circcirc algorithm to find intersections between two circles
    [xout,yout] = circcirc(circle1_center_xy(1,1),circle1_center_xy(1,2),circle1_radius,circle2_center_xy(1,1),circle2_center_xy(1,2),circle2_radius);

    if ~isnan(xout)
        % intersection points were found! To be an intersection, the point must
        % be on both arc1 and arc2

        % Which point(s) to keep?
        intersection_points = [xout', yout'];

        % Is the point duplicated?
        if isequal(intersection_points(1,:),intersection_points(2,:))
            intersection_points = intersection_points(1,:);
        end
    else
        intersection_points = [nan nan];
    end
end
end % ends fcn_INTERNAL_intersectCircleCircle

%% fcn_INTERNAL_intersectCircleArc
function intersection_points = fcn_INTERNAL_intersectCircleArc(circle_parameters, arc_parameters)

intersection_points = fcn_INTERNAL_intersectCircleCircle(circle_parameters, arc_parameters);

if ~any(isnan(intersection_points))
    % Keep only the closet point on the arc. Returns nan if there are none
    [~, intersection_points] = fcn_INTERNAL_findPointsInArc(arc_parameters,intersection_points);
end

end % ends fcn_INTERNAL_intersectCircleArc

%% fcn_INTERNAL_findPointsInArc
function [points_in_arc, closest_point_to_start_of_arc] = fcn_INTERNAL_findPointsInArc(arc_parameters,points_to_test)
% Get the arc fit details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy                = arc_parameters(1,1:2);
arc_radius                   = arc_parameters(1,3);
arc_start_angle_in_radians   = arc_parameters(1,4);
arc_end_angle_in_radians     = arc_parameters(1,5);
arc_start_xy                 = arc_center_xy + arc_radius*[cos(arc_start_angle_in_radians) sin(arc_start_angle_in_radians)];
arc_end_xy                   = arc_center_xy + arc_radius*[cos(arc_end_angle_in_radians) sin(arc_end_angle_in_radians)];
arc_is_counter_clockwise     = arc_parameters(1,7);

% Are the intersections within the arc range that we were given? To
% check this, we use the three points on each arc - the start, the
% intersection, and the end to calculate the arc direction. We then
% check to see if it is the same as the given direction for that arc -
% if it is, the point is on the arc. The way we search is to initialize
% the potential arc intersection points to the circle intersections,
% and remove any of the arc intersection points that are not on both of
% the arcs.
points_in_arc = points_to_test;

if ~isequal(arc_is_counter_clockwise,0)
    for ith_row = 1:length(points_to_test(:,1))
        % Check arc
        intersection_is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(arc_start_xy, points_to_test(ith_row,:), arc_end_xy,-1);
        if arc_is_counter_clockwise ~= intersection_is_counterClockwise
            points_in_arc(ith_row,:) = [nan nan];
        end
    end
else
    arc_angle_in_radians_1_to_3 = zeros(length(points_to_test(:,1)),1);
    for ith_row = 1:length(points_to_test(:,1))
        [~, arc_angle_in_radians_1_to_3(ith_row,1), ~, ~, ~] = fcn_geometry_arcAngleFrom3Points(arc_start_xy, points_to_test(ith_row,:), arc_end_xy, -1);
    end
    if length(points_to_test(:,1))>1
        if abs(arc_angle_in_radians_1_to_3(1)) < abs(arc_angle_in_radians_1_to_3(2))
            points_in_arc(2,:) = [nan nan];
        else
            points_in_arc(1,:) = [nan nan];
        end
    end
end

% Remove the nan rows
points_in_arc = points_in_arc(~isnan(points_in_arc(:,1)),:);

% Make sure something remains
if isempty(points_in_arc)
    points_in_arc = [nan nan];
end

% Remove repeats
points_in_arc = unique(points_in_arc,'rows');

% At this point, there may be 2 intersection points. Are there 2?
% NOTE: the points may contain nan values, and if so, then nan is
% returned
closest_point_to_start_of_arc = fcn_INTERNAL_findPointClosestToArcStart(arc_parameters, points_in_arc);

end % Ends fcn_INTERNAL_findPointsInArc

%% fcn_INTERNAL_findPointClosestToArcStart
function closest_point = fcn_INTERNAL_findPointClosestToArcStart(arc_parameters, points_to_test)

% Get the arc fit details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy                = arc_parameters(1,1:2);
arc_radius                   = arc_parameters(1,3);
arc_start_angle_in_radians   = arc_parameters(1,4);
arc_end_angle_in_radians     = arc_parameters(1,5);
arc_start_xy                 = arc_center_xy + arc_radius*[cos(arc_start_angle_in_radians) sin(arc_start_angle_in_radians)];
arc_end_xy                   = arc_center_xy + arc_radius*[cos(arc_end_angle_in_radians) sin(arc_end_angle_in_radians)];

N_points = length(points_to_test(:,1));
arc_angles = nan(N_points,1);

% Calculate arc angles for each point
for ith_point = 1:N_points
    if ~any(isnan(points_to_test(ith_point,:)))
        arc_angles(ith_point,1)  = fcn_geometry_arcAngleFrom3Points(arc_start_xy, points_to_test(ith_point,:), arc_end_xy,(-1));
    end
end

[~, min_point_index] = min(arc_angles);

closest_point = points_to_test(min_point_index,:);
end % Ends fcn_INTERNAL_findPointClosestToArcStart

%% fcn_INTERNAL_intersectCircleSegment
function intersection_points = fcn_INTERNAL_intersectCircleSegment(circle_parameters, segment_parameters)

intersection_points = fcn_INTERNAL_intersectLineCircle(segment_parameters, circle_parameters);

if ~all(isnan(intersection_points))
    % Keep only the closest point on the segment. Returns nan if there are none
    [~, intersection_points] = fcn_INTERNAL_findPointsInSegment(segment_parameters,intersection_points);

end
end


%% fcn_INTERNAL_intersectArcGeomtries
function  intersection_points = fcn_INTERNAL_intersectArcGeomtries(secondFitType, first_arc_parameters, secondFitType_parameters)
switch lower(secondFitType)
    case 'circle'
        intersection_points = fcn_INTERNAL_intersectArcCircle(first_arc_parameters, secondFitType_parameters);
    case 'arc'
        intersection_points = fcn_INTERNAL_intersectArcArc(first_arc_parameters, secondFitType_parameters);
    case 'line'
        intersection_points = fcn_INTERNAL_intersectArcLine(first_arc_parameters, secondFitType_parameters);
    case 'segment'        
        intersection_points = fcn_INTERNAL_intersectArcSegment(first_arc_parameters,secondFitType_parameters);
    case 'spiral'
        warning('on','backtrace');
        warning('An error will be thrown at this point due to missing code.');
        error('fcn_geometry_intersectionGeom case is not yet ready for spiral case');
    otherwise
        warning('on','backtrace');
        warning('An error will be thrown at this point due to missing code.');
        error('fcn_geometry_intersectionGeom is not yet ready for curves from fit type: %s',firstFitType);
end
end % Ends fcn_INTERNAL_intersectArcGeomtries


%% fcn_INTERNAL_intersectArcCircle
function intersection_points = fcn_INTERNAL_intersectArcCircle(firstArc_parameters, circle_parameters)

% First, calculate circle to circle intersections
intersection_points = fcn_INTERNAL_intersectCircleCircle(firstArc_parameters, circle_parameters);

if ~all(isnan(intersection_points))
    % Keep only the closest point to start of the arc
    [~, intersection_points] = fcn_INTERNAL_findPointsInArc(firstArc_parameters,intersection_points);
end
end % Ends fcn_INTERNAL_intersectArcCircle

%% fcn_INTERNAL_intersectArcArc
function intersection_points = fcn_INTERNAL_intersectArcArc(firstArc_parameters, secondArc_parameters)

% First, calculate circle to circle intersections
intersection_points = fcn_INTERNAL_intersectCircleCircle(firstArc_parameters, secondArc_parameters);

if ~all(isnan(intersection_points))
    % Keep only the points on the second arc
    intersection_points = fcn_INTERNAL_findPointsInArc(secondArc_parameters,intersection_points);

    if ~all(isnan(intersection_points))
        % Keep only the closest point to start of the first arc
        [~, intersection_points] = fcn_INTERNAL_findPointsInArc(firstArc_parameters,intersection_points);
    end
end
end % Ends fcn_INTERNAL_intersectArcArc

%% fcn_INTERNAL_intersectArcLine
function intersection_points = fcn_INTERNAL_intersectArcLine(arc_parameters, line_parameters)

% First, calculate circle to line intersections
intersection_points = fcn_INTERNAL_intersectLineCircle(line_parameters, arc_parameters);

if ~all(isnan(intersection_points))

    % Keep only the closest point on the arc
    [~, intersection_points] = fcn_INTERNAL_findPointsInArc(arc_parameters,intersection_points);
end
end % Ends fcn_INTERNAL_intersectArcLine


%% fcn_INTERNAL_intersectArcSegment
function  intersection_points = fcn_INTERNAL_intersectArcSegment(arc_parameters, segment_parameters)

% First, calculate arc to line intersections
intersection_points = fcn_INTERNAL_intersectLineCircle(segment_parameters, arc_parameters);

if ~all(isnan(intersection_points))
    % Keep only the points on the segment
    intersection_points = fcn_INTERNAL_findPointsInSegment(segment_parameters,intersection_points);

    if ~all(isnan(intersection_points))
        % Keep only the closest point on the arc
        [~, intersection_points] = fcn_INTERNAL_findPointsInArc(arc_parameters,intersection_points);
    end
end
end % Ends fcn_INTERNAL_intersectArcSegment


%% fcn_INTERNAL_intersectLineGeomtries
function  intersection_points = fcn_INTERNAL_intersectLineGeometries(secondFitType, first_line_parameters, secondFitType_parameters)
switch lower(secondFitType)
    case 'circle'
        intersection_points = fcn_INTERNAL_intersectLineCircle(first_line_parameters, secondFitType_parameters);
    case 'arc'
        intersection_points = fcn_INTERNAL_intersectLineArc(first_line_parameters, secondFitType_parameters);
    case 'line'
        intersection_points = fcn_INTERNAL_intersectLineLine(first_line_parameters, secondFitType_parameters);
    case 'segment'        
        intersection_points = fcn_INTERNAL_intersectLineSegment(first_line_parameters,secondFitType_parameters);
    case 'spiral'
        warning('on','backtrace');
        warning('An error will be thrown at this point due to missing code.');
        error('fcn_geometry_intersectionGeom case is not yet ready for spiral case');
    otherwise
        warning('on','backtrace');
        warning('An error will be thrown at this point due to missing code.');
        error('fcn_geometry_intersectionGeom is not yet ready for curves from fit type: %s',firstFitType);
end
end % Ends fcn_INTERNAL_intersectLineGeomtries

%% fcn_INTERNAL_intersectLineCircle
function intersection_points = fcn_INTERNAL_intersectLineCircle(line_parameters, circle_parameters)

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
line_unit_tangent_vector     = line_parameters(1,1:2);
line_base_point_xy           = line_parameters(1,3:4);

% Get the line circle details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
circle_center_xy                = circle_parameters(1,1:2);
circle_radius                   = circle_parameters(1,3);

if line_unit_tangent_vector(1)==0
    slope = inf;
    intercept = line_base_point_xy(1);
else
    slope = line_unit_tangent_vector(2)/line_unit_tangent_vector(1);
    intercept = line_base_point_xy(2) - slope*line_base_point_xy(1); % b = y - m*x
end
% Use MATLAB's linecirc algorithm to find intersections
[xout,yout] = linecirc(slope,intercept,circle_center_xy(1,1),circle_center_xy(1,2),circle_radius);

if ~all(isnan(xout))
    % intersection points were found!
    intersection_points = [xout', yout'];

    % Keep only the points that are not nan
    intersection_points = intersection_points(~isnan(intersection_points(:,1)),:);

else
    intersection_points = [nan nan];
end

% Remove repeats
if isequal(length(intersection_points(:,1)),2)
    if isequal(intersection_points(1,:),intersection_points(2,:))
        intersection_points = unique(intersection_points,'rows');
    end
end
end % Ends fcn_INTERNAL_intersectLineCircle


%% fcn_INTERNAL_intersectLineArc
function intersection_points = fcn_INTERNAL_intersectLineArc(line_parameters, arc_parameters)

intersection_points = fcn_INTERNAL_intersectLineCircle(line_parameters, arc_parameters); 

if ~any(isnan(intersection_points))
    % Keep only the closet point on the arc. Returns nan if there are none
    [~, intersection_points] = fcn_INTERNAL_findPointsInArc(arc_parameters,intersection_points);
end

end % Ends fcn_INTERNAL_intersectLineArc

%% fcn_INTERNAL_intersectSegmentLine
function intersection_points = fcn_INTERNAL_intersectLineLine(first_segment_parameters, second_segment_parameters)

% Get the first segment fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_segment_unit_tangent_vector   = first_segment_parameters(1,1:2);
first_segment_base_point_xy         = first_segment_parameters(1,3:4);
first_segment_s_start               = first_segment_parameters(1,5);
first_segment_s_end                 = first_segment_parameters(1,6);
first_segment_start_xy              = first_segment_base_point_xy + first_segment_unit_tangent_vector*first_segment_s_start;
first_segment_end_xy                = first_segment_base_point_xy + first_segment_unit_tangent_vector*first_segment_s_end;

% Get the second segment fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_segment_unit_tangent_vector   = second_segment_parameters(1,1:2);
second_segment_base_point_xy         = second_segment_parameters(1,3:4);
second_segment_s_start               = second_segment_parameters(1,5);
second_segment_s_end                 = second_segment_parameters(1,6);
second_segment_start_xy              = second_segment_base_point_xy + second_segment_unit_tangent_vector*second_segment_s_start;
second_segment_end_xy                = second_segment_base_point_xy + second_segment_unit_tangent_vector*second_segment_s_end;

% For segment to segment interesection
flag_search_type = 4;

[~,intersection_points] = fcn_Path_findProjectionHitOntoPath([first_segment_start_xy; first_segment_end_xy],second_segment_start_xy,second_segment_end_xy,flag_search_type,(-1));

end % Ends fcn_INTERNAL_intersectLineLine

%% fcn_INTERNAL_intersectSegmentLine
function intersection_points = fcn_INTERNAL_intersectLineSegment(first_segment_parameters, second_segment_parameters)

% Get the first segment fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_segment_unit_tangent_vector   = first_segment_parameters(1,1:2);
first_segment_base_point_xy         = first_segment_parameters(1,3:4);
first_segment_s_start               = first_segment_parameters(1,5);
first_segment_s_end                 = first_segment_parameters(1,6);
first_segment_start_xy              = first_segment_base_point_xy + first_segment_unit_tangent_vector*first_segment_s_start;
first_segment_end_xy                = first_segment_base_point_xy + first_segment_unit_tangent_vector*first_segment_s_end;

% Get the second segment fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_segment_unit_tangent_vector   = second_segment_parameters(1,1:2);
second_segment_base_point_xy         = second_segment_parameters(1,3:4);
second_segment_s_start               = second_segment_parameters(1,5);
second_segment_s_end                 = second_segment_parameters(1,6);
second_segment_start_xy              = second_segment_base_point_xy + second_segment_unit_tangent_vector*second_segment_s_start;
second_segment_end_xy                = second_segment_base_point_xy + second_segment_unit_tangent_vector*second_segment_s_end;

% For segment to segment interesection
flag_search_type = 3;

[~,intersection_points] = fcn_Path_findProjectionHitOntoPath([first_segment_start_xy; first_segment_end_xy],second_segment_start_xy,second_segment_end_xy,flag_search_type,(-1));

end % Ends fcn_INTERNAL_intersectLineSegment

%% fcn_INTERNAL_intersectSegmentGeomtries
function  intersection_points = fcn_INTERNAL_intersectSegmentGeomtries(secondFitType, first_segment_parameters, secondFitType_parameters)
switch lower(secondFitType)
    case 'circle'
        intersection_points = fcn_INTERNAL_intersectSegmentCircle(first_segment_parameters, secondFitType_parameters);
    case 'arc'
        intersection_points = fcn_INTERNAL_intersectSegmentArc(first_segment_parameters, secondFitType_parameters);
    case 'line'
        % error('need to write this')
        intersection_points = fcn_INTERNAL_intersectSegmentLine(first_segment_parameters, secondFitType_parameters);
    case 'segment'       
        intersection_points = fcn_INTERNAL_intersectSegmentSegment(first_segment_parameters,secondFitType_parameters);
    case 'spiral'
        warning('on','backtrace');
        warning('An error will be thrown at this point due to missing code.');
        error('fcn_geometry_intersectionGeom case is not yet ready for spiral case');
    otherwise
        warning('on','backtrace');
        warning('An error will be thrown at this point due to missing code.');
        error('fcn_geometry_intersectionGeom is not yet ready for curves from fit type: %s',firstFitType);
end
end % Ends fcn_INTERNAL_intersectSegmentGeomtries

%% fcn_INTERNAL_intersectSegmentCircle
function intersection_points = fcn_INTERNAL_intersectSegmentCircle(segment_parameters, circle_parameters)
intersection_points = fcn_INTERNAL_intersectLineCircle(segment_parameters, circle_parameters); 

if ~any(isnan(intersection_points))
    % Keep only the closet point on the segment. Returns nan if there are none
    [~, intersection_points] = fcn_INTERNAL_findPointsInSegment(segment_parameters,intersection_points);
end
end % Ends fcn_INTERNAL_intersectSegmentCircle

%% fcn_INTERNAL_findPointsInSegment
function [points_to_test, closest_point_to_start_of_segment] = fcn_INTERNAL_findPointsInSegment(segment_parameters, points_to_test)
% Get the segment fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_unit_tangent_vector   = segment_parameters(1,1:2);
segment_base_point_xy         = segment_parameters(1,3:4);
segment_s_start               = segment_parameters(1,5);
segment_s_end                 = segment_parameters(1,6);
% segment_start_xy              = segment_base_point_xy + segment_unit_tangent_vector*segment_s_start;
% segment_end_xy                = segment_base_point_xy + segment_unit_tangent_vector*segment_s_end;

% Are the intersections within the segment's range?
N_intersections = length(points_to_test(:,1));

% Find distances to the segment's start point
vectors_from_segment_base_point = points_to_test - ones(N_intersections,1)*segment_base_point_xy;

% Distances are dot product with the line's vector
s_distances = sum(ones(N_intersections,1)*segment_unit_tangent_vector.*vectors_from_segment_base_point,2);
    
good_distance_indicies = intersect(find(s_distances>=segment_s_start),find(s_distances<=segment_s_end));
if isempty(good_distance_indicies)
    points_to_test = [nan nan];
    closest_point_to_start_of_segment = [nan nan];
else
    points_to_test = points_to_test(good_distance_indicies,:);

    good_s_distances = s_distances(good_distance_indicies,1);
    closest_to_start = min(good_s_distances);
    closest_point_to_start_of_segment = segment_base_point_xy + segment_unit_tangent_vector*closest_to_start;
end

end % Ends fcn_INTERNAL_findPointsInSegment


%% fcn_INTERNAL_intersectSegmentArc
function intersection_points = fcn_INTERNAL_intersectSegmentArc(segment_parameters, arc_parameters)
intersection_points = fcn_INTERNAL_intersectLineCircle(segment_parameters, arc_parameters); 

if ~all(isnan(intersection_points))

    % Keep only the points on the arc
    intersection_points = fcn_INTERNAL_findPointsInArc(arc_parameters,intersection_points);

    if ~all(isnan(intersection_points))

        % Keep only the closet point on the segment. Returns nan if there are none
        [~, intersection_points] = fcn_INTERNAL_findPointsInSegment(segment_parameters,intersection_points);
    end
end
end % Ends fcn_INTERNAL_intersectSegmentArc

%% fcn_INTERNAL_intersectSegmentLine
function intersection_points = fcn_INTERNAL_intersectSegmentLine(first_segment_parameters, second_segment_parameters)

% Get the first segment fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_segment_unit_tangent_vector   = first_segment_parameters(1,1:2);
first_segment_base_point_xy         = first_segment_parameters(1,3:4);
first_segment_s_start               = first_segment_parameters(1,5);
first_segment_s_end                 = first_segment_parameters(1,6);
first_segment_start_xy              = first_segment_base_point_xy + first_segment_unit_tangent_vector*first_segment_s_start;
first_segment_end_xy                = first_segment_base_point_xy + first_segment_unit_tangent_vector*first_segment_s_end;

% Get the second segment fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_segment_unit_tangent_vector   = second_segment_parameters(1,1:2);
second_segment_base_point_xy         = second_segment_parameters(1,3:4);
second_segment_s_start               = second_segment_parameters(1,5);
second_segment_s_end                 = second_segment_parameters(1,6);
second_segment_start_xy              = second_segment_base_point_xy + second_segment_unit_tangent_vector*second_segment_s_start;
second_segment_end_xy                = second_segment_base_point_xy + second_segment_unit_tangent_vector*second_segment_s_end;

% For segment to segment interesection
flag_search_type = 1;

[~,intersection_points] = fcn_Path_findProjectionHitOntoPath([first_segment_start_xy; first_segment_end_xy],second_segment_start_xy,second_segment_end_xy,flag_search_type,(-1));

end % Ends fcn_INTERNAL_intersectSegmentLine

%% fcn_INTERNAL_intersectSegmentSegment
function intersection_points = fcn_INTERNAL_intersectSegmentSegment(first_segment_parameters, second_segment_parameters)

% Get the first segment fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
first_segment_unit_tangent_vector   = first_segment_parameters(1,1:2);
first_segment_base_point_xy         = first_segment_parameters(1,3:4);
first_segment_s_start               = first_segment_parameters(1,5);
first_segment_s_end                 = first_segment_parameters(1,6);
first_segment_start_xy              = first_segment_base_point_xy + first_segment_unit_tangent_vector*first_segment_s_start;
first_segment_end_xy                = first_segment_base_point_xy + first_segment_unit_tangent_vector*first_segment_s_end;

% Get the second segment fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
second_segment_unit_tangent_vector   = second_segment_parameters(1,1:2);
second_segment_base_point_xy         = second_segment_parameters(1,3:4);
second_segment_s_start               = second_segment_parameters(1,5);
second_segment_s_end                 = second_segment_parameters(1,6);
second_segment_start_xy              = second_segment_base_point_xy + second_segment_unit_tangent_vector*second_segment_s_start;
second_segment_end_xy                = second_segment_base_point_xy + second_segment_unit_tangent_vector*second_segment_s_end;

% For segment to segment interesection
flag_search_type = 0;

[~,intersection_points] = fcn_Path_findProjectionHitOntoPath([first_segment_start_xy; first_segment_end_xy],second_segment_start_xy,second_segment_end_xy,flag_search_type,(-1));

end % Ends fcn_INTERNAL_intersectSegmentSegment

