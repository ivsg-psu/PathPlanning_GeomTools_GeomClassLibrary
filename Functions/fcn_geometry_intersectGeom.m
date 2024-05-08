function intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, varargin) 
%% fcn_geometry_intersectGeom
% 
% This function finds the intersection point(s) of two geometries, such as
% line-arc, line-circle, line segment-arc, arc-arc, arc-circle, arc-line, etc.
% 
% For segments, arcs, and spirals that have multiple intersection points,
% return the point that is nearest to the start of the segment, arc, or
% spiral, with "nearest" meaning in station distance, not geometric
% position. In other words, returns the intersection that is encountered
% first when traversing the "from" geometry starting at its "start"
% position.
%
% For circles and lines that are the "starting" geometry, the intersection
% point is returned that is the first point encountered in the 2nd geometry
% if the 2nd geometry is a segment, arc, or spiral.
%
% For circles intersecting with circles or lines, both points are returned.
%
% FORMAT: 
%
% intersection_points = fcn_geometry_intersectGeom(firstFitType,  firstFitType_parameters, secondFitType,  secondFitType_parameters, fig_num) 
% 
% INPUTS:
%
% firstFitType: The string input of the first geometry
%
% firstFitType_parameters: The parameters of the first geometry
%
% secondFitType: The string input of the second geometry
%
% secondFitType_parameters: The parameters of the second geometry
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
% This function was written on 2024_05_02 by Aneesh Batchu
% Questions or comments? abb6486@psu.edu or sbrennan@psu.edu

% Revision History
% 2024_05_02 - Aneesh Batchu
% -- wrote the code 
% 2024_05_06 - Aneesh Batchu 
% -- Fixed BUGS in Arc-Arc, line-arc, arc-line, line segment-arc
% intersection cases. Added a few conditional statements to remove NaNs
% from the potential_intersection_points to compute the true "intersection
% points"


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

    case 'arc'

        % Calculate needed values from parameter sets
        % Get the arc fit details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
        arc1_center_xy                = firstFitType_parameters(1,1:2);
        arc1_radius                   = firstFitType_parameters(1,3);
        arc1_start_angle_in_radians   = firstFitType_parameters(1,4);
        arc1_end_angle_in_radians     = firstFitType_parameters(1,5);
        arc1_start_xy                 = arc1_center_xy + arc1_radius*[cos(arc1_start_angle_in_radians) sin(arc1_start_angle_in_radians)];
        arc1_end_xy                   = arc1_center_xy + arc1_radius*[cos(arc1_end_angle_in_radians) sin(arc1_end_angle_in_radians)];
        arc1_is_counter_clockwise     = firstFitType_parameters(1,7);

        switch lower(secondFitType)

            case 'line'
                % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
                line_unit_tangent_vector     = secondFitType_parameters(1,1:2);
                line_base_point_xy           = secondFitType_parameters(1,3:4);
                % line_s_start               = clean_line_parameters(1,5);
                % line_s_end                 = clean_line_parameters(1,6);
                % line_start_xy              = line_base_point_xy + line_unit_tangent_vector*line_s_start;
                % line_end_xy                = line_base_point_xy + line_unit_tangent_vector*line_s_end;

                if line_unit_tangent_vector(1)==0
                    slope = inf;
                    intercept = line_base_point_xy(1);
                else
                    slope = line_unit_tangent_vector(2)/line_unit_tangent_vector(1);
                    intercept = line_base_point_xy(2) - slope*line_base_point_xy(1); % b = y - m*x
                end
                % Use MATLAB's linecirc algorithm to find intersections
                [xout,yout] = linecirc(slope,intercept,arc1_center_xy(1,1),arc1_center_xy(1,2),arc1_radius);

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
                        intersection_is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(arc1_start_xy, two_intersection_points(ith_row,:), arc1_end_xy,-1);
                        if arc1_is_counter_clockwise == intersection_is_counterClockwise
                            potential_intersection_points(ith_row,:) = two_intersection_points(ith_row,:);
                        else
                            potential_intersection_points(ith_row,:) = [nan nan];
                        end
                    end

                    if isequal(potential_intersection_points(1,:),potential_intersection_points(2,:))
                        intersection_points = potential_intersection_points(1,:);
                    else
                        % Find which point is closest to the line's start point
                        vectors_from_line_base_point = potential_intersection_points - [1;1]*line_base_point_xy;

                        % Distances are dot product with the line's vector
                        distances = sum([1; 1]*line_unit_tangent_vector.*vectors_from_line_base_point,2).^0.5;
                        if distances(1)<distances(2)
                            intersection_points = potential_intersection_points(1,:);
                        elseif distances(2)<distances(1)
                            intersection_points = potential_intersection_points(2,:);
                        else
                            nan_indices = any(isnan(potential_intersection_points), 2);
                            potential_intersection_points = potential_intersection_points(~nan_indices,:);
                            if ~isempty(potential_intersection_points)
                                intersection_points = potential_intersection_points;
                            else
                                intersection_points = [nan nan];
                            end
                        end
                    end

                else
                    intersection_points = [nan nan];
                end
                

            case 'arc'

                % Get the arc fit details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
                arc2_center_xy                = secondFitType_parameters(1,1:2);
                arc2_radius                   = secondFitType_parameters(1,3);
                arc2_start_angle_in_radians   = secondFitType_parameters(1,4);
                arc2_end_angle_in_radians     = secondFitType_parameters(1,5);
                arc2_start_xy                 = arc2_center_xy + arc2_radius*[cos(arc2_start_angle_in_radians) sin(arc2_start_angle_in_radians)];
                arc2_end_xy                   = arc2_center_xy + arc2_radius*[cos(arc2_end_angle_in_radians) sin(arc2_end_angle_in_radians)];
                arc2_is_counter_clockwise     = secondFitType_parameters(1,7);

                % Use MATLAB's circcirc algorithm to find intersections between two circles
                [xout,yout] = circcirc(arc1_center_xy(1,1),arc1_center_xy(1,2),arc1_radius,arc2_center_xy(1,1),arc2_center_xy(1,2),arc2_radius);

                % % Check results of above
                % if 1==1
                %     figure(233);
                %     clf;
                %     hold on;
                %     grid on;
                %     axis equal
                % 
                %     % Plot the circles
                %     fcn_geometry_plotCircle(arc1_center_xy, arc1_radius,'g-',(233));
                %     fcn_geometry_plotCircle(arc2_center_xy, arc2_radius,'r-',(233));
                % 
                %     intersections = [xout', yout'];
                %     plot(intersections(:,1),intersections(:,2),'k.','MarkerSize',20);
                % end

                if ~isnan(xout)
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
                        intersection_points = potential_arc_intersection_points(1,:);
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
                            intersection_points = potential_arc_intersection_points(1,:);
                        elseif arc_angle_point2<arc_angle_point1
                            intersection_points = potential_arc_intersection_points(2,:);
                        else
                            nan_indices = any(isnan(potential_arc_intersection_points), 2);
                            potential_arc_intersection_points = potential_arc_intersection_points(~nan_indices,:);
                            if ~isempty(potential_arc_intersection_points)
                                 intersection_points = potential_arc_intersection_points;
                            else
                                intersection_points = [nan nan];
                            end
                        end
                    end

                else
                    intersection_points = [nan nan];
                end
                
            case 'spiral'
                error('fcn_geometry_intersectionGeom case is not yet ready for spiral case');
            otherwise
                warning('on','backtrace');
                warning('An error will be thrown at this point due to missing code.');
                error('fcn_geometry_intersectionGeom is not yet ready for curves from fit type: %s',firstFitType);

        end

    case 'line'

        % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
        line_unit_tangent_vector     = firstFitType_parameters(1,1:2);
        line_base_point_xy           = firstFitType_parameters(1,3:4);
        % line_s_start               = clean_line_parameters(1,5);
        % line_s_end                 = clean_line_parameters(1,6);
        % line_start_xy              = line_base_point_xy + line_unit_tangent_vector*line_s_start;
        % line_end_xy                = line_base_point_xy + line_unit_tangent_vector*line_s_end;
        switch lower(secondFitType)
            case 'circle'
                circle_center_xy                = secondFitType_parameters(1,1:2);
                circle_radius                   = secondFitType_parameters(1,3);

                if line_unit_tangent_vector(1)==0
                    slope = inf;
                    intercept = line_base_point_xy(1);
                else
                    slope = line_unit_tangent_vector(2)/line_unit_tangent_vector(1);
                    intercept = line_base_point_xy(2) - slope*line_base_point_xy(1); % b = y - m*x
                end
                % Use MATLAB's linecirc algorithm to find intersections
                [xout,yout] = linecirc(slope,intercept,circle_center_xy(1,1),circle_center_xy(1,2),circle_radius);

                if ~isnan(xout)
                    % intersection points were found!

                    % Which point(s) to keep?
                    intersection_points = [xout', yout'];
                    if isnan(intersection_points(1,:))
                        intersection_points = intersection_points(2,:);
                    elseif isnan(intersection_points(2,:))
                        intersection_points = intersection_points(1,:);
                    end

                else
                    intersection_points = [nan nan];
                end

            case 'arc'
                % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
                arc_center_xy                = secondFitType_parameters(1,1:2);
                arc_radius                   = secondFitType_parameters(1,3);
                arc_start_angle_in_radians   = secondFitType_parameters(1,4);
                arc_end_angle_in_radians     = secondFitType_parameters(1,5);
                arc_start_xy                 = arc_center_xy + arc_radius*[cos(arc_start_angle_in_radians) sin(arc_start_angle_in_radians)];
                arc_end_xy                   = arc_center_xy + arc_radius*[cos(arc_end_angle_in_radians) sin(arc_end_angle_in_radians)];
                % arc_is_circle                = clean_arc_parameters(1,6);
                arc_is_counter_clockwise     = secondFitType_parameters(1,7);

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
                        intersection_points = potential_intersection_points(1,:);
                    else
                        % Find which point is closest to the line's start point
                        vectors_from_line_base_point = potential_intersection_points - [1;1]*line_base_point_xy;

                        % Distances are dot product with the line's vector
                        distances = sum([1; 1]*line_unit_tangent_vector.*vectors_from_line_base_point,2).^0.5;
                        if distances(1)<distances(2)
                            intersection_points = potential_intersection_points(1,:);
                        elseif distances(2)<distances(1)
                            intersection_points = potential_intersection_points(2,:);
                        else
                            nan_indices = any(isnan(potential_intersection_points), 2);
                            potential_intersection_points = potential_intersection_points(~nan_indices,:);
                            if ~isempty(potential_intersection_points)
                                intersection_points = potential_intersection_points;
                            else
                                intersection_points = [nan nan];
                            end


                        end
                    end

                else
                    intersection_points = [nan nan];
                end

            case 'spiral'
                error('fcn_geometry_intersectionGeom case is not yet ready for spiral case');

            otherwise
                warning('on','backtrace');
                warning('An error will be thrown at this point due to missing code.');
                error('fcn_geometry_intersectionGeom is not yet ready for curves from fit type: %s',firstFitType);

        end

    case 'line segment'
        % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
        line_unit_tangent_vector     = firstFitType_parameters(1,1:2);
        line_base_point_xy           = firstFitType_parameters(1,3:4);
        line_s_start               = firstFitType_parameters(1,5);
        line_s_end                 = firstFitType_parameters(1,6);
        line_start_xy              = line_base_point_xy + line_unit_tangent_vector*line_s_start;
        line_end_xy                = line_base_point_xy + line_unit_tangent_vector*line_s_end;

        switch lower(secondFitType)
            case 'arc'
                % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
                arc_center_xy                = secondFitType_parameters(1,1:2);
                arc_radius                   = secondFitType_parameters(1,3);
                arc_start_angle_in_radians   = secondFitType_parameters(1,4);
                arc_end_angle_in_radians     = secondFitType_parameters(1,5);
                arc_start_xy                 = arc_center_xy + arc_radius*[cos(arc_start_angle_in_radians) sin(arc_start_angle_in_radians)];
                arc_end_xy                   = arc_center_xy + arc_radius*[cos(arc_end_angle_in_radians) sin(arc_end_angle_in_radians)];
                % arc_is_circle                = secondFitType_parameters(1,6);
                arc_is_counter_clockwise     = secondFitType_parameters(1,7);
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
                        intersection_points = potential_intersection_points(1,:);
                    else
                        % Find which point is closest to the line's start point
                        vectors_from_line_base_point = potential_intersection_points - [1;1]*line_base_point_xy;

                        % Distances are dot product with the line's vector
                        distances = sum([1; 1]*line_unit_tangent_vector.*vectors_from_line_base_point,2).^0.5;
                        % line_base_point_xy_matrix = [1; 1]*line_base_point_xy;
                        % distances = (sum((two_intersection_points - line_base_point_xy_matrix).^2, 2)).^0.5;
                        if distances(1)<distances(2)

                            distance = distances(1);
                            intersection_point_index = 1;
                            % intersection_points = potential_intersection_points(1,:);
                        elseif distances(2)<distances(1)
                            distance = distances(2);
                            intersection_point_index = 2;
                            % intersection_points = potential_intersection_points(2,:);
                        else
                            nan_indices = any(isnan(distances), 2);
                            distance = distances(~nan_indices,:);
                            if ~isempty(distance)
                                intersection_point_index = ~nan_indices;
                            else
                                distance = nan;
                            end
                        end

                        % Find the vector for start point to the end point
                        distance_btw_line_start_point_and_line_end_point = (sum((line_start_xy - line_end_xy).^2, 2)).^0.5;
                        
                        % Check if the distance of the intersection points
                        % are less than the length of the line segment
                        flag_compare_dist_of_intersection_points_to_length = distance <= distance_btw_line_start_point_and_line_end_point;

                        if ~isequal(flag_compare_dist_of_intersection_points_to_length, 0)
                            intersection_points = potential_intersection_points(intersection_point_index,:);
                        else
                            intersection_points = [nan nan];
                        end

                    end

                else
                    intersection_points = [nan nan];
                end

            case 'circle'
                circle_center_xy                = secondFitType_parameters(1,1:2);
                circle_radius                   = secondFitType_parameters(1,3);

                if line_unit_tangent_vector(1)==0
                    slope = inf;
                    intercept = line_base_point_xy(1);
                else
                    slope = line_unit_tangent_vector(2)/line_unit_tangent_vector(1);
                    intercept = line_base_point_xy(2) - slope*line_base_point_xy(1); % b = y - m*x
                end
                % Use MATLAB's linecirc algorithm to find intersections
                [xout,yout] = linecirc(slope,intercept,circle_center_xy(1,1),circle_center_xy(1,2),circle_radius);

                if ~isnan(xout)
                    % intersection points were found!

                    % Which point(s) to keep?
                    two_intersection_points = [xout', yout'];
    
                    % Find which point is closest to the line's start point
                    % vectors_from_line_base_point = two_intersection_points - [1;1]*line_base_point_xy;

                    % Distances are dot product with the line's vector
                    % distances = sum([1; 1]*line_unit_tangent_vector.*vectors_from_line_base_point,2).^0.5;
                    
                    line_base_point_xy_matrix = [1; 1]*line_base_point_xy;

                    distances = (sum((two_intersection_points - line_base_point_xy_matrix).^2, 2)).^0.5;


                    % Find the vector for start point to the end point
                    distance_btw_line_start_point_and_line_end_point = (sum((line_start_xy - line_end_xy).^2, 2)).^0.5;

                    % Check if the distance of the intersection points
                    % are less than the length of the line segment
                    flag_compare_dist_of_intersection_points_to_length = distances <= distance_btw_line_start_point_and_line_end_point;

                    if ~isequal(flag_compare_dist_of_intersection_points_to_length, [0;0])
                        intersection_points = two_intersection_points(flag_compare_dist_of_intersection_points_to_length,:);
                    else
                        intersection_points = [nan nan];
                    end

                    % if isnan(intersection_points(1,:))
                    %     intersection_points = intersection_points(2,:);
                    % elseif isnan(intersection_points(2,:))
                    %     intersection_points = intersection_points(1,:);
                    % end

                else
                    intersection_points = [nan nan];
                end



            case 'spiral'
                error('fcn_geometry_intersectionGeom case is not yet ready for spiral case');

            otherwise
                warning('on','backtrace');
                warning('An error will be thrown at this point due to missing code.');
                error('fcn_geometry_intersectionGeom is not yet ready for curves from fit type: %s',firstFitType);
        end
    case 'circle'

        circle_center_xy                = firstFitType_parameters(1,1:2);
        circle_radius                   = firstFitType_parameters(1,3);

        switch lower(secondFitType)
            case 'arc'
                % Get the arc fit details from arc2 parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
                arc2_center_xy                = secondFitType_parameters(1,1:2);
                arc2_radius                   = secondFitType_parameters(1,3);
                arc2_start_angle_in_radians   = secondFitType_parameters(1,4);
                arc2_end_angle_in_radians     = secondFitType_parameters(1,5);
                arc2_start_xy                 = arc2_center_xy + arc2_radius*[cos(arc2_start_angle_in_radians) sin(arc2_start_angle_in_radians)];
                arc2_end_xy                   = arc2_center_xy + arc2_radius*[cos(arc2_end_angle_in_radians) sin(arc2_end_angle_in_radians)];
                arc2_is_counter_clockwise     = secondFitType_parameters(1,7);

                % Use MATLAB's circcirc algorithm to find intersections between two circles
                [xout,yout] = circcirc(circle_center_xy(1,1),circle_center_xy(1,2),circle_radius,arc2_center_xy(1,1),arc2_center_xy(1,2),arc2_radius);

                if ~isnan(xout)
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

                        % % Check arc1
                        % intersection_is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(arc1_start_xy, circle_intersection_points(ith_row,:), arc1_end_xy,-1);
                        % if arc1_is_counter_clockwise ~= intersection_is_counterClockwise
                        %     potential_arc_intersection_points(ith_row,:) = [nan nan];
                        % end

                        % Check arc2
                        intersection_is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(arc2_start_xy, circle_intersection_points(ith_row,:), arc2_end_xy,-1);
                        if arc2_is_counter_clockwise ~= intersection_is_counterClockwise
                            potential_arc_intersection_points(ith_row,:) = [nan nan];
                        end

                    end

                    if isequal(potential_arc_intersection_points(1,:),potential_arc_intersection_points(2,:))
                        intersection_points = potential_arc_intersection_points(1,:);
                    else
                        
                        if isnan(potential_arc_intersection_points(2,:))
                            intersection_points = potential_arc_intersection_points(1,:);
                        elseif isnan(potential_arc_intersection_points(1,:))
                            intersection_points = potential_arc_intersection_points(2,:);
                        else
                            intersection_points = potential_arc_intersection_points;
                        end
                        % % Find which point is closest to arc1's start point
                        % if ~any(isnan(potential_arc_intersection_points(1,:)))
                        %     arc_angle_point1  = fcn_geometry_arcAngleFrom3Points(arc1_start_xy, potential_arc_intersection_points(1,:), arc1_end_xy,(-1));
                        % else
                        %     arc_angle_point1 = nan;
                        % end
                        % if ~any(isnan(potential_arc_intersection_points(2,:)))
                        %     arc_angle_point2  = fcn_geometry_arcAngleFrom3Points(arc1_start_xy, potential_arc_intersection_points(2,:), arc1_end_xy,(-1));
                        % else
                        %     arc_angle_point2 = nan;
                        % end
                        % 
                        % if arc_angle_point1<arc_angle_point2
                        %     intersection_points = potential_arc_intersection_points(1,:);
                        % else
                        %     intersection_points = potential_arc_intersection_points(2,:);
                        % end
                    end

                else
                    intersection_points = [nan nan];
                end
        end

    case 'spiral'
        error('fcn_geometry_intersectionGeom case is not yet ready for spiral case');
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
    fcn_geometry_plotGeometry(lower(firstFitType),firstFitType_parameters);
    fcn_geometry_plotGeometry(lower(secondFitType),secondFitType_parameters);
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

