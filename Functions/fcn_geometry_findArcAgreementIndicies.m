function [indicies_in_station_agreement, flag_is_a_circle, start_angle_in_radians, end_angle_in_radians] = ...
    fcn_geometry_findArcAgreementIndicies(points, circleCenter, circleRadius, index_source_point, station_tolerance, varargin)
% fcn_geometry_findArcAgreementIndicies
%
% Given a set of points, a circle center and radius, and the index of a
% source point within the set of points, finds the indicies of the points
% that are within a station_tolerance distance from each other, as measured
% in arc length along the circle, around the source point. This is useful
% to determine sets of points that fall on or around a circle perimeter,
% centered near a test point, that form a "contiguous" arc but perhaps not
% completely around the circle. Contiguous is defined as having an arc
% distance, when the points are projected onto the circle, less than or
% equal to the station tolerance.
%
% FORMAT:
% 
% [indicies_in_station_agreement, flag_is_a_circle, start_angle_in_radians, end_angle_in_radians] = ...
%    fcn_geometry_findArcAgreementIndicies(points, circleCenter, circleRadius, index_source_point, station_tolerance, (fig_num));
% 
% INPUTS:
%      points: a Nx2 vector where N is the number of points, but at least 2 rows. 
%
%      circleCenter: an 1x2 vector containing the [X Y] location of the
%      circle used to define the arc
%
%      circleRadius: an 1x1 scalar defining the radius of the circle used
%      to define the arc
%
%      index_source_point: the index of the point that will be the "root"
%      point for searching for a contiguous arc among all the points. This
%      is needed for the common situations where point data may include
%      many different arc portions along the same circle, and thus one must
%      define which arc is desired.
%
%      station_tolerance: the projection distance between the points in a
%      curve fit, along the direction of the line, that indicate whether a
%      point "belongs" to the arc (if distance is less than or equal
%      to the tolerance), or is "outside" the fit (if distance is greater
%      than the tolerance). 
% 
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      indicies_in_station_agreement: the indicies on the points that are
%      in a contiguous arc that contains the source_point
%
%      flag_is_a_circle: a flag that is 1 if the arc extends completely
%      around the perimeter thus forming a complete circle. It is 0
%      otherwise.
% 
%      start_angle_in_radians: the start angle of the arc
%
%      end_angle_in_radians: the end angle of the arc.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_circleCenterFrom3Points
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_findArcAgreementIndicies
% for a full test suite.
%
% This function was written on 2024_01_09 - S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_01_09 - S. Brennan
% -- first write of the code
% 2024_01_12 - S. Brennan
% -- fixed typo in test script name


%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==6 && isequal(varargin{end},-1))
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
    debug_fig_num = 234343; %#ok<NASGU>
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
        narginchk(5,6);

        % Check the points input to be length greater than or equal to 2
        fcn_DebugTools_checkInputsToFunctions(...
            points, '2column_of_numbers',[2 3]);

        % Check the circleCenter input is [2 x 1]
        fcn_DebugTools_checkInputsToFunctions(circleCenter, '2column_of_numbers',1);

        % Check the circleRadius input is a positive single number
        fcn_DebugTools_checkInputsToFunctions(circleRadius, 'positive_1column_of_numbers',1);

        % Check the index_source_point input is a positive integer
        fcn_DebugTools_checkInputsToFunctions(index_source_point, 'positive_1column_of_integers',1);

        % Check the station_tolerance input is a positive single number
        fcn_DebugTools_checkInputsToFunctions(station_tolerance, 'positive_1column_of_numbers',1);

    end
end


% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (6<=nargin) && (0==flag_max_speed)
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

% Method used here is to get the angles of all the points, and the
% angle of the "source" point, e.g. the point that starts the fit
% (the first point of the fitted source points). These angles are
% padded on front/back, as the arc may connect all the around to a
% circle either before the "source" point or after. The angles are
% then used to find a change in angle (e.g. arc length), which is
% used to find distance along the circle subtended by each change
% in angle. These distances are then converted into a distance
% vector which represents the station points around the
% circumference. These station points are then checked for
% continuity (similar to the line-fitting code) to find
% contiguous arcs.

% % Find arc of test segment
% [~, arc_angle_in_radians_1_to_3, ~, ~, start_angle_in_radians] = ...
%     fcn_geometry_arcAngleFrom3Points(test_source_points(1,:), test_source_points(2,:), test_source_points(3,:));
% start_angle_in_radians = mod(start_angle_in_radians,2*pi);
% end_angle_in_radians = start_angle_in_radians + arc_angle_in_radians_1_to_3;

% Find angles of the points, relative to center
shifted_points = (points - circleCenter);
point_angles = atan2(shifted_points(:,2),shifted_points(:,1));

% Make sure all the angles are positive
point_angles = mod(point_angles,2*pi);

% Sort the angles in ascending order
[sorted_angles, sorted_index] = sort(point_angles);

% Find where the source point index moved after sorting. This is
% used to determine where points need to be padded
sorted_source_index = find(sorted_index == index_source_point);

% Pad the indicies and the point angles to allow wrap-around. This is so
% that the arc, if it continues before or after the index_source_point,
% will still be found. In effect, it makes the index_source_point be the
% first angle, one of the angles in the middle, and the last angle. This
% way, any arc will be found that contains this index_source_point.
padded_sorted_indicies = sorted_index;
padded_sorted_angles = sorted_angles;

% Pad the front end
if sorted_source_index~=length(sorted_index)
    padded_sorted_angles = [sorted_angles(sorted_source_index+1:end)-2*pi; padded_sorted_angles];
    padded_sorted_indicies = [sorted_index(sorted_source_index+1:end); padded_sorted_indicies];
end

% Pad the back end
if sorted_source_index~=1
    padded_sorted_angles = [padded_sorted_angles; sorted_angles(1:sorted_source_index)+2*pi];
    padded_sorted_indicies = [padded_sorted_indicies; sorted_index(1:sorted_source_index)];
end
% padded_sorted_angles   = sorted_angles(padded_sorted_indicies);
padded_sorted_source_index = find(padded_sorted_indicies(1:end-1) == index_source_point);

if isempty(padded_sorted_source_index) || length(padded_sorted_source_index)>1 
    error('Unexpected index found.');
end

% Find the angle differences
diff_angles = diff(padded_sorted_angles);

if any(diff_angles<0)
    error('negative angles detected - this should not happen!');
end

% Find station distance increments
diff_distances = diff_angles*circleRadius;

% Find stations
station_distances_of_points_in_radial_agreement = [0; cumsum(diff_distances)];

% Sort the station distances and find those in agreement with station
% tolerance
indicies_in_station_agreement_sorted = ...
    fcn_geometry_findPointsInSequence(...
    station_distances_of_points_in_radial_agreement, ...
    padded_sorted_source_index, ...
    station_tolerance, -1);

indicies_in_input_form = padded_sorted_indicies(indicies_in_station_agreement_sorted);
start_angle_in_radians = padded_sorted_angles(indicies_in_station_agreement_sorted(1));
end_angle_in_radians   = padded_sorted_angles(indicies_in_station_agreement_sorted(end));

indicies_in_station_agreement = unique(indicies_in_input_form,'stable');

% Check if any repeated - if so, it is a circle!
flag_is_a_circle = 0;
if length(indicies_in_input_form) ~= length(indicies_in_station_agreement)
    flag_is_a_circle = 1;
    indicies_in_station_agreement = sorted_index;
    start_angle_in_radians = sorted_angles(1);
    end_angle_in_radians   = sorted_angles(end);

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

    % Plot the results in point space
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end        

    hold on;
    grid on;
    title('Testing for contiguous arc');
    xlabel('X [meters]');
    ylabel('Y [meters]')

    % Plot the input points
    plot(points(:,1),points(:,2),'k.','MarkerSize',20);

    % Plot the circle fit
    if flag_is_a_circle==1
        title('Circle fit');
    else
        title('Arc fit');
    end

    fcn_geometry_plotCircle(circleCenter, circleRadius, 'b-',fig_num)
    plot(points(index_source_point,1),points(index_source_point,2),'bo','MarkerSize',15);

    % Plot the results
    % [indicies_in_station_agreement, flag_is_a_circle, start_angle_in_radians, end_angle_in_radians]
    plot(points(indicies_in_station_agreement,1),points(indicies_in_station_agreement,2),'r.','MarkerSize',15);

    % Label the points
    text(points(indicies_in_station_agreement(1),1),    points(indicies_in_station_agreement(1),2),    sprintf('Start angle: %.3f deg', start_angle_in_radians*180/pi));
    text(points(indicies_in_station_agreement(end),1),  points(indicies_in_station_agreement(end),2),  sprintf('End angle: %.3f deg', end_angle_in_radians*180/pi));

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
