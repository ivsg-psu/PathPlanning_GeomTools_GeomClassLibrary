function [agreement_indicies, station_distances] = fcn_geometry_findAgreementsOfPointsToLineVector( points, unit_projection_vector, base_point_index, transverse_tolerance, station_tolerance, varargin)
% fcn_geometry_findAgreementsOfPointsToLineVector
% Given a set of XY points, a unit projection vector and a base point that the vector is
% attached to, find the indicies of the points that are within a
% transverse_tolerance distance away from the vector. Then, among these
% points, finds the indicies that are within a station_tolerance distance
% from each other in a group that includes the base point.
%
% Format:
% agreement_indicies = fcn_geometry_findAgreementsOfPointsToLineVector( points, unit_projection_vector, base_point_index, transverse_tolerance, station_tolerance, (fig_num))
%
% INPUTS:
%      points: a Nx2 vector where N is the number of points, but at least 2
%      rows.
%
%      unit_projection_vector: a 1x2 vector that is a unit vector in the
%      direction to check.
%
%      base_point_index: the index of the point to use as a base point,
%      among the input points.
%
%      transverse_tolerance: the orthogonal distance between the points and
%      the linear vector fit that indicate whether a point "belongs" to the
%      fit. A point belongs to the fit if the transverse distance is less
%      than or equal to the transverse_tolerance
%
%      station_tolerance: the projection distance between the points in a
%      vector fit, along the direction of the line, that indicate whether a
%      point "belongs" to the linear fit of a line segment or not. A line
%      segment is considered continous if station interval distance is less
%      than or equal to the station_tolerance. Set station_tolerance to an
%      empty value, [], to avoid line segment fitting and instead find pure
%      line fits along the vector.
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      agreement_indicies: the indicies of the points that are within
%      agreement of the best-fit parameters, given the transverse and
%      station tolerance settings. If a station_tolerance is specified, the
%      indicies are sorted in station-increasing order.
% 
%      station_distances: a [1 x 2] vector of the minimum and maximum
%      station distances, relative to the base point, that are in agreement
%
% DEPENDENCIES:
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_findPointsInSequence
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_findAgreementsOfPointsToLineVector
% for a full test suite.

% This function was written on 2024_01_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2023_01_14 - S. Brennan
% -- wrote the code
% 2023_01_15 - S. Brennan
% -- added station distances output

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
        narginchk(5,6);

        % Check the points input to be length greater than or equal to 2
        fcn_DebugTools_checkInputsToFunctions(...
            points, '2column_of_numbers',[2 3]);

        % Check the unit_projection_vector input to be [1x2]
        fcn_DebugTools_checkInputsToFunctions(...
            unit_projection_vector, '2column_of_numbers',[1 1]);

        % Check the base_point_index input to be [1x1] positive integer or
        % a [1x2] point
        if length(base_point_index)==1
            fcn_DebugTools_checkInputsToFunctions(...
                base_point_index, 'positive_1column_of_integers',[1 1]);
            if base_point_index<1 || base_point_index>length(points(:,1))
                error('The input base_point_index must be in the length of points, namely between: %.0d and %.0d',1,length(points(:,1)));
            end
        else
            % Check that it is a 1x2
            fcn_DebugTools_checkInputsToFunctions(...
                base_point_index, '2column_of_numbers',[1 1]);
        end

        % Check the transverse_tolerance input is a positive single number
        fcn_DebugTools_checkInputsToFunctions(transverse_tolerance, 'positive_1column_of_numbers',1);

        % Check the station_tolerance input is a positive single number
        if ~isempty(station_tolerance)
            fcn_DebugTools_checkInputsToFunctions(station_tolerance, 'positive_1column_of_numbers',1);
        end
    end
end

% Does user want to specify fig_num?
flag_do_plots = 0;
if 6<= nargin && 0==flag_max_speed
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

% Save the original point length
N_points = length(points(:,1));

% Do we need to ammend a point to the points? This occurs when the user
% enters a base point that is not in the points list originally, via the
% base_point_index

% Check the base_point_index input to be [1x1] positive integer or
% a [1x2] point
flag_was_amended = 0;
if length(base_point_index)==1
    base_point = points(base_point_index,:);
else
    % Amend the points array
    flag_was_amended = 1;
    points = [points; base_point_index];
    base_point = base_point_index;
    base_point_index = N_points+1;
end


% Find the indicies in transverse agreement. This is done by finding
% the lateral distances of the points to the fit line, and then keeping
% only the points within the lateral distance.
unit_orthogonal_vector = unit_projection_vector*[0 1; -1 0];
base_projection_vectors = points - base_point;
lateral_distances = sum(unit_orthogonal_vector.*base_projection_vectors,2);

indicies_in_lateral_agreement = find((abs(lateral_distances)<=transverse_tolerance)');

% Next find the indicies in station agreement. To do this, find the
% station coordinates for every point as predicted by their tangent
% distances. Then sort the points by this distance, keeping track of
% how indicies move when sorted. Then check differences in station
% change from point to point, cutting off the fit if and when these
% differences are larger than the station tolerance.

if ~isempty(station_tolerance)
    % Find station distances as projections of tangent vectors
    station_distances_of_points_in_lateral_agreement = ...
        sum(base_projection_vectors(indicies_in_lateral_agreement,:).* unit_projection_vector,2);

    % Find where the base point is in the rearranged list
    base_point_index_in_lateral_agreement = find(indicies_in_lateral_agreement == base_point_index);

    % Sort the station distances and find those in agreement with station
    % tolerance
    indicies_in_station_agreement = ...
        fcn_geometry_findPointsInSequence(...
        station_distances_of_points_in_lateral_agreement, ...
        base_point_index_in_lateral_agreement, ...
        station_tolerance, -1);

    % Find indicies in both lateral and station agreement
    if flag_was_amended
        indicies_in_station_agreement = indicies_in_station_agreement(indicies_in_station_agreement~=base_point_index_in_lateral_agreement);
    end

    indicies_in_both_lateral_and_station_agreement = indicies_in_lateral_agreement(indicies_in_station_agreement);
 
    % station_distances_of_points_in_lateral_and_station_agreement = ...
    %     sum(base_projection_vectors(indicies_in_both_lateral_and_station_agreement,:).* unit_projection_vector,2);
    % 
    % station_distances = [min(station_distances_of_points_in_lateral_and_station_agreement) max(station_distances_of_points_in_lateral_and_station_agreement)];
else
    indicies_in_both_lateral_and_station_agreement = indicies_in_lateral_agreement;
    % station_distances = [-inf inf];
end

station_distances_of_points_in_lateral_and_station_agreement = ...
    sum(base_projection_vectors(indicies_in_both_lateral_and_station_agreement,:).* unit_projection_vector,2);

station_distances = [min(station_distances_of_points_in_lateral_and_station_agreement) max(station_distances_of_points_in_lateral_and_station_agreement)];
agreement_indicies = indicies_in_both_lateral_and_station_agreement;



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

    title('Points and agreements, plotted in point-space');
    xlabel('X [meters]');
    ylabel('Y [meters]')


    % Plot the input points
    plot(points(:,1),points(:,2),'k.','MarkerSize',30)

    % Plot the unit vector
    quiver(base_point(:,1),base_point(:,2),unit_projection_vector(:,1),unit_projection_vector(:,2),0,'g','Linewidth',6);

    % Plot the fit
    if ~any(isinf(station_distances))
        fit_points = [unit_projection_vector*station_distances(1); unit_projection_vector*station_distances(2)] + ones(2,1)*base_point;
        plot(fit_points(:,1),fit_points(:,2),'-','Linewidth',3,'Color',[0 0.5 0]);
    else
        distances_base_point = sum((points - base_point).^2,2).^0.5;
        longest_distance = max(distances_base_point);
        fit_points = [-unit_projection_vector*longest_distance; unit_projection_vector*longest_distance]+ ones(2,1)*base_point;
        plot(fit_points(:,1),fit_points(:,2),'-','Linewidth',3,'Color',[0 0.5 0]);

    end

    % Plot the boundaries
    base_projection_vectors = points - base_point;

    station_distances_of_points_in_lateral_agreement = ...
        sum(base_projection_vectors(indicies_in_lateral_agreement,:).* unit_projection_vector,2);
    max_station = max(station_distances_of_points_in_lateral_agreement);
    min_station = min(station_distances_of_points_in_lateral_agreement);
    boundary_start = base_point + unit_projection_vector*min_station;
    boundary_end   = base_point + unit_projection_vector*max_station;

    boundary_left  = [boundary_start; boundary_end] + ones(2,1)*unit_orthogonal_vector*transverse_tolerance;
    boundary_right = [boundary_start; boundary_end] - ones(2,1)*unit_orthogonal_vector*transverse_tolerance;
    plot(boundary_left(:,1),boundary_left(:,2),'-','Color',[0.7 0.7 0.7],'LineWidth',3);
    plot(boundary_right(:,1),boundary_right(:,2),'-','Color',[0.7 0.7 0.7],'LineWidth',3);


    % Plot the results
    plot(points(indicies_in_lateral_agreement,1),points(indicies_in_lateral_agreement,2),'m.','MarkerSize',30)
    plot(points(indicies_in_both_lateral_and_station_agreement,1),points(indicies_in_both_lateral_and_station_agreement,2),'c.-','MarkerSize',20,'LineWidth',3)

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



