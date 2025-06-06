function [flag_is_similar, minimum_distance_to_each_point, indicies_of_nearest_reference_points, mean_error, max_error, std_dev_error] = fcn_geometry_comparePointsToPoints(reference_points_XY, test_points_XY, varargin)
%% fcn_geometry_comparePointsToPoints
% Finds how "close" one set of points is to another.
%
% Given a reference set of points and a testing set of points, finds the
% minimum distance between each point in the testing set to the nearest
% point in the reference set. If an optional threshold is given, returns
% flag_is_similar = 1 (default) if all the distances are within a
% threshold. 
%
% Format:
% [flag_is_similar, ...
%  minimum_distance_to_each_point,...
%  indicies_of_nearest_reference_points, ...
%  mean_error, max_error, std_dev_error] = ...
% fcn_geometry_comparePointsToPoints(reference_points_XY, test_points_XY, ...
% (threshold), (fig_num))
%
% INPUTS:
%
%      reference_points_XY: an [Nx2] set of points that serve as the
%      reference, with N>=1.
%
%      test_points_XY: an [Mx2] set of points to compare to the reference.
%
%      (OPTIONAL INPUTS)
%
%      threshold: the maximum allowable distance, in meters, between
%      reference_points_XY and test_points_XY, where if any points are
%      larger than this distance, then flag_is_similar = 0. Otherwise, if
%      all test points are within the threshold distance, returns
%      flag_is_similar = 1;
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. 
%      NOTE: for plotting, the points are plotted so that point-pairs with
%      smallest distance are plotted first, and largest distance are
%      plotted last. The colors for each point-pair start mild and then get
%      more "red" as errors get larger, with maximum red values for any
%      points at or above the threshold. If no threshold is given, then the
%      maximum distance is made to be red, with all distances less than
%      this scaled accordingly on the color map.
%
% OUTPUTS:
%
%      flag_is_similar: a value of "true" if all points are within the
%      threshold (default if no threshold given), or "false" if any point
%      is larger than the threshold
%
%      minimum_distance_to_each_point: an [Mx1] set of the distances from
%      each test point to the closest reference point.
%
%      indicies_of_nearest_reference_points: an [Mx1] set of the indices of
%      the N reference points that were the closest point to each of the M
%      test points.
%
%      mean_error: the average error between the reference points and the
%      test points.
%
%      max_error: the maximum error between the reference points and the
%      test points
%
%      std_dev_error: the standard deviation in error between the reference
%      points and the test points
%
%
% DEPENDENCIES:
%      
%      (none)
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_comparePointsToPoints
% for a full test suite.
%
% This function was written on 2024_04_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_04_14 - S. Brennan
% -- wrote the code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{end},-1))
    flag_do_debug = 0; % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % Flag to plot the results for debugging
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
        narginchk(2,4);

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

% Does user want to specify threshold?
threshold = [];
if (3<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        threshold = temp;
    end
end


% Does user want to specify fig_num?
flag_do_plots = 0;
if (0==flag_max_speed) &&  (4<= nargin)
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
flag_is_similar = true;

NtestPoints = length(test_points_XY(:,1));
minimum_distance_to_each_point       = nan(NtestPoints,1);
indicies_of_nearest_reference_points = nan(NtestPoints,1);


for ith_point = 1:NtestPoints
    current_test_point = test_points_XY(ith_point,:);
    distances_to_reference_points = sum((reference_points_XY - current_test_point).^2,2).^0.5;
    [min_dist, min_index] = min(distances_to_reference_points);
    minimum_distance_to_each_point(ith_point,1) = min_dist;
    indicies_of_nearest_reference_points(ith_point,1) = min_index;
end

% Make sure all the points "pass"
if ~isempty(threshold) && any(minimum_distance_to_each_point>threshold)
    flag_is_similar = false;
end

% find mean, max, and std
mean_error    = mean(minimum_distance_to_each_point,"all","omitmissing");
max_error     = max(minimum_distance_to_each_point,[],"all","omitmissing");
std_dev_error = std(minimum_distance_to_each_point, 0,"all","omitmissing");

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
    xlabel('X [meters]');
    ylabel('Y [meters]')

    % Get the colormap
    error_colormap = turbo;
    Ncolors = length(error_colormap(:,1));

    % Plot the inputs
    plot(reference_points_XY(:,1),reference_points_XY(:,2),'.','Color',[0 1 0],'MarkerSize',20);
    plot(test_points_XY(:,1),test_points_XY(:,2),'.','Color',[0 0 1 ],'MarkerSize',20);

    % Sort the points by distance
    [sorted_distances,sorted_indicies] = sort(minimum_distance_to_each_point);

    % Figure out which max error to use
    if isempty(threshold)
        max_allowable_error = sorted_distances(end);
    else
        max_allowable_error = threshold;
    end

    for ith_test_point = 1:NtestPoints
        current_index = sorted_indicies(ith_test_point);
        current_distance = minimum_distance_to_each_point(current_index);
        current_test_point = test_points_XY(current_index,:);
        closest_ref_point  = reference_points_XY(indicies_of_nearest_reference_points(current_index),:);
        line_points = [current_test_point; closest_ref_point];

        % Get color
        error_ratio = min(current_distance/max_allowable_error,1);
        current_color_index = round(error_ratio*(Ncolors-1))+1;
        current_color = error_colormap(current_color_index,:);

        plot(line_points(:,1),line_points(:,2),'-','LineWidth',3,'Color',current_color);
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
