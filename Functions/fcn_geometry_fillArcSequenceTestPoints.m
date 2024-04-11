function [test_points, circleCenters, trueStartPointsOfArcs, arcStartIndicies, namedCurveTypes] = fcn_geometry_fillArcSequenceTestPoints(arc_pattern, M, sigma, varargin)
%% fcn_geometry_fillArcSequenceTestPoints
% given a Ax2 matrix containing [curvature station_length] in each row,
% where A are the number of interconnected arcs to form. This function then
% creates a sequence of test points that follows each respective arc for
% the specified station length, in sequence
%
% [test_points] = fcn_geometry_fillArcSequenceTestPoints(seed_points, M, sigma)
%
% INPUTS:
%      arc_pattern: a Ax2 matrix containing [curvature station_length] in
%      each row, where A is the number of arcs.
%
%      M: the number of test points to generate per unit
%      distance.
%
%      sigma: athe standard deviation in points
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      test_points: a list of test points used to test regression fitting
%
%      circleCenters: the [x y] location of each circle
%      center
% 
%      trueStartPointsOfArcs: the true (not noise injected) [x y] location where each arc starts.
%
%      arcStartIndicies: the indicies where each segment started
%
%      namedCurveTypes: a cell array of names indicating the type of curve
%      for each test point sequence. Names include:
%
%           'line': for line segments
%
%           'arc': for arc segments
% 
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_circleCenterFrom3Points
%      fcn_geometry_arcDirectionFrom3Points
%      fcn_geometry_calcUnitVector
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_fillArcSequenceTestPoints
% for a full test suite.
%
% This function was written on 2024_03_31 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2024_03_31 - S. Brennan
% -- wrote the code


%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{end},-1))
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
    debug_fig_num = 34838;
    figure(debug_fig_num);
    clf;
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

if (0==flag_max_speed)
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(3,4);

        % Check the points input to be length greater than or equal to 2
        fcn_DebugTools_checkInputsToFunctions(...
            arc_pattern, '2column_of_numbers');

    end
end

% Does user want to show the plots?
flag_do_plots = 0;
if (4 == nargin) && (0==flag_max_speed)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve for the circle
N_segments = length(arc_pattern(:,1));

% Find radius and arc angles for each input
circleRadius = 1./arc_pattern(:,1);

% Find if the arcs are counterclockwise or clockwise
is_counterClockwise = 2*(circleRadius>0)-1;

% Save the circle centers that we encounter
circleCenters = nan(N_segments,2);
trueStartPointsOfArcs = nan(N_segments,2);
arcStartIndicies = nan(N_segments,1);

current_angle = 0;
current_position = [0 0];
test_points = current_position;
for ith_segment = 1:N_segments
    trueStartPointsOfArcs(ith_segment,:) = current_position;
    arcStartIndicies(ith_segment) = length(test_points(:,1));
    total_arc_length = arc_pattern(ith_segment,2);

    % Find current radius
    current_radius = abs(circleRadius(ith_segment));

    % Find the circle center
    if ~isinf(current_radius)
        namedCurveTypes{ith_segment} = 'arc';
        current_circle_center = current_position + current_radius*[cos(current_angle) sin(current_angle)]*[0 is_counterClockwise(ith_segment,1); -is_counterClockwise(ith_segment,1) 0];
        circleCenters(ith_segment,:) = current_circle_center;

        % Find the number points
        projection_distances = (0:(1/M):total_arc_length)';

        N_points = length(projection_distances);

        % Convert from distances to angles, and then angles to points
        angles = current_angle + is_counterClockwise(ith_segment)*projection_distances./current_radius;
        if 1==is_counterClockwise(ith_segment,1)
            unit_radial_vectors = [cos(angles-pi/2) sin(angles-pi/2)];
        else
            unit_radial_vectors = [cos(angles+pi/2) sin(angles+pi/2)];
        end
        orthogonal_distances = randn(N_points,1)*sigma + ones(N_points,1)*current_radius;
        perturbed_points = ones(N_points,1)*current_circle_center + orthogonal_distances.*unit_radial_vectors;

        % Repeat for the last point, which may not be included in the list
        % above, and should also NOT have any noise added
        final_angle = current_angle + is_counterClockwise(ith_segment)*total_arc_length./current_radius;
        if 1==is_counterClockwise(ith_segment,1)
            unit_radial_vectors = [cos(final_angle-pi/2) sin(final_angle-pi/2)];
        else
            unit_radial_vectors = [cos(final_angle+pi/2) sin(final_angle+pi/2)];
        end
        final_point = current_circle_center + current_radius.*unit_radial_vectors;
    else
        namedCurveTypes{ith_segment} = 'line';
        circleCenters(ith_segment,:) = [nan nan];
        projection_distances = (0:(1/M):total_arc_length)';
        N_points = length(projection_distances);
        unit_tangential_vector = [cos(current_angle) sin(current_angle)];
        unit_radial_vector = [cos(current_angle+pi/2) sin(current_angle+pi/2)];
        orthogonal_distances = randn(N_points,1)*sigma;
        perturbed_points = current_position + projection_distances*unit_tangential_vector + orthogonal_distances.*unit_radial_vector;
        final_point = current_position + total_arc_length.*unit_tangential_vector;
    end

    % Convert to test points
    test_points = [test_points; perturbed_points]; %#ok<AGROW>

    % Save where we are for the next loop
    current_angle = current_angle + is_counterClockwise(ith_segment)*total_arc_length./current_radius;
    current_position = final_point;

    if flag_do_debug
        figure(34345);
        hold on;
        grid on;
        axis equal

        plot(trueStartPointsOfArcs(:,1), trueStartPointsOfArcs(:,2), 'g.','MarkerSize',30);
        plot(circleCenters(ith_segment,1), circleCenters(ith_segment,2), 'r+','MarkerSize',30);
        plot(test_points(:,1), test_points(:,2), 'k.','MarkerSize',20);
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
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end       
    axis equal;
    hold on;
    grid on;

    % Get the color ordering?
    try
        color_ordering = orderedcolors('gem12');
    catch
        color_ordering = colororder;
    end

    N_colors = length(color_ordering(:,1));


    % Plot the circle centers
    plot(circleCenters(:,1), circleCenters(:,2), 'r+','MarkerSize',30);

    % Plot the results
    plot(test_points(:,1),test_points(:,2),'b.','MarkerSize',10);

    % Plot the segments in different colors
    for ith_row = 1:length(arc_pattern(:,1))
        current_color = color_ordering(mod(ith_row,N_colors)+1,:);
        plot(trueStartPointsOfArcs(ith_row,1),trueStartPointsOfArcs(ith_row,2),'.','MarkerSize',50,'Color',current_color);
        if ith_row < length(arc_pattern(:,1))
            range = (arcStartIndicies(ith_row):arcStartIndicies(ith_row+1));
        else
            range = (arcStartIndicies(ith_row):length(test_points(:,1)));
        end
        plot(test_points(range,1),test_points(range,2),'o','MarkerSize',5,'Color',current_color);
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







