function [agreement_indicies, flag_is_a_circle, start_angle_in_radians, end_angle_in_radians] = ...
    fcn_geometry_findAgreementsOfPointsToArc(points, base_point_index, circleCenter, circleRadius, transverse_tolerance, varargin) 
% fcn_geometry_findAgreementsOfPointsToArc
% Given a set of XY points, a base_point_index where an arc is rooted, a
% circleCenter, and a circleRadius, finds the indicies of the points that
% are within a transverse_tolerance distance away from the circle radius,
% while keeping within station_tolerance distance from each point within
% the cluster centered at the base_point_index.
%
% Format:
% agreement_indicies = fcn_geometry_findAgreementsOfPointsToArc(points, circleCenter, circleRadius, transverse_tolerance, (fig_num))
%
% INPUTS:
%
%      points: a Nx2 vector where N is the number of points, but at least 2
%      rows.
%
%      base_point_index: the index of the point to use as a base point,
%      among the input points.
%
%      circleCenter: a 1x2 vector that is a the [x y] coordinates of the
%      circle center to test.
%
%      circleRadius: the radius of the circle to test
%
%      transverse_tolerance: the orthogonal distance between the points and
%      the linear vector fit that indicate whether a point "belongs" to the
%      fit. A point belongs to the fit if the transverse distance is less
%      than or equal to the transverse_tolerance
%
%      (OPTIONAL INPUTS)
%
%      station_tolerance: the projection distance between the points in a
%      curve fit, along the direction of the line, that indicate whether a
%      point "belongs" to the circle fit (if distance is less than or equal
%      to the tolerance), or is "outside" the fit (if distance is greater
%      than the tolerance). If left empty, then a circle fit is performed
%      using only the transverse_tolerance.
%
%      flag_force_circle_fit: specify that only circle fits are allowed.
%      Can be set to:
%
%          0: search for both arcs and circles (default)
%
%          1: search only for circle fits. If station tolerance is given,
%          then a circle will only be returned if all points meet both the
%          transverse and station tolerances. Note: to force a circle fit
%          irregardless of station tolerance, set station_tolerance to an
%          empty value, e.g. station_tolerance = [];

%
%      threshold_to_check_arc: a positive integer that represents the
%      minimum number of points that must be in circle agreement before an
%      arc is checked. The circle testing is very fast; however, the arc
%      testing is slow. If the threshold is known beforehand, then entering
%      this threshold significantly speeds up the analysis.
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
%      flag_is_a_circle: 1 if the result forms a circle, 0 otherwise
%  
%      start_angle_in_radians: the angle that the arc starts at, between 0
%      and 2pi
% 
%      end_angle_in_radians: the angle that the arc ends at, between 0
%      and 2pi
% 
% DEPENDENCIES:
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_findAgreementsOfPointsToCircle
%      fcn_geometry_findArcAgreementIndicies
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_findAgreementsOfPointsToArc
% for a full test suite.

% This function was written on 2024_01_15 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2023_01_15/16 - S. Brennan
% -- wrote the code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==9 && isequal(varargin{end},-1))
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
        narginchk(5,9);

        % Check the points input to be length greater than or equal to 2
        fcn_DebugTools_checkInputsToFunctions(...
            points, '2column_of_numbers',[2 3]);

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

        % Check the circleCenter input to be [1x2]
        fcn_DebugTools_checkInputsToFunctions(...
            circleCenter, '2column_of_numbers',[1 1]);

        % Check the circleRadius input to be [1x1] positive number
        fcn_DebugTools_checkInputsToFunctions(...
            circleRadius, 'positive_1column_of_numbers',1);

        % Check the transverse_tolerance input is a positive single number
        fcn_DebugTools_checkInputsToFunctions(transverse_tolerance, 'positive_1column_of_numbers',1);

    end
end

% Does user want to specify threshold_to_check_arc?
station_tolerance = [];
if 6<= nargin
    temp = varargin{1};
    if ~isempty(temp)
        station_tolerance = temp;
    end
end

% Does user want to specify threshold_to_check_arc?
flag_force_circle_fit = 0;
if 7<= nargin
    temp = varargin{2};
    if ~isempty(temp)
        flag_force_circle_fit = temp;
    end
end


% Does user want to specify threshold_to_check_arc?
threshold_to_check_arc = 1;
if 8<= nargin
    temp = varargin{3};
    if ~isempty(temp)
        threshold_to_check_arc = temp;
    end
end

% Does user want to specify fig_num?
flag_do_plots = 0;
if 0==flag_max_speed && 9<=nargin
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
indicies_in_transverse_agreement = fcn_geometry_findAgreementsOfPointsToCircle(points, circleCenter, circleRadius, transverse_tolerance, -1); 

% Find the indicies in station agreement. Only do this if the transverse
% count is larger than current best count (this needs to be debugged)
flag_bad_fit_found = 0;
if length(indicies_in_transverse_agreement)>threshold_to_check_arc 
    if ~isempty(station_tolerance)
        % Grab only the points in radial agreement
        points_in_radial_agreement = points(indicies_in_transverse_agreement,:);

        % Find index of the source point
        index_source_point = find(indicies_in_transverse_agreement == base_point_index,1);

        % Was an index given that does not fit?
        if ~isempty(index_source_point)
            % Call a function to find which points are in station agreement
            % and output the results in sorted order
            [indicies_in_station_agreement, flag_is_a_circle, start_angle_in_radians, end_angle_in_radians] = ...
                fcn_geometry_findArcAgreementIndicies(points_in_radial_agreement, circleCenter, circleRadius, index_source_point, station_tolerance, -1);


            % Calculate agreement with both vectors
            % agreement_indicies_binary_form = indicies_in_transverse_agreement.*points_within_angle;
            agreement_indicies = indicies_in_transverse_agreement(indicies_in_station_agreement);
        else
            flag_bad_fit_found = 1;
        end

        if (flag_force_circle_fit==1) && (flag_is_a_circle==0)
            flag_bad_fit_found = 1;
        end
    else
        % Force a circle output because station is empty
        agreement_indicies = indicies_in_transverse_agreement;
        flag_is_a_circle = 1;
        start_angle_in_radians = 0;
        end_angle_in_radians = 2*pi;
    end
else
    flag_bad_fit_found = 1;
end

if 1==flag_bad_fit_found
    agreement_indicies = [];
    flag_is_a_circle = -1;
    start_angle_in_radians = nan;
    end_angle_in_radians = nan;
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

    xlabel('X [meters]');
    ylabel('Y [meters]')

    if flag_is_a_circle==1
        title('Circle fit');
    elseif flag_is_a_circle==0
        title('Arc fit');
    elseif flag_is_a_circle==-1
        title('Bad fit');
    end
    
    % Plot the input points
    plot(points(:,1),points(:,2),'k.','MarkerSize',30)

    % Plot the boundaries of the circle fit
    fcn_geometry_plotCircle(circleCenter, circleRadius, 'b-',fig_num); 
    fcn_geometry_plotCircle(circleCenter, circleRadius-transverse_tolerance, 'r-',fig_num); 
    fcn_geometry_plotCircle(circleCenter, circleRadius+transverse_tolerance, 'r-',fig_num); 
    plot(circleCenter(:,1),circleCenter(:,2),'b+','MarkerSize',15);

    if ~isempty(agreement_indicies)
        % Plot the radial agreement points
        plot(points(agreement_indicies,1),points(agreement_indicies,2),'r.','MarkerSize',15);

        % Label the points
        text(points(agreement_indicies(1),1),points(agreement_indicies(1),2),sprintf('Start angle: %.3f deg', start_angle_in_radians*180/pi));
        text(points(agreement_indicies(end),1),points(agreement_indicies(end),2),sprintf('End angle: %.3f deg', end_angle_in_radians*180/pi));
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



