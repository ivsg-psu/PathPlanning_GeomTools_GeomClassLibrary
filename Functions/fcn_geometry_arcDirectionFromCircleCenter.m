function is_counterClockwise = fcn_geometry_arcDirectionFromCircleCenter(points, circleCenter, varargin)
%% fcn_geometry_arcDirectionFromCircleCenter
% calculates whether or not a set of points make an arc in the clockwise or
% counter-clockwise diretion. The direction is calculated by performing a
% cross product between the vectors from each adjacent pair of points
% toward the circleCenter. The funtion returns the sign of the result.
%
% FORMAT:
%
% is_counterClockwise = fcn_geometry_arcDirectionFromCircleCenter(points, circleCenter,(fig_num))
%
% INPUTS:
%
%      points: a Nx2 vectors of [X Y] points, usually the ones used to
%      estimate the circle
%
%      circleCenter: a 1x2 vectors of the [X Y] position of the circle
%      center.
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      is_counterClockwise: a scalar containing 1 if the sum vote of points
%      creates a counter-clockwise arc (positive) or -1 if it would create
%      a clockwise arc (e.g. negative).
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      cross
%      sign
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_arcDirectionFromCircleCenter
% for a full test suite.
%
% This function was written on 2024_04_02 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_04_02 - sbrennan@psu.edu
% -- original write of the code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==3 && isequal(varargin{end},-1))
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
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(2,3);

        % % Check the points1 input
        % fcn_DebugTools_checkInputsToFunctions(...
        %     points, '2column_of_numbers');
        %
        % N_points = length(points(:,1));
        %
        % % Check the points2 input
        % fcn_DebugTools_checkInputsToFunctions(...
        %     points2, '2column_of_numbers',[N_points N_points]);
        %
        % % Check the points2 input
        % fcn_DebugTools_checkInputsToFunctions(...
        %     points2, '2column_of_numbers',[N_points N_points]);
    end
end

N_points = length(points(:,1));

% Does user want to show the plots?
flag_do_plots = 0;
if (0==flag_max_speed) && (3 <= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
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

% The code below is set up to be vectorized. The method is to create a
% vector from points1 to points2 and to points3. If the
projection_points1_to_points2 = points(2:end,:) - points(1:end-1,:);
projection_points_to_center   = circleCenter - points(1:end-1,:);

% Do the cross product from segment 1-2 to segment 1-3
cross_product = cross([projection_points1_to_points2 zeros(N_points-1,1)],[projection_points_to_center zeros(N_points-1,1)],2);

vote_clockwise = sign(cross_product(:,3))>0;
is_counterClockwise = sum(vote_clockwise)>(N_points/2);

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

    % Prep the figure
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end

    hold on;
    grid on;
    grid minor;
    axis equal;

    % Plot the inputs
    plot(points(:,1),points(:,2),'k.','MarkerSize',30);  % Plot points
    plot(circleCenter(1,1),circleCenter(1,2),'r+','MarkerSize',20);  % Plot circleCenter

    anti_votes = sign(cross_product(:,3))<0;
    quiver(points(vote_clockwise,1),points(vote_clockwise,2),circleCenter(1)-points(vote_clockwise,1),circleCenter(2)-points(vote_clockwise,2),0,'Color',[0 0 1]);
    quiver(points(anti_votes,1),points(anti_votes,2),circleCenter(1)-points(anti_votes,1),circleCenter(2)-points(anti_votes,2),0,'Color',[1 0 0]);


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

