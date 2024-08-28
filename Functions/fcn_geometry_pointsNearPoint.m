function nearbyIndicies = fcn_geometry_pointsNearPoint(anchorPoint, pointsToSearch, searchRadius, varargin)
%fcn_geometry_pointsNearPoint  finds points near a given point
%
% Given a 2D or 3D list of pointsToSearch and an anchorPoint, finds the
% indicies of pointsToSearch that are within searchRadius of distance to
% the anchorPoint. The output is a [Mx1] vector where M is the number of
% indicies, sorted in the same order as pointsToSearch, that are within the
% searchRadius of anchorPoint. If no points are within the range, it
% returns an empty vector.
%
% The method to use is to calculate the distance from the anchorPoint to
% all other points. Note: this is a very slow process.
%
% FORMAT:
%
%       nearbyIndicies = fcn_geometry_pointsNearPoint(anchorPoint, pointsToSearch, searchRadius, (fig_num))
%
% INPUTS:
%
%      anchorPoint: the point as a [1x2] or [1x3] vector of [X Y (Z)] to
%      search around
%
%      pointsToSearch: a [Nx2] or [Nx3] vector of points to check
%
%      searchRadius: a [1x1] scalar of [searchRadius] where searchRadius is
%      the query distance that determines search criteria for "nearby", in
%      meters
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      nearbyIndicies: the indicies of the pointsToSearch that are "nearby"
%      the anchorPoint
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%
%       script_test_fcn_geometry_pointsNearPoint
%
% This function was written on 2024_08_26 by Sean Brennan
% Questions or comments? sbrennan@psu.edu

% Revision History
% 2024_08_26 S. Brennan
% -- started writing function

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{end},-1))
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_PLOTCV2X_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_PLOTCV2X_FLAG_CHECK_INPUTS");
    MATLABFLAG_PLOTCV2X_FLAG_DO_DEBUG = getenv("MATLABFLAG_PLOTCV2X_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_PLOTCV2X_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_PLOTCV2X_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_PLOTCV2X_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_PLOTCV2X_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 999978; %#ok<NASGU>
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

if 0 == flag_max_speed
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(3,4);

    end
end

% Does user want to specify fig_num?
flag_do_plots = 0;
fig_num = []; % Initialize the figure number to be empty
if (0==flag_max_speed) && (4 <= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end


%% Write main code for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the output
Npoints = length(pointsToSearch(:,1));

% Loop through all the points, finding qualifying agreement indicies
allDistances = real(sum((ones(Npoints,1)*anchorPoint - pointsToSearch).^2,2).^0.5);
nearbyIndicies = find(allDistances<=searchRadius);


%% Any debugging?
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

% before opeaning up a figure, lets start to capture the frames for an
% animation if the user has entered a name for the mov file
if flag_do_plots == 1

    % Prepare the figure
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end

    hold on;
    axis equal
    grid minor;
    xlabel('X [meters]');
    ylabel('Y [meters]')

    pointDimension = length(pointsToSearch(1,:));
    if 2==pointDimension
        % Plot the anchorPoint
        plot(anchorPoint(:,1), anchorPoint(:,2), 'b.','MarkerSize',20);

        % Plot the pointsToSearch
        plot(pointsToSearch(:,1),pointsToSearch(:,2),'k.','MarkerSize',5);

        % Plot the region
        anglesAlongRange = linspace(0,2*pi,100)';
        regionBoundary = anchorPoint + searchRadius*[cos(anglesAlongRange) sin(anglesAlongRange)];
        regionShape = polyshape(regionBoundary(:,1),regionBoundary(:,2),'Simplify',false,'KeepCollinearPoints',true);
        current_color = [0 1 0];
        plot(regionShape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',1,'EdgeAlpha',0);

        % Plot the points in agreement
        plot(pointsToSearch(nearbyIndicies,1), pointsToSearch(nearbyIndicies,2), 'go','MarkerSize',10,'LineWidth',2);
        
        % Make axis slightly larger?
        if flag_rescale_axis
            temp = axis;
            %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
            axis_range_x = temp(2)-temp(1);
            axis_range_y = temp(4)-temp(3);
            percent_larger = 0.3;
            axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
        end

    elseif 3==pointDimension
        % Plot the anchorPoint
        plot3(anchorPoint(:,1), anchorPoint(:,2), anchorPoint(:,3), 'b.','MarkerSize',20);

        % Plot the pointsToSearch
        plot3(pointsToSearch(:,1), pointsToSearch(:,2), pointsToSearch(:,3), 'k.','MarkerSize',5);

         % Plot the points in agreement
        plot3(pointsToSearch(nearbyIndicies,1), pointsToSearch(nearbyIndicies,2), pointsToSearch(nearbyIndicies,3),  'go','MarkerSize',10,'LineWidth',2);
        
        view(3);

        % Make axis slightly larger?
        if flag_rescale_axis
            temp = axis;
            %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
            axis_range_x = temp(2)-temp(1);
            axis_range_y = temp(4)-temp(3);
            axis_range_z = temp(6)-temp(5);
            percent_larger = 0.3;
            axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y,  temp(5)-percent_larger*axis_range_z, temp(6)+percent_larger*axis_range_z]);
        end

    else
        warning('on','backtrace');
        warning('Unknown dimension encountered when plotting.')
        error('Unknown dimension encountered when plotting.')
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
