function intersectionPoints = fcn_geometry_findBreakpointsOfAdjacentSegments(curveStartPoint, sortedHoughSegmentEndPoints, curveEndPoint, varargin)
%% fcn_geometry_findBreakpointsOfAdjacentSegments
%
% This function takes the start point of the curve, end point of the curve
% and all the endpoints of segments as the inputs and outputs the
% intersection points to fit a continuous curve to all the hough segment
% domains.
% 
% FORMAT: 
%
% fcn_geometry_findIntersectionPointsOfAdjacentSegments(curveStartPoint, curveEndPoint, sortedHoughSegmentEndPoints, fig_num)
% 
% INPUTS:
%
% curveStartPoint: The start point of the curve
%
% curveEndPoint: The end point of the curve
%
% sortedHoughSegmentEndPoints: The endpoints of the sorted Regression
% Domains "or" the second output of "fcn_geometry_sortRegressionDomains"
% function
%
% (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
% intersectionPoints: This is a [length(sortedHoughSegmentEndPoints(:,1))-2
% X 2] intersection points matrix of all the adjacent segments
% 
%
% DEPENDENCIES:
%
%   NONE
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_findIntersectionPointsOfAdjacentSegments
% for a full test suite.
%
% This function was written on 2024_02_29 by Aneesh Batchu
% Questions or comments? abb6486@psu.edu or sbrennan@psu.edu

% Revision History
% 2024_03_02 
% -- wrote the code - Aneesh Batchu
% 2024_03_11 
% -- Functionalized the code - Aneesh Batchu
% 2024_03_14
% -- Replaced this function "fcn_geometry_findIntersectionOfSegments" with
% "fcn_Path_findProjectionHitOntoPath"

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
        narginchk(3,4);

        % Check the tolerance input is a positive single number
        % fcn_DebugTools_checkInputsToFunctions(tolerance, 'positive_1column_of_numbers',1);

    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (4<= nargin)
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


if isempty(curveStartPoint) && isempty(curveEndPoint)
    curveStartPoint = sortedHoughSegmentEndPoints(1,1:2);
    curveEndPoint = sortedHoughSegmentEndPoints(end,3:4);
end

% Pre-allocating the intersectionPoints matrix for speed calculations
intersectionPoints = zeros(length(sortedHoughSegmentEndPoints(:,1))-2,2);

% This loop calls the "fcn_geometry_findIntersectionOfSegments" function to
% find the break points (intersection points) of the adjacent segments
for i = 1:length(sortedHoughSegmentEndPoints(:,1)) - 1
    % fprintf(1,'Intersection point: \n');
    wall_start = sortedHoughSegmentEndPoints(i,1:2);
    wall_end   = sortedHoughSegmentEndPoints(i,3:4);
    sensor_vector_start = sortedHoughSegmentEndPoints(i+1,1:2);
    sensor_vector_end   = sortedHoughSegmentEndPoints(i+1,3:4);
    fig_debugging = [];
    flag_search_type =4;
    % [~,location] = ...
    %     fcn_geometry_findIntersectionOfSegments(...
    %     wall_start, wall_end,sensor_vector_start,sensor_vector_end,...
    %     flag_search_type,fig_debugging);
    [~,location] = fcn_Path_findProjectionHitOntoPath(...
        [wall_start; wall_end],sensor_vector_start,sensor_vector_end,...
        flag_search_type,fig_debugging);

    % disp(location);
    intersectionPoints(i,:) = location;
    
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

    plotPoints = [curveStartPoint; intersectionPoints; curveEndPoint];
    

    % figure(fig_debugging)
    % hold on
    % plot(plotPoints(:,1), plotPoints(:,2), '--', 'LineWidth',5, 'Color', [1 1 0]);

    figure(fig_num)
    grid on
    grid minor
    box on
    hold on
    plot(plotPoints(:,1), plotPoints(:,2), '-', 'LineWidth',2, 'Color', [0 1 0], 'DisplayName','Fitted Curve');
    plot(intersectionPoints(:,1), intersectionPoints(:,2), '.', 'MarkerSize',40, 'Color', [1 0 0],'DisplayName','Intersection Points');
    plot(intersectionPoints(:,1), intersectionPoints(:,2), 'o', 'MarkerSize',30, 'LineWidth', 4, 'Color', [0 1 1],'DisplayName','Intersection Points');
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

