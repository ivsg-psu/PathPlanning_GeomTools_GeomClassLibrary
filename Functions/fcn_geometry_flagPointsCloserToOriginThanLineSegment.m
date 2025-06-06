function [point_flags] = ...
    fcn_geometry_flagPointsCloserToOriginThanLineSegment(...
    segment_points, ...
    test_points, ...
    varargin)
% fcn_geometry_flagPointsCloserToOriginThanLineSegment
% Genrates a column vector that is 1 if associated point is closer to
% origin (inside) a line segment, 0 otherwise. Points must be strictly
% inside to count (e.g. cannot be ON the line segment).
%
% The method used is to fit a line to the line segment to calculate the
% slope and intercept. It then applies the same slope to the test points,
% and calculates their intercepts. If the test-point intercepts are closer
% to the origin than the line segment, then they are "within" the segment.
% For vertical line segments, a special test is done to see if the x-axis
% coordinate is closer to the orgin than the line segment.
% 
% FORMAT: 
%
% [point_flags] = ...
%     fcn_geometry_flagPointsCloserToOriginThanLineSegment(...
%     segment_points, ...
%     test_points, ...
%     (fig_num))
%
% INPUTS:
%      segment_points: a 2x2 vector where the first row is the [x y]
%      location of the start of the segment, the second row is the [x y]
%      location of the end of the segment
%
%      test_points: a Nx2 vector where N is the number of points that are
%      being checked (N must be at least 1);
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%      point_flags: a vector (Nx1) representing whether associated points
%      are closer to the origin than the line segment
%
% EXAMPLES:
%      
%      % BASIC example
%      segment_points = [2 3; 4 5];
%      test_points  = [1 1];
%      [slope,intercept] = fcn_geometry_find_slope_intercept_from_N_points(points)
% 
% See the script: script_test_fcn_geometry_flagPointsCloserThanLineSegment
% for a full test suite.
%
% This function was written on 2020_06_25 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2020_06_25 
% -- wrote the code
% 2021_05_26
% -- Improved the comments, prepped for geometry class
% 2021_05_31
% -- Fixed bug where points beyond the origin were also included
% incorrectly
% 2024_01_17 - Aneesh Batchu
% -- added max speed options
% 2024_04_15 - S. Brennan
% -- added plot handles to legend entries to keep from causing warnings 


%% Debugging and Input checks
% flag_check_inputs = 1; % Set equal to 1 to check the input arguments
% flag_do_plot = 0;      % Set equal to 1 for plotting
% flag_do_debug = 0;     % Set equal to 1 for debugging

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

if 0==flag_max_speed
    % Should we check the inputs?
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(2,3);

        % Check the segment_points input to have at least 2 or more rows
        fcn_DebugTools_checkInputsToFunctions(...
            segment_points, '2column_of_numbers',[2 3]);

        % Check the test_points input to have at least 2 or more rows
        fcn_DebugTools_checkInputsToFunctions(...
            test_points, '2column_of_numbers');

    end
end

% % Does user want to show the plots?
% if 3 == nargin
%     fig_num = varargin{end};
%     figure(fig_num);
%     flag_do_plot = 1;
% else
%     if flag_do_debug
%         fig = figure; 
%         fig_num = fig.Number;
%         flag_do_plot = 1;
%     end
% end

% Does user want to show the plots?
flag_do_plot = 0;
if (0==flag_max_speed) && (3 == nargin) 
    temp = varargin{1};
    if ~isempty(temp)
        fig_num = temp;
        figure(fig_num);
        flag_do_plot = 1;
    end
else
    if flag_do_debug
        fig = figure; 
        fig_num = fig.Number;
        flag_do_plot = 1;
    end
end

%% Main code starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% The code below just calculates the slope and intercept of the line
% segment created by the segment points:
%     y = m*x + b
% then calculates the intercept related to the line that passes through
% each test point. If the intercept is closer to the origin, then it is
% inside the line segment. We have to consider a special case for when the
% line is vertial

% Fill in slope and intercept from the segment points
[slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(segment_points);

% Check to see if the line is vertical? (different result for each)
% If line is vertical, just compare x values to the first x-value in the
% line segment
if isinf(slope)  
    point_flags = ...
        (abs(test_points(:,1)) < abs(segment_points(1,1))) & ...
        (sign(test_points(:,1))== sign(segment_points(1,1)));
    
else  % Not infinite, so the result will be an ordinary line
    % the b value is just b = y - m*x. We take the absolute value so that
    % this approach will work for negative intercepts as well.
    point_flags = ...
        (abs(test_points(:,2) - slope*test_points(:,1)) < abs(intercept)) & ...
        (sign(test_points(:,2) - slope*test_points(:,1)) == sign(intercept));
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
if flag_do_plot
    figure(fig_num);
    hold on;
    grid on;
    axis equal;
    
    % START BY PLOTTING THE LINE SEGMENT
    % Create an x vector of 100 points connecting lowest and highest points
    % in x, and plot the line fit
    x = linspace(min(segment_points(:,1)),max(segment_points(:,1)),100)';
    y = x*slope + intercept;
    if isinf(slope)  % The result is a vertical line        
        y = linspace(min(segment_points(:,2)),max(segment_points(:,2)),100)';
    end    
    h_lineSegment = plot(x,y,'b-');
    
    % PLOT THE ORIGIN LINE
    y_origin = x*slope;
    if isinf(slope)  % The result is a vertical line
        x = x*0;
        y_origin = linspace(min(segment_points(:,2)),max(segment_points(:,2)),100)';
    end    
    h_originCutoff = plot(x,y_origin,'g-');
    
    
    % NOW PLOT THE POINTS
    h_testPoints = plot(test_points(:,1),test_points(:,2),'k.');
    
    % Circle inside points with green
    inside_points = test_points(point_flags>0,:);
    if isempty(inside_points)
        inside_points = [nan nan];
    end
    h_inside = plot(inside_points(:,1),inside_points(:,2),'go');

    % Circle outside points with red
    outside_points = test_points(point_flags==0,:);
    if isempty(outside_points)
        outside_points = [nan nan];
    end
    h_outside = plot(outside_points(:,1),outside_points(:,2),'ro');
    
    legend([h_lineSegment, h_originCutoff,  h_testPoints, h_inside, h_outside], 'Line segment','Origin cut-off','Test points','Inside','Outside');

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends function

% 
% function [slope,intercept] = fcn_geometry_find_slope_intercept_from_N_points(points)
% % fcn_geometry_find_slope_intercept_from_N_points
% % Finds the slope and intercept of a line connecting two points
% % NOTE: this was exerpted from a much more robust stand-alone function of
% % the same name. See that function for details.
% 
% Npoints = length(points(:,1));
% 
% % Fill in X and Y
% X = points(:,1);
% Y = points(:,2);
% 
% % Check to see if the result is going to be singular. This happens if all
% % the x values are the same, e.g. the line is vertical
% if all(X == X(1))  % Are all the x values the same?
%     slope = inf;
%     intercept = inf;
% else  % The result will be an ordinary line
%     result = [X ones(Npoints,1)]\Y;    
%     slope = result(1,1);
%     intercept = result(2,1);
% end
% end

