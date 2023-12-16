function [point_flags] = ...
    fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(...
    segment_points, ...
    test_points, ...
    varargin)
% fcn_geometry_flagPointsFurtherFromOriginThanLineSegment
% Genrates a column vector that is 1 if associated point is further from
% origin (outside) a line segment, 0 otherwise. Points must be strictly
% outside to count (e.g. cannot be ON the line segment).
%
% The method used is to fit a line to the line segment to calculate the
% slope and intercept. It then applies the same slope to the test points,
% and calculates their intercepts. If the test-point intercepts are closer
% to the origin than the line segment, then they are "within" the segment.
% For vertical line segments, a special test is done to see if the x-axis
% coordinate is closer to the orgin than the line segment.
%
% Note: the points must be a tolerance beyond the line to be "further", and
% this tolerance is a factor of the numerical precision of the environment,
% eps, typically 10*eps where eps = 2.2204e-16.
% 
% FORMAT: 
%
% [point_flags] = ...
%     fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(...
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
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%      point_flags: a vector (Nx1) representing whether associated points
%      are closer to the origin than the line segment
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_flagPointsFurtherFromOrigin.m
% for a full test suite.
%
% This function was written on 2021_05_31 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2021_05_31
% -- Wrote the code via mods to fcn_geometry_flagPointsCloserToOriginThanLineSegment

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

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

% Should we check the inputs?
if flag_check_inputs == 1
    % Are there the right number of inputs?
    narginchk(2,3);
    
    % Check the segment_points input to have at least 2 or more rows
    fcn_geometry_checkInputsToFunctions(...
        segment_points, '2column_of_numbers',[2 3]);
    
    % Check the test_points input to have at least 2 or more rows
    fcn_geometry_checkInputsToFunctions(...
        test_points, '2column_of_numbers');
    
end

% Does user want to show the plots?
if 3 == nargin
    fig_num = varargin{end};
    figure(fig_num);
    flag_do_plot = 1;
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
        (abs(test_points(:,1)) > (abs(segment_points(1,1))+10*eps)) & ...
        (sign(test_points(:,1))== sign(segment_points(1,1)));
    
else  % Not infinite, so the result will be an ordinary line
    % the b value is just b = y - m*x. We take the absolute value so that
    % this approach will work for negative intercepts as well.
    point_flags = ...
        (abs(test_points(:,2) - slope*test_points(:,1)) > (abs(intercept)+10*eps)) & ...
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
    plot(x,y,'b-');
    
    % PLOT THE ORIGIN LINE
    y_origin = x*slope;
    if isinf(slope)  % The result is a vertical line
        x = x*0;
        y_origin = linspace(min(segment_points(:,2)),max(segment_points(:,2)),100)';
    end    
    plot(x,y_origin,'g-');
    
    
    % NOW PLOT THE POINTS
    plot(test_points(:,1),test_points(:,2),'k.');
    
    % Circle inside points with red
    inside_points = test_points(point_flags==0,:);
    plot(inside_points(:,1),inside_points(:,2),'ro');

    % Circle outside points with green
    outside_points = test_points(point_flags>0,:);
    plot(outside_points(:,1),outside_points(:,2),'go');
    
    legend('Line segment','Origin cut-off','Test points','Inside','Outside');

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

