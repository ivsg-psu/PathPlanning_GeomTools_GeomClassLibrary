
function [point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points, varargin)
% fcn_geometry_flagPointsCloserToOriginThanLineSegment
% Genrates a column vector that is 1 if associated point is closer to
% origin (inside) a line segment, 0 otherwise. Points must be strictly
% inside to count (e.g. cannot be ON the line segment)
% 
% Format: 
% [point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points, test_points, varargin)
%
% INPUTS:
%      segment_points: a 2x2 vector where the first row is the [x y]
%      location of the start of the segment, the second row is the [x y]
%      location of the end of the segment
%
%      test_points: a Nx2 vector where N is the number of points that are
%      being checked (N must be at least 1);
%
% OUTPUTS:
%      point_flags: a vector (Nx1) representing whether associated points
%      are closer to the origin than the line segment
%
% Examples:
%      
%      % BASIC example
%      segment_points = [2 3; 4 5];
%      test_points  = [1 1];
%      [slope,intercept] = fcn_geometry_find_slope_intercept_from_N_points(points)
% 
% See the script: script_test_fcn_geometry_flagPointsCloserToOriginThanLineSegment
% for a full test suite.
%
% This function was written on 2020_06_25 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2020_06_25 - wrote the code
%

flag_do_debug = 1; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

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
% Are the input vectors the right shape?
Npoints = length(test_points(:,1));

if flag_check_inputs == 1
    % Are there the right number of inputs?
    if nargin < 2 || nargin > 3
        error('Incorrect number of input arguments')
    end
    
    if Npoints<1
        error('The points vector must have at least 1 row, with each row representing a different (x y) point');
    end
    
    if length(test_points(1,:))~=2
        error('The test_points vector must have 2 columns, with column 1 representing the x portions of the points, column 2 representing the y portions.');
    end
    
    if ~isequal(size(segment_points),[2 2])
       error('The segment_points vector must have 2 rows and 2 columns, with column 1 representing the x portions of the points, column 2 representing the y portions. The first row is for the start point, the second is for the end point of the segment.');
    end 
    
end

% Does user want to show the plots?
if 3 == nargin
    fig_num = varargin{1};
    figure(fig_num);
    flag_do_debug = 1;
else
    if flag_do_debug
        fig = figure; 
        fig_num = fig.Number;
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
 
% The code below just calculates the slope and intercept of the line segment:
%     y = m*x + b
% then calculates the intercept related to the line that passes through each test point.
% If the intercept is closer to the origin, then it is inside the line
% segment. We have to consider a special case for when the line is vertial

% Fill in slope and intercept
[slope,intercept] = fcn_geometry_find_slope_intercept_from_N_points(segment_points);

% Check to see if the line is vertical? (different result for each)
if isinf(slope)  % If line is vertical, just compare x values to the first x-value in the line segment
    point_flags = abs(test_points(:,1)) < abs(segment_points(1,1));
else  % The result will be an ordinary line
    % the b value is just b = y - m*x. We take the absolute value so that
    % this approach will work for negative intercepts as well.
    point_flags = abs(test_points(:,2) - slope*test_points(:,1))    < abs(intercept);
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
if flag_do_debug
    figure(fig_num);
    hold on;
    grid on;
    
    % START BY PLOTTING THE LINE SEGMENT
    % Create an x vector of 100 points connecting lowest and highest points
    % in x, and plot the line fit
    x = linspace(min(segment_points(:,1)),max(segment_points(:,1)),100)';
    y = x*slope + intercept;
    if isinf(slope)  % The result is a vertical line
        y = linspace(min(segment_points(:,2)),max(segment_points(:,2)),100)';
    end    
    plot(x,y,'b-');
    
    % NOW PLOT THE POINTS
    inside_points = test_points(point_flags>0,:);
    plot(inside_points(:,1),inside_points(:,2),'g.');

    outside_points = test_points(point_flags==0,:);
    plot(outside_points(:,1),outside_points(:,2),'r.');

end
end

function [slope,intercept] = fcn_geometry_find_slope_intercept_from_N_points(points)
% fcn_geometry_find_slope_intercept_from_N_points
% Finds the slope and intercept of a line connecting two points
% NOTE: this was exerpted from a much more robust stand-alone function of
% the same name. See that function for details.

Npoints = length(points(:,1));

% Fill in X and Y
X = points(:,1);
Y = points(:,2);

% Check to see if the result is going to be singular. This happens if all
% the x values are the same, e.g. the line is vertical
if all(X == X(1))  % Are all the x values the same?
    slope = inf;
    intercept = inf;
else  % The result will be an ordinary line
    result = [X ones(Npoints,1)]\Y;    
    slope = result(1,1);
    intercept = result(2,1);
end
end

