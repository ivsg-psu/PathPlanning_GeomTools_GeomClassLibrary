
function [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,varargin)
% fcn_geometry_fitSlopeInterceptNPoints
% Finds the slope and intercept of a line connecting two points
% Format: 
% [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points)
%
% INPUTS:
%      points: a Nx2 vector where N is the number of points, but at least 2. 
%
% OUTPUTS:
%      slope: a scalar (1x1) representing the slope connecting the two
%      points
%      intercept: a scalar (1x1) representing the y-axis intercept of the
%      line fit
%
% DEPENDENCIES:
%
%      fcn_geometry_checkInputsToFunctions
%      fcn_geometry_plotCircle
%
% EXAMPLES:
%      
%      % BASIC example
%      points = [2 3; 4 5];
%      [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points)
% 
% See the script: script_test_fcn_geometry_fitSlopeInterceptNPoints
% for a full test suite.
%
% NOTE: This function does NOT work for fitting all lines.
%
% This function was written on 2020_06_25 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2020_06_25 
% -- wrote the code
% 2021_05_24
% -- revised name to fcn_geometry_fitSlopeInterceptNPoints
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

if flag_check_inputs == 1
    % Are there the right number of inputs?
    if nargin < 1 || nargin > 2
        error('Incorrect number of input arguments')
    end
    
    % Check the points input to be length greater than or equal to 2
    fcn_geometry_checkInputsToFunctions(...
        points, '2column_of_numbers',[2 3]);
    
end

Npoints = length(points(:,1));

% Does user want to show the plots?
if 2 == nargin
    fig_num = varargin{1};
    figure(fig_num);
    flag_do_debug = 1;
else
    if flag_do_debug
        fig = figure; 
        fig_num = fig.Number;
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
 
% The code below sets up the problem as a regression form. Since the
% equation of a line is:
%     y = m*x + b
% then:
%     y1 = m*x1 + b
%     y2 = m*x2 + b
%     (etc)
% or:
%     Y = m*X+b
% or:
%     Y = [X 1]*[m b]';
% which allows one to solve for [m b] via matrix multiplication via
%
%    result = ([X 1]'*[X 1])\([X 1]*Y)
%
% where result = [m b]
% 

% Fill in X and Y
X = points(:,1);
Y = points(:,2);

% Check to see if the result is going to be singular. This happens if all
% the x values are the same, e.g. the line is vertical
if all(X == X(1))  % Are all the x values the same?
    slope = inf;
    intercept = inf;
else  % The result will be an ordinary line
    
    % If X is square already, no need to do the transpose calculations
    if Npoints == 2
        result = [X ones(Npoints,1)]\Y;
    else
        A = [X ones(Npoints,1)];
        result = (A'*A)\(A'*Y);
    end
    
    slope = result(1,1);
    intercept = result(2,1);
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
    plot(points(:,1),points(:,2),'r.');
    
    % Create an x vector of 100 points connecting lowest and highest points
    % in x, and plot the line fit
    x = linspace(min(points(:,1)),max(points(:,1)),100)';
    y = x*slope + intercept;

    if isinf(slope)  % The result is a vertical line
        y = linspace(min(points(:,2)),max(points(:,2)),100)';
    end
    
    plot(x,y,'b.');
    
end

