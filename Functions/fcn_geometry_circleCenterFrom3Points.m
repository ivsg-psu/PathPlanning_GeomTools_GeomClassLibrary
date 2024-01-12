function [centers,radii] = fcn_geometry_circleCenterFrom3Points(points1,varargin)
% fcn_geometry_circleCenterFrom3Points calculates the center of a circle from
% three points given as vectors in x and y
%
% FORMAT:
%
% [centers,radii] = fcn_geometry_circleCenterFrom3Points(points,(fig_num))
%
%  OR
%
% [centers,radii] = fcn_geometry_circleCenterFrom3Points(points1, points2, points2,(fig_num))
%
% INPUTS:
%
%      points: a Nx2 vector where N is at least 3. If N = 3, a circle will be
%      fit between these threee points, if N = 4 or more, then one circle
%      will be fit to the first three points, another cicle to the next
%      three points, etc.
%
%      OR
%
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      centers: an [(N-2)x1] vector of the centers of the circles, in [x y]
%
%      radii: the radius of each the circles, as an [(N-2)x1] vector
%
% DEPENDENCIES:
%
%      fcn_geometry_checkInputsToFunctions
%      fcn_geometry_plotCircle
%
% EXAMPLES:
%      
%      % BASIC example
%      points = [0 0; 1 4; 0.5 -1];
%      [centers,radii] = fcn_geometry_circleCenterFrom3Points(points,1)
% 
%      % ADVANCED example that uses vectors of x and y
%      points = [0 0; 1 4; 0.5 -1; -1 4];
%      [xc,yc,radii] = fcn_geometry_circleCenterFrom3Points(points,1)
%
%      % ADVANCED example that lets user select N points 
%      figure(1); clf; grid on; axis equal;
%      points = ginput; % Get arbitrary N points until user hits return
%      [xc,yc,radii] = fcn_geometry_circleCenterFrom3Points(points,1)     
%
% See the script: script_test_fcn_geometry_circleCenterFrom3Points
% for a full test suite.
%
% This function was written on 2020_03_20 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2020_03_20 - wrote the code
% 2020_05_22 - added more comments, particularly to explain inputs more
% clearly
% 2021_05_23 
% -- merged previous function into geometry class
% -- automated input argument checking
% -- changed from x,y separate inputs into points inputs
% 2023_12_27
% -- added external environment test
% -- added speed-up wherein if fig_num set to -1, it skips plotting, input
% checking, debug modes.
% 2024_01_11
% -- added conditioning test for the matrix inversion, to avoid errors for
% colinear points

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{3},-1)) || (nargin==2 && isequal(varargin{1},-1))
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


if flag_max_speed==0
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(1,4);

        % Check the points input
        fcn_DebugTools_checkInputsToFunctions(...
            points1, '2column_of_numbers');
    end
end

% Is the user giving separated input point vectors?
flag_use_separated_point_inputs = 0;
if nargin>2
    flag_use_separated_point_inputs = 1;
    temp = varargin{1};
    if ~isempty(temp)
        points2 = temp;
    else
        error('Expected 2nd input to be a point type')
    end

    if nargin>=3
        temp = varargin{2};
        if ~isempty(temp)
            points3 = temp;
        else
            error('Expected 3rd input to be a point type')
        end
    end

    N_points = length(points1(:,1));

    if flag_check_inputs
        % Check the points2 input
        fcn_DebugTools_checkInputsToFunctions(...
            points2, '2column_of_numbers',[N_points N_points]);

        % Check the points3 input
        fcn_DebugTools_checkInputsToFunctions(...
            points3, '2column_of_numbers',[N_points N_points]);
    end
end

% Does user want to show the plots?
flag_do_plots = 0;
if 0==flag_max_speed
    if (2 == nargin || 4 == nargin)
        temp = varargin{end};
        if ~isempty(temp)
            fig_num = temp;
            flag_do_plots = 1;
        end
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


if flag_use_separated_point_inputs
    centers = zeros(N_points,2);
    radii   = zeros(N_points,1);
    for ith_row = 1:length(points1(:,1))
        if flag_max_speed ==1
            [center,radius] = fcn_geometry_circleCenterFrom3Points([points1(ith_row,:); points2(ith_row,:); points3(ith_row,:);],-1);
        else
            [center,radius] = fcn_geometry_circleCenterFrom3Points([points1(ith_row,:); points2(ith_row,:); points3(ith_row,:);]);
        end
        centers(ith_row,:) = center;
        radii(ith_row)   = radius;
    end
else
    % The code below is set up to be vectorized if there are more than one
    % solution. Since the code is quite different looking for each, they are
    % separated out. However, it may be that the N-solution case works for N is
    % equal to 1. This was not tested.

    % Do some pre-calculations
    num_solutions = length(points1(:,1))-2; % This is the number of solutions to expect
    r_squared = sum(points1.^2,2); % These are the radii-squared of points from origin
    diff_points = diff(points1);
    diff_x = diff_points(:,1);
    diff_y = diff_points(:,2);
    diff_rsquared = diff(r_squared);

    if 1 == num_solutions % Expecting just one solution. No need for big A, b matrices
        % solve for the center point
        A = [diff_x diff_y];
        b = 1/2*diff_rsquared;

    else % Simultaneous solutions to be calculated - create big A and b matrices
        % Construct the A-matrix and b matrix that will create the regressor.
        % Start by filling A and b matrices up with zeros (see notes for
        % explanation of iputs)
        A = zeros(2*num_solutions,2*num_solutions);
        b = zeros(2*num_solutions,1);

        % Fill in the non-zero portions of the matrix, which will be 1 per each
        % of the N solutions
        for i_solution = 1:num_solutions
            A(1+2*(i_solution-1):2+2*(i_solution-1),1+2*(i_solution-1)) = ...
                diff_x(i_solution:i_solution+1);
            A(1+2*(i_solution-1):2+2*(i_solution-1),2+2*(i_solution-1)) = ...
                diff_y(i_solution:i_solution+1);
            b(1+2*(i_solution-1):2+2*(i_solution-1),1) = ...
                1/2*diff_rsquared(i_solution:i_solution+1);
        end
    end

    % Check the conditioning of the matrix
    warning('on','backtrace');
    if abs(det(A))>1E-8
        % Solve the center points
        centers = A\b;
        centers = reshape(centers,2,length(centers(:,1))/2);
        centers = centers'; % Make it into a column vector
    else
        centers = [inf inf];
    end


    % NOTE: the following line is the slowest in the code. It can be sped
    % up if we do not take the square root
    radii = sum((points1(1:num_solutions,:)-centers).^2,2).^0.5;
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
    
    %Prep the figure
    figure(fig_num);
    hold on % allow multiple plot calls
    
    if flag_use_separated_point_inputs
        plot(points1(:,1),points1(:,2),'go');  % Plot all the input points
        plot(points2(:,1),points2(:,2),'bo');  % Plot all the input points
        plot(points3(:,1),points3(:,2),'ro');  % Plot all the input points
    else
        plot(points1(:,1),points1(:,2),'ro');  % Plot all the input points
    end
    plot(centers(:,1),centers(:,2),'g+'); % Plot all the circle centers

    axis equal;
    grid on; grid minor;

    % plot all the circle fits    
    for i_fit = 1:length(centers(:,1))             
        fcn_geometry_plotCircle(centers(i_fit,:),radii(i_fit,1));
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

