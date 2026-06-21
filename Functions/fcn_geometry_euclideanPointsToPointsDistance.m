function [dist] = ...
    fcn_geometry_euclideanPointsToPointsDistance(...
    points1,...
    points2,...
    varargin)
%% fcn_geometry_euclideanPointsToPointsDistance 
% 
% This calculates the distance(s) between a [Nxd] vector of points,
% POINTS1, and another [Nxd] vector of points, POINTS2, where d is the
% dimension.
%
% FORMAT:
%
% [DIST] = fcn_geometry_euclideanPointsToPointsDistance(POINTS1,POINTS2,(figNum))
%
% INPUTS:
%
%      POINTS1: an Nx2 or Nx3 series of xy or xyz points 
%      in the form: [x1 y1 z1; x2 y2 z2; ... ; xn yn zn]
%
%      POINTS2: an Nx2 or Nx3 series of xy or xyz points 
%      in the form: [x1 y1 z1; x2 y2 z2; ... ; xn yn zn]
%
%      (OPTIONAL INPUTS)
% 
%      figNum: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      DIST: an N x  1 vector of distances [d1; d2; ... ; dn], where N is
%      the number of point sets
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%         pt1 = [1 1 5; 5 3 64; 7 2 -2];
%         pt2 = [0 -3 -6; 34 1 17; 18 7 0];
%         dist=fcn_geometry_euclideanPointsToPointsDistance(pt1,pt2);
%
% See the script: script_test_fcn_geometry_euclideanPointsToPointsDistance
% for a full test suite.
%
% This function was written on 2018_11_17 by Seth Tau
% Questions or comments? sat5340@psu.edu 

% REVISION HISTORY:
%
% 2021_05_28 by Sean Brennan, sbrennan@psu.edu
% - In fcn_geometry_euclideanPointsToPointsDistance
%   % * revised function to prep for geometry class
%   % * rewrote function to use vector sum
%   % * added plotting option
% 
% 2021_05_28 by Sean Brennan, sbrennan@psu.edu
% - In fcn_geometry_euclideanPointsToPointsDistance
%   % * fixed comments, added debugging option
% 
% 2024_01_17 by Aneesh Batchu, abb6486@psu.edu
%  - In fcn_geometry_euclideanPointsToPointsDistance
%   % * added max speed options 
% 
% 2024_01_17 by Aneesh Batchu, abb6486@psu.edu
%  - In fcn_geometry_euclideanPointsToPointsDistance
%   % * Updated the function headers to "standard" form 

%% Debugging and Input checks
% flag_check_inputs = 1; % Set equal to 1 to check the input arguments
% flag_do_plot = 0;      % Set equal to 1 for plotting
% flag_do_debug = 0;     % Set equal to 1 for debugging

% Check if flag_max_speed set. This occurs if the figNum variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 3; % The largest Number of argument inputs to the function
flag_max_speed = 0;
if (nargin==MAX_NARGIN && isequal(varargin{end},-1))
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

%% check input arguments?
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
        narginchk(2,MAX_NARGIN);

        % Check the points1 input
        fcn_DebugTools_checkInputsToFunctions(...
            points1, '2or3column_of_numbers');

        % Use number of rows in points1 to calculate Npoints
        Npoints = length(points1(:,1));

        % Check the points2 input, forcing length to match points1
        fcn_DebugTools_checkInputsToFunctions(...
            points2, '2or3column_of_numbers',Npoints);
    end
end
    

% % Does user want to show the plots?
% flag_do_plot = 0;
% if (0==flag_max_speed) && (3 == nargin) 
%     temp = varargin{1};
%     if ~isempty(temp)
%         figNum = temp;
%         figure(figNum);
%         flag_do_plot = 1;
%     end
% else
%     if flag_do_debug
%         fig = figure; 
%         figNum = fig.Number;
%         flag_do_plot = 1;
%     end
% end

% Does user want to show the plots?
flag_do_plots = 0; % Default is to show plots
figNum = []; % Empty by default
if (0==flag_max_speed) && (MAX_NARGIN == nargin)
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        figNum = temp;
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

dist = sum((points1-points2).^2,2).^0.5;

%% Plot results?
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
    % Set up the figure
    figure(figNum);
    clf
    hold on;
    grid on; grid minor;
        
    midpoints = (points1+points2)/2;
    for ith_point=1:Npoints
        % 2D plot?
        if length(midpoints(1,:))==2
            % Plot the points
            xdata = [points1(ith_point,1) points2(ith_point,1)];
            ydata = [points1(ith_point,2) points2(ith_point,2)];
            plot(xdata,ydata,'.-','Linewidth',3,'Markersize',20);
            
            % Label the midpoints
            text(midpoints(ith_point,1),midpoints(ith_point,2),sprintf('d - %.1f',dist(ith_point,1)));
        else
            % Plot the points
            xdata = [points1(ith_point,1) points2(ith_point,1)];
            ydata = [points1(ith_point,2) points2(ith_point,2)];
            zdata = [points1(ith_point,3) points2(ith_point,3)];
            plot3(xdata,ydata,zdata,'.-','Linewidth',3,'Markersize',20);

            % Label the midpoints
            text(midpoints(ith_point,1),midpoints(ith_point,2),midpoints(ith_point,3),sprintf('d - %.1f',dist(ith_point,1)));
            
            % Set to 3D view
            view(3);
        end
        
    end
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end


end % Ends the function



