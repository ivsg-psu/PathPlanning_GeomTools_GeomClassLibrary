function [dist] = fcn_geometry_euclideanPointToPointsDistance(pt1,pt2, varargin)
% fcn_geometry_euclideanPointToPointsDistance
% Finds distance(s) from one point to another point(s)
%
% Format: 
% [dist] = fcn_geometry_euclideanPointToPointsDistance(pt1,pt2)
%
% [DIST] = fcn_geometry_euclideanPointToPointsDistance(PT1,PT2)
% returns:
% DIST: an n-by-1 vector of distances [d1; d2; ... ; dn], where n is the
% number of point sets
%
% with inputs:
% PT1: an n-by-2or3 series of xy or xyz points [x1 y1 z1; x2 y2 z2; ... ; xn yn zn]
% PT2: an n-by-2or3 series of xy or xyz points [x1 y1 z1; x2 y2 z2; ... ; xn yn zn]
%
% Example:
% pt1 = [1 1 5; 5 3 64; 7 2 -2];
% pt2 = [0 -3 -6; 34 1 17; 18 7 0];
% dist=fcn_general_calculation_euclidean_point_to_point_distance(pt1,pt2);
%
% INPUTS:
%      pt1,pt2: [NxM set of points]
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      dist: distances between the point(s)
%
% DEPENDENCIES:
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%      
% See the script: fcn_geometry_euclideanPointToPointsDistance
% for a full test suite.
%
% This function was written on 2018_11_17 by Seth Tau
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2018_11_17 - S. Tau
% -- original write of the code
% 2024_01_18 - S. Brennan
% -- Vectorized the code

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
        narginchk(2,3);

        % Check the pt1 input to have 2 or 3 columns
        fcn_DebugTools_checkInputsToFunctions(pt1, '2or3column_of_numbers',[2 3]);

        N_points = length(pt1(:,1));

        % Check the pt2 input to have 2 or 3 columns        
        fcn_DebugTools_checkInputsToFunctions(pt2, '2or3column_of_numbers', N_points);
    end
end

% Does user want to specify fig_num?
flag_do_plots = 0;
if (0==flag_max_speed) && (3<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp; %#ok<NASGU>
        figure(fig_num);
        flag_do_plots = 1;
    end
end

if flag_do_debug
    fig_debug = 8484; %#ok<NASGU>
else
    fig_debug = []; %#ok<NASGU>
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
dist = sum((pt1-pt2).^2,2).^0.5;

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
   % Set up the figure
    figure(fig_num);
    clf
    hold on;
    grid on; grid minor;
        
    midpoints = (pt1+pt2)/2;
    for ith_point=1:N_points
        % 2D plot?
        if length(midpoints(1,:))==2
            % Plot the points
            xdata = [pt1(ith_point,1) pt2(ith_point,1)];
            ydata = [pt1(ith_point,2) pt2(ith_point,2)];
            plot(xdata,ydata,'.-','Linewidth',3,'Markersize',20);
            
            % Label the midpoints
            text(midpoints(ith_point,1),midpoints(ith_point,2),sprintf('d - %.1f',dist(ith_point,1)));
        else
            % Plot the points
            xdata = [pt1(ith_point,1) pt2(ith_point,1)];
            ydata = [pt1(ith_point,2) pt2(ith_point,2)];
            zdata = [pt1(ith_point,3) pt2(ith_point,3)];
            plot3(xdata,ydata,zdata,'.-','Linewidth',3,'Markersize',20);

            % Label the midpoints
            text(midpoints(ith_point,1),midpoints(ith_point,2),midpoints(ith_point,3),sprintf('d - %.1f',dist(ith_point,1)));
            
            % Set to 3D view
            view(3);
        end
        
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



