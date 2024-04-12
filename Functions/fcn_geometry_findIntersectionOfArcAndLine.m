function intersectionPoints = fcn_geometry_findIntersectionOfArcAndLine(lineStruct, arcStruct, tolerance, varargin)
%% fcn_geometry_findIntersectionOfArcAndLine
%
% This function finds the intersection point(s) of a line segment and an
% arc by taking the vector regression line segment, regression arc,
% tolerance, and fig_num as the inputs.
%
%
% FORMAT: 
%
% intersectionPoints = fcn_geometry_findIntersectionOfArcAndLine(lineStruct, arcStruct, tole, fig_num)
% 
% INPUTS:
%
% lineStruct: This is a cell struct array of a vector regression line
% segment. This cell struct consists of the first end point, last end
% point, fit type, and fit parameters of a vector regression line segment
%
% arcStruct: TThis is a cell struct array of a regression arc. This cell
% struct consists of the first end point, last end point, fit type, and fit
% parameters of a regression arc
%
% tolerance: This tolerance is used to determine whether there are one or
% two intersection points when a line segment intersects with an arc. If
% the distance between the center and the intersection point falls within
% the range of the arc's radius plus or minus the tolerance, only one
% intersection point is outputted.
%
%
% (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS: 
%
% intersectionPoints: Based on the intersection of a line segment with an
% arc, the intersectionPoint(s) is(are) outputted. 
%
%
% DEPENDENCIES:
%
%   No Dependencies
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_findIntersectionOfArcAndLine
% for a full test suite.
%
% This function was written on 2024_04_11 by Aneesh Batchu
% Questions or comments? abb6486@psu.edu or sbrennan@psu.edu

% Revision History
% 2024_04_11 - Aneesh Batchu

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
if (0==flag_max_speed) && (3<= nargin)
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



% Location of S (Start Point of Line Segment) 
pointS = lineStruct.firstEndPoint;

% Location of E (End Point of Line Segement)
pointE = lineStruct.lastEndPoint;

% Location of C (Arc center)
pointC = arcStruct.fitParameters(1,1:2);

% Calculate vector EC
vectorEC = pointC - pointE;

% Calculate vector ES
vectorES = pointS - pointE;

% Calculate the unit vector of vectorES
vectorES_magnitude = sum(vectorES.^2,2).^0.5;
unit_vectorES = vectorES./vectorES_magnitude;

% Calculate unit orthogonal vector of vectorES
unit_orthogonal_vectorES = unit_vectorES*[0 1; -1 0];

% Find the distance between the arc center and start point of the line
% segment by calculating the dot product of VectorEC and unit_orthogonal_vectorES
dist_btw_pointC_and_pointS = dot(vectorEC, unit_orthogonal_vectorES);

% Radius of the regression Arc
radiusArc = arcStruct.fitParameters(1,3);

if dist_btw_pointC_and_pointS > radiusArc
    disp('No Intersection Between the Line and Arc')
    intersectionPoints = [nan, nan]; 
end

% tolerance = 0.1;
% One intersection point, if the distance between arc center and start
% point of segment is equal to radius +/- tolerance
if (abs(dist_btw_pointC_and_pointS) <= (radiusArc + tolerance)) && (abs(dist_btw_pointC_and_pointS) >= (radiusArc - tolerance))

    intersectionPoints = pointE + dot(vectorEC, unit_vectorES)*unit_vectorES;
    
end


% Two intersection points, if the distance between arc center and start
% point of segment is less than radius of the Arc - tolerance
if abs(dist_btw_pointC_and_pointS) < (radiusArc - tolerance)

    % Center point of the two intersection points
    % centerPoint_of_intersectionPoints = vectorES + dot(vectorEC, unit_vectorES)*unit_vectorES;
    centerPoint_of_intersectionPoints = pointE + dot(vectorEC, unit_vectorES)*unit_vectorES;
    plot(centerPoint_of_intersectionPoints(1,1), centerPoint_of_intersectionPoints(1,2), '.', 'Color', 'g', 'MarkerSize',30);

    % Opposite side of right angle triangle: dist_btw_pointC_and_pointS
    dist_btw_centerPoint_and_pointC = dot(vectorEC, unit_orthogonal_vectorES);

    % Adjacent side of the right angle triangle
    dist_btw_centerPoint_and_intersectionPoints = abs((radiusArc^2 - dist_btw_centerPoint_and_pointC^2))^0.5;
    
    %
    intersectionPoint1 = centerPoint_of_intersectionPoints + dist_btw_centerPoint_and_intersectionPoints*unit_vectorES; 

    intersectionPoint2 = centerPoint_of_intersectionPoints - dist_btw_centerPoint_and_intersectionPoints*unit_vectorES; 

    intersectionPoints = [intersectionPoint1; intersectionPoint2]; 

 
    
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

    figure(fig_num)

    % Plot Arc Center
    plot(pointC(1,1), pointC(1,2), 'c.','MarkerSize',30);

    % Plot the vectorEC in GREEN
    quiver(pointE(1,1), pointE(1,2), vectorEC(1,1), vectorEC(1,2), 'Color', 'g', 'LineWidth', 4);

    % Plot the vectorES in BLUE
    quiver(pointE(1,1), pointE(1,2), vectorES(1,1), vectorES(1,2), 'Color', 'b', 'LineWidth', 4);

    % Plot the unit orthogonal vector ES (unit_vectorES) in YELLOW
    quiver(pointE(1,1), pointE(1,2), unit_vectorES(1,1), unit_vectorES(1,2), 'Color', 'y', 'LineWidth', 3);

    % Plot the unit orthogonal vector ES (unit_orthogonal_vectorES) in RED
    quiver(pointE(1,1), pointE(1,2), unit_orthogonal_vectorES(1,1), unit_orthogonal_vectorES(1,2), 'Color', 'r', 'LineWidth', 3);

    if ~isnan(intersectionPoints)
        % Plot the intersectionPoint(s) in RED
        if size(intersectionPoints,1) == 1
            plot(intersectionPoints(1,1), intersectionPoints(1,2), '.', 'Color', 'r', 'MarkerSize',40);
        elseif size(intersectionPoints, 1) == 2
            plot(intersectionPoints(:,1), intersectionPoints(:,2), '.', 'Color', 'r', 'MarkerSize',40)
        end
    end
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
