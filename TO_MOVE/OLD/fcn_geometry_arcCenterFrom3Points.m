function [arcCenter, arcRadius, theta_start, theta_end] = fcn_geometry_arcCenterFrom3Points(point1, point2, point3, varargin)
% fcn_geometry_arcCenterFrom3Points calculates the center of an arc from
% three points given as vectors in x and y
%
% FORMAT:
%
% [arcCenter, arcRadius, theta_start, theta_end] = fcn_geometry_arcCenterFrom3Points(points1, points2, points2,(fig_num))
%
% INPUTS
%
% threePoints - point1 = [x1,y1], point2 = [x2,y2], point3 = [x3,y3]  
%
% (OPTIONAL INPUTS)
%
% fig_num: a figure number to plot the results
%
% OUTPUTs
%
% arcCenter - center of the fitted arc
% arcRadius - radius of the fitted arc
%
% DEPENDENCIES: 
% 
%       (None)
%
% This function was written on 2020_03_24 by A. Batchu
% Questions or comments? abb6486@psu.edu 

% Revision history:
% 2023_12_24 - wrote the code

flag_do_plots = 0; % Flag to plot
flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 34838;
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

if flag_check_inputs
    % Are there the right number of inputs?
    narginchk(3,4);

end

% Does user want to show the plots?
if 4 == nargin
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        figure(fig_num);
        flag_do_plots = 1; 
    end
else
    if flag_do_debug
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

% Input must contain three points
threePoints = [point1; point2; point3];

% Calculating the center of an arc
%
% 1- Find any two chords of the circle - Calculate the distance between two
% point pairs.
% 2- Find the mid-point of the two-point pairs.
% 3- Find the perpendicular bisector of those two chords.
% 4-The intersection of those two perpendicular bisectors is the center of
% the circle.
%
% These formulas are derived based on the geometric properties mentioned
% above.

Denominator = 2 * ( point1(1,1)*(point2(1,2)-point3(1,2)) + point2(1,1)*(point3(1,2)-point1(1,2)) + point3(1,1)*(point1(1,2)-point2(1,2)) );

xcoor_center = ( (point1(1,1)^2+point1(1,2)^2)*(point2(1,2)-point3(1,2)) + ...
    (point2(1,1)^2+point2(1,2)^2)*(point3(1,2)-point1(1,2)) + ...
    (point3(1,1)^2+point3(1,2)^2)*(point1(1,2)-point2(1,2)) ) / Denominator;

ycoor_center = ( (point1(1,1)^2+point1(1,2)^2)*(point3(1,1)-point2(1,1)) + ...
    (point2(1,1)^2+point2(1,2)^2)*(point1(1,1)-point3(1,1)) + ...
    (point3(1,1)^2+point3(1,2)^2)*(point2(1,1)-point1(1,1)) ) / Denominator;

arcCenter = [xcoor_center, ycoor_center];

% Circle radius is obtained by calculating the distance between the center
% and any point on the circle.

arcRadius = ( (point1(1,1)-xcoor_center)^2 + (point1(1,2)-ycoor_center)^2 )^0.5;

% Arc start and end angle
theta_start = atan2(point1(1,2) - arcCenter(1,2), point1(1,1) - arcCenter(1,1));
theta_end = atan2(point3(1,2) - arcCenter(1,2), point3(1,1) - arcCenter(1,1));

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

    % Plot the points and the fitted circle
    theta = linspace(theta_start, theta_end, 100);

    % theta = linspace(0, 2*pi, 100);
    circle_xcoor = arcCenter(1,1) + arcRadius*cos(theta);
    circle_ycoor = arcCenter(1,2) + arcRadius*sin(theta);

    % plot the circle
    figure(fig_num)
    plot(circle_xcoor, circle_ycoor, 'r', 'LineWidth',2);
    hold on

    % Plot the threePoints on the circle
    scatter(threePoints(:,1), threePoints(:,2), 'b','filled');

    % Plot the center
    scatter(arcCenter(1,1), arcCenter(1,2), 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');

    % Details
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Fitted Circle Passing Through Three Points');
    legend('Fitted Circle', 'Points', 'Center');

    axis equal
    box on
    grid on

    % texts = 0, for no text on figure
    texts = 1;

    if texts
        % Three points
        text(point1(1,1), point1(1,2), '  Point 1');
        text(point2(1,1), point2(1,2), '  Point 2');
        text(point3(1,1), point3(1,2), '  Point 3');
        % Center
        text(arcCenter(1,1), arcCenter(1,2), ['  Center (' num2str(arcCenter(1,1)) ', ' num2str(arcCenter(1,2)) ')' ])
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
