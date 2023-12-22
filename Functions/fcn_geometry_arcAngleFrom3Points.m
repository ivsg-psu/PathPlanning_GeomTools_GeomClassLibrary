function [arc_angle_in_radians_1_to_2, arc_angle_in_radians_1_to_3, circle_centers, radii, start_angles_in_radians]  = fcn_geometry_arcAngleFrom3Points(points1, points2, points3, varargin)
% fcn_geometry_arcAngleFrom3Points calculates the angle between three
% points such that the angle starts at point1, passes through point2, and
% ends at point3. The angles are presented as positive.
%
% FORMAT:
%
% [arc_angle_in_radians_1_to_2, arc_angle_in_radians_1_to_3, circle_centers, radii, start_angles_in_radians]  = fcn_geometry_arcAngleFrom3Points(points1, points2, points3,(fig_num))
%
% INPUTS:
%
%      point1, point2, point3: a Nx2 vectors of point pairings. Note: this
%      function is vectorized so that multiple points can be entered
%      simultaneously
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      arc_angle_in_radians_1_to_2: an [Nx1] vector containing the angle
%      subtending from point1 to point2, passing toward point3
%
%      arc_angle_in_radians_1_to_3: an [Nx1] vector containing the angle
%      subtending from point1 to point3, passing through point2
%
%      circle_centers: the centers of the circles of the arcs, as [Nx2]
%      vector
%
%      radii: the radii of the arcs as an [Nx1] vector
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      cross
%      sign
%
% EXAMPLES:   
%
% See the script: script_test_fcn_geometry_arcAngleFrom3Points
% for a full test suite.
%
% This function was written on 2023_12_19 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2023_12_19 - sbrennan@psu.edu
% -- original write of the code

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
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

if flag_check_inputs    
    % Are there the right number of inputs?
    narginchk(3,4);
    
    % Check the points1 input
    fcn_DebugTools_checkInputsToFunctions(...
        points1, '2column_of_numbers');

    N_points = length(points1(:,1));

    % Check the points2 input
    fcn_DebugTools_checkInputsToFunctions(...
        points2, '2column_of_numbers',[N_points N_points]);

    % Check the points3 input
    fcn_DebugTools_checkInputsToFunctions(...
        points3, '2column_of_numbers',[N_points N_points]);
else
    N_points = length(points1(:,1));
end

% Does user want to show the plots?
flag_do_plots = 0;
if 4 == nargin
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
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

% Find circle center and radius for each set of 3 points
[circle_centers, radii] = fcn_geometry_circleCenterFrom3Points(points1, points2, points3);

% Find if the arcs are counterclockwise or clockwise
is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(points1, points2, points3);

% Find the projection vectors from the centers of the circles to each point
% set
arc_angle_in_radians_1_to_2 = fcn_geometry_findAngleUsing2PointsOnCircle(...
    circle_centers,...
    radii,...
    points1,...
    points2,...
    is_counterClockwise);

arc_angle_in_radians_1_to_3 = fcn_geometry_findAngleUsing2PointsOnCircle(...
    circle_centers,...
    radii,...
    points1,...
    points3,...
    is_counterClockwise);

% Calculate the start angles
% Find the projection vectors from the centers of the circles to the start
projection_to_points1  = points1 - circle_centers;

% Calculate the angle
start_angles_in_radians = atan2(projection_to_points1(:,2),projection_to_points1(:,1));




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

    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end        
    axis equal;
    hold on;
    grid on;

    % Plot the input points
    plot(points1(:,1),points1(:,2),'g.','MarkerSize',20);
    plot(points2(:,1),points2(:,2),'b.','MarkerSize',20);
    plot(points3(:,1),points3(:,2),'r.','MarkerSize',20);
    
    % Plot the circles
    fcn_geometry_plotCircle(circle_centers,radii);

    % Plot the start and end points
    plot(points1(:,1),points1(:,2),'go');
    text(points1(:,1),points1(:,2),'Start');

    plot(points2(:,1),points2(:,2),'b.','MarkerSize',10);
    text(points2(:,1),points2(:,2),'Mid');

    plot(points3(:,1),points3(:,2),'rx');
    text(points3(:,1),points3(:,2),'End');

    for i=1:length(circle_centers(:,1))
        % Plot unit vectors
        plot([circle_centers(i,1); points1(i,1)],...
            [circle_centers(i,2); points1(i,2)],'g','LineWidth',3);

        plot([circle_centers(i,1); points3(i,1)],...
            [circle_centers(i,2); points3(i,2)],'r','LineWidth',3);

        % Plot the angle value
        location = circle_centers(i,:);
        if is_counterClockwise(i,1)>0
            text(location(1,1),location(1,2),sprintf(':  %.1f deg counterclockwise',arc_angle_in_radians_1_to_3(i,1)*180/pi));
        else
            text(location(1,1),location(1,2),sprintf(':  %.1f deg clockwise',arc_angle_in_radians_1_to_3(i,1)*180/pi));
        end
    end

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
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





