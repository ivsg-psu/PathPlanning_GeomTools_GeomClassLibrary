function is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(points1, points2, points3,varargin)
% fcn_geometry_arcDirectionFrom3Points calculates whether or not 3 points
% make an arc in the clockwise or counter-clockwise diretion. The direction
% is calculated by performing a cross product between the vectors from
% points1 to points3, versus the vectors from points1 to points2. The
% funtion returns the sign of the result.
%
% FORMAT:
%
% is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(point1, point2, point3,(fig_num))
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
%      is_counterClockwise: an [Nx1] vector containing 1 if the
%      corresponding row of points creates a counter-clockwise arc
%      (positive) or -1 if it would create a clockwise arc (e.g. negative).
%      It returns 0 if the direction is undefined (if points are repeated,
%      colinear, etc.).
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      cross
%      sign
%
% EXAMPLES:   
%
% See the script: script_test_fcn_geometry_arcDirectionFrom3Points
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

    % Check the points2 input
    fcn_DebugTools_checkInputsToFunctions(...
        points2, '2column_of_numbers',[N_points N_points]);
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
 
% The code below is set up to be vectorized. The method is to create a
% vector from points1 to points2 and to points3. If the 
projection_points1_to_points2 = points2 - points1;
projection_points1_to_points3 = points3 - points1;

cross_product = cross([projection_points1_to_points2 zeros(N_points,1)], [projection_points1_to_points3 zeros(N_points,1)]);

is_counterClockwise = sign(cross_product(:,3));

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
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end        

    hold on % allow multiple plot calls
    grid on;
    grid minor;
    axis equal;
    
    % Plot the inputs
    plot(points1(:,1),points1(:,2),'g.','MarkerSize',30);  % Plot all the points1
    plot(points2(:,1),points2(:,2),'b.','MarkerSize',20);  % Plot all the points2
    plot(points3(:,1),points3(:,2),'r.','MarkerSize',10);  % Plot all the points3

    axis equal;
    grid on; grid minor;

    % Calculate the location for text
    half_point = (points1+points3)/2;
    half_toward_point2 = (half_point+points2)/2;
    text_locations = half_toward_point2; % +half_ortho_distance.*unit_orthogonal_points1_to_points3;
    difference = sum((half_point-points2).^2,2).^0.5;

    % plot all the results    
    for i_fit = 1:N_points    
        % Plot the connecting lines
        plot([points1(i_fit,1) points3(i_fit,1)],[points1(i_fit,2) points3(i_fit,2)],'b-','LineWidth',3);

        % Plot the midpoint lines
        plot([half_point(i_fit,1) points2(i_fit,1)],[half_point(i_fit,2) points2(i_fit,2)],'b-','LineWidth',3);

        % Label the result

        nudge = 0.1*difference(i_fit);
        text(text_locations(i_fit,1)+nudge,text_locations(i_fit,2)+nudge,sprintf('%.0d',is_counterClockwise(i_fit)));
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

