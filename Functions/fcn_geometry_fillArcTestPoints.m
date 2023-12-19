function [test_points, true_circle_centers, true_circle_radii] = fcn_geometry_fillArcTestPoints(seed_points, M, sigma, varargin)
% fcn_geometry_fillArcTestPoints
% given N points, with N>=2, creates a set of M points per unit distance
% between these points randomly distributed with variance sigma.
%
% [test_points] = fcn_geometry_fillArcTestPoints(seed_points, M, sigma)
%
% INPUTS:
%      seed_points: a Nx2 vector where N is the number of points, but at
%      least 2.
%
%      M: the number of test points to generate per unit
%      distance.
%
%      sigma: athe standard deviation in points
%
% OUTPUTS:
%
%      test_points: a list of test points used to test regression fitting
%
%      true_circle_centers, true_circle_radii: the true values of the
%      circle equations used to fill the points
% 
% DEPENDENCIES:
%
%      none
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_fillArcTestPoints
% for a full test suite.
%
% NOTE: This function does NOT work for fitting all lines.
%
% This function was written on 2023_12_17 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2023_12_17 
% -- wrote the code using fcn_geometry_fillLineTestPoints as starter

flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 34838;
    figure(debug_fig_num);
    clf;
else
    debug_fig_num = []; 
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

if flag_check_inputs == 1
    % Are there the right number of inputs?
    narginchk(3,4);

    % Check the points input to be length greater than or equal to 2
    fcn_geometry_checkInputsToFunctions(...
        seed_points, '2column_of_numbers',[2 3]);

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
N_segments = length(seed_points(:,1)) -2;

% Find circle center and radius for each set of 3 points
[circleCenter, circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points,debug_fig_num);

true_circle_centers = circleCenter;
true_circle_radii   = circleRadius;

test_points = [];
for ith_point = 1:N_segments
    projection_to_first_point = seed_points(ith_point,:) - circleCenter(ith_point,:);
    projection_to_second_point = seed_points(ith_point+1,:) - circleCenter(ith_point,:);
    projection_to_third_point = seed_points(ith_point+2,:) - circleCenter(ith_point,:);

    unit_project_to_first_point = fcn_geometry_calcUnitVector(projection_to_first_point);
    unit_project_to_second_point = fcn_geometry_calcUnitVector(projection_to_second_point);
    unit_project_to_third_point = fcn_geometry_calcUnitVector(projection_to_third_point);

    cross_product1_to_2 = cross([unit_project_to_first_point 0], [unit_project_to_second_point 0]);
    cross_product2_to_3 = cross([unit_project_to_second_point 0], [unit_project_to_third_point 0]);

    dot_product1_to_2 = sum(unit_project_to_first_point.*unit_project_to_second_point,2);
    dot_product2_to_3 = sum(unit_project_to_second_point.*unit_project_to_third_point,2);
    
    cross_angle1_to_2_sign = asin(cross_product1_to_2(3));
    cross_angle2_to_3_sign = asin(cross_product2_to_3(3));
    dot_angle1_to_2_magnitude = acos(dot_product1_to_2);
    dot_angle2_to_3_magnitude = acos(dot_product2_to_3);


    if ith_point<N_segments
        cross_angle = sign(cross_angle1_to_2_sign)*dot_angle1_to_2_magnitude;
    else
        cross_angle = sign(cross_angle1_to_2_sign)*dot_angle1_to_2_magnitude + sign(cross_angle2_to_3_sign)*dot_angle2_to_3_magnitude;
    end

    angle_start = atan2(unit_project_to_first_point(1,2),unit_project_to_first_point(1,1));
    % angle_end = angle_start + cross_angle;
    total_distance = abs(cross_angle)*circleRadius(ith_point); % The arc distance

    % Find the number points
    projection_distances = (0:(1/M):total_distance)';
    N_points = length(projection_distances);

    % Convert to angles, and then to points
    angles = angle_start + sign(cross_angle)*projection_distances./circleRadius(ith_point);
    unit_radial_vectors = [cos(angles) sin(angles)];
    orthogonal_distances = randn(N_points,1)*sigma + ones(N_points,1)*circleRadius(ith_point);
    perturbed_points = ones(N_points,1)*circleCenter(ith_point,:) + orthogonal_distances.*unit_radial_vectors;

    if flag_do_debug
        figure(debug_fig_num);
        hold on;
        grid on;
        axis equal

        plot(circleCenter(ith_point,1), circleCenter(ith_point,2), 'r+','MarkerSize',30);
        plot(seed_points(ith_point,1), seed_points(ith_point,2), 'g.','MarkerSize',20);
        plot(seed_points(ith_point+1,1), seed_points(ith_point+1,2), 'b.','MarkerSize',20);
        plot(seed_points(ith_point+2,1), seed_points(ith_point+2,2), 'r.','MarkerSize',20);

        plot(perturbed_points(:,1), perturbed_points(:,2), 'k.','MarkerSize',20);
    end

    % Convert to test points
    test_points = [test_points; perturbed_points]; %#ok<AGROW>
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
    figure(fig_num);
    hold on;
    grid on;

    % Plot the input points
    plot(seed_points(:,1),seed_points(:,2),'r.','MarkerSize',20);
    
    % Plot the results
    plot(test_points(:,1),test_points(:,2),'b.','MarkerSize',10);

    % Make axis slightly larger
    temp = axis;
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

    
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







