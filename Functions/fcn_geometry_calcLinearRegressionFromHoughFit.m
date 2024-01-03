function [regression_fit_line_segment, domain_box] = fcn_geometry_calcLinearRegressionFromHoughFit(source_points, associated_points_in_domain, varargin)
% fcn_geometry_calcLinearRegressionFromHoughFit
% Given a set of points that are matched via a Hough vote, finds the linear
% regression fit line and domain box. 
% 
% Format: 
% [regression_fit_line_segment, domain_box] = fcn_geometry_calcLinearRegressionFromHoughFit(source_points,associated_points_in_domain, (fig_num))
%
% INPUTS:
%      source_points: a 2x2 matrix in the format: 
% 
%      [start_point_x  start_point_y; end_point_x end_point_y]
%
%      of the points used to "aim" the Hough voting vector. These points
%      are used to determine the direction of the resulting segment fit.
%
%      associated_points_in_domain: an Nx2 list of points that should be
%      fit with regression, identified as within the domain according to
%      Hough voting.
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot the results.
%
% OUTPUTS:
%
%      regression_fit_line_segment: a 2x2 matrix of the format:
%
%       [start_point_x  start_point_y; end_point_x end_point_y]
%
%      these points are the start and end of the segment that is on the
%      best-fit line fit of the associated_points_in_domain. The start
%      point is aligned to be perfectly orthogonal to the "lowest" point in
%      the associated points in domain
%
%      domain_box: the box that encloses the 2-standard-deviation interval
%      around the line segment.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_fitSlopeInterceptNPoints
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_calcLinearRegressionFromHoughFit
% for a full test suite.
%
% This function was written on 2023_12_15 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2023_12_14 
% -- wrote the code


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

if flag_check_inputs == 1
    % Are there the right number of inputs?
    narginchk(2,3);
    
    % Check the source_points input to be length exactly equal to 2
    fcn_DebugTools_checkInputsToFunctions(...
        source_points, '2column_of_numbers',[2 2]);
    
    % Check the associated_points_in_domain input to be length greater than or equal to 2
    fcn_DebugTools_checkInputsToFunctions(...
        associated_points_in_domain, '2column_of_numbers',[2 3]);    
end


% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if 3<= nargin
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
 % Find the root point
base_point_of_domain = source_points(1,:);
end_point_of_domain  = source_points(2,:);
unit_tangent_vector_of_domain = fcn_geometry_calcUnitVector(end_point_of_domain - base_point_of_domain);
% unit_orthogonal_vector_of_domain = unit_tangent_vector_of_domain*[0 1; -1 0];


% Find domain ranges related to point-to-point fit
projection_vectors_from_base = associated_points_in_domain - base_point_of_domain;
tangent_distances = sum(projection_vectors_from_base.*unit_tangent_vector_of_domain,2);
% orthogonal_distances = sum(projection_vectors_from_base.*unit_orthogonal_vector_of_domain,2);

% Sort the points in direction of point-to-point projection
[~,sorted_indicies] = sort(tangent_distances);
sorted_points_in_domain = associated_points_in_domain(sorted_indicies,:);

% Find regression-fit line
[slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(sorted_points_in_domain);

% Convert slope into projection and orthogonal vectors
projection_vector_angle = atan2(unit_tangent_vector_of_domain(:,2),unit_tangent_vector_of_domain(:,1));
regression_vector_angle = atan2(slope,1);

if abs(projection_vector_angle - regression_vector_angle)>pi/2
    if flag_do_debug
        fprintf(1,'Correcting angle...\n')
    end
    regression_vector_angle = regression_vector_angle + pi;
end

if 1==0
    if abs(projection_vector_angle - regression_vector_angle)>(10*pi/180)
        warning('Projection angle is (degrees): %.3f',projection_vector_angle*180/pi);
        warning('Regression angle is (degrees): %.3f',regression_vector_angle*180/pi);
        error('Regression angle and fitted angle are more than 10 degrees different. An error likely occurred!');
    end
end

unit_converted_projection_vector = [cos(regression_vector_angle) sin(regression_vector_angle)];
unit_converted_orthogonal_vector = unit_converted_projection_vector*[0 1; -1 0];


% Find new base point that uses regression fit
y_value_on_line = slope*base_point_of_domain(1,1)+intercept;
point_on_line = [base_point_of_domain(1,1) y_value_on_line];


% Project the base point onto the line segment using the dot
% product to create a new base point
projection_to_lowest_point = sorted_points_in_domain(1,:) - point_on_line;
transverse_distance_to_lowest_point = sum(projection_to_lowest_point.*unit_converted_projection_vector,2);
min_regression_point = point_on_line + transverse_distance_to_lowest_point*unit_converted_projection_vector;

projection_to_highest_point = sorted_points_in_domain(end,:) - point_on_line;
transverse_distance_to_highest_point = sum(projection_to_highest_point.*unit_converted_projection_vector,2);
max_regression_point = point_on_line + transverse_distance_to_highest_point*unit_converted_projection_vector;

% Find the transverse distances and standard deviations of
% transverse distances
projections = sorted_points_in_domain - point_on_line;
orthogonal_distances = sum(projections.*unit_converted_orthogonal_vector,2);
std_dev_transverse_distance = std(orthogonal_distances);

% Save point-to-point best-fit line segment
regression_fit_line_segment = [min_regression_point; max_regression_point];

% Calculate the domain box

% min_tangent_distance = min(tangent_distances); % Can speed this up by using indicies above
% max_tangent_distance = max(tangent_distances); % Can speed this up by using indicies above

max_orthogonal_distance = 5*std_dev_transverse_distance; % max(abs(orthogonal_distances));

domain_box = ...
    [...
    min_regression_point - max_orthogonal_distance*unit_converted_orthogonal_vector;
    max_regression_point - max_orthogonal_distance*unit_converted_orthogonal_vector;
    max_regression_point + max_orthogonal_distance*unit_converted_orthogonal_vector;
    min_regression_point + max_orthogonal_distance*unit_converted_orthogonal_vector;
    min_regression_point - max_orthogonal_distance*unit_converted_orthogonal_vector;
    ];


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
    % Plot the results of the point fit
    figure(fig_num);
    hold on;
      
    % Plot the input points
    plot(source_points(1,1),source_points(1,2),'g.','MarkerSize',30);
    plot(source_points(2,1),source_points(2,2),'r.','MarkerSize',30);
    h_plot = plot(associated_points_in_domain(:,1),associated_points_in_domain(:,2),'.','MarkerSize',10);
    current_color = get(h_plot,'Color');



    % Plot the line fit
    plot(regression_fit_line_segment(:,1),regression_fit_line_segment(:,2),'-','LineWidth',3,'Color',current_color);

    % Plot the domain
    domainShape = polyshape(domain_box(:,1),domain_box(:,2),'Simplify',false,'KeepCollinearPoints',true);
    plot(domainShape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',1,'EdgeAlpha',0);

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


