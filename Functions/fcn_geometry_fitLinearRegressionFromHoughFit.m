function [regression_domain, std_dev_orthogonal_distance] = fcn_geometry_fitLinearRegressionFromHoughFit(Hough_domain, varargin)
% fcn_geometry_fitLinearRegressionFromHoughFit Given a domain containting a
% set of points that are matched via a Hough vote, finds the regression fit
% vector and domain box. 
% 
% NOTE: the vector fit is not the same as a least-squares linear
% regression, which minimizes sum-of-squares of the VERTICAL errors between
% a line fit and the respective points. The method here is similar to
% total-least-squares. Namely, the vector is found that approximately
% minimizes the sum-of-squares distance between the vector and the
% orthogonal projection to each point. Thus, the fit is minimizing
% ORTHOGONAL errors. Unlike linear regression, this method works for data
% aligned vertically.
% 
% Format: 
% [regression_fit_line_segment, domain_box] = fcn_geometry_fitLinearRegressionFromHoughFit(Hough_domain, (fig_num))
%
% INPUTS:
%      Hough_domain: a structure that records details of the domain of
%      fitting. See fcn_geometry_fillEmptyDomainStructure for details.
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      regression_domain: a structure that records details of the domain of
%      fitting. See fcn_geometry_fillEmptyDomainStructure for details.
%
%      std_dev_orthogonal_distance: the standard deviation in the point
%      fit, as measured in the transverse direction (orthogonal to the line
%      fit). E.g., this is the total-least-squares standard deviation.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_calcUnitVector
%      fcn_geometry_fitVectorToNPoints
%      fcn_geometry_fitSlopeInterceptNPoints 
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_fitLinearRegressionFromHoughFit
% for a full test suite.
%
% This function was written on 2023_12_15 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2023_12_14 - S Brennan
% -- wrote the code
% 2024_01_03 - S. Brennan
% -- added fast mode option
% -- added environmental variable options
% 2024_01_12 - S. Brennan
% -- fixed output angles to 0 to 2*pi range
% 2024_01_13 - S. Brennan
% -- replaced linear least-squares-y regression with vector regression
% -- added std_dev_transverse_distance as an output
% -- fixed bug doe to wrap-around errors when doing angle alignment check
% 2024_01_15 - S. Brennan
% -- switched inputs and outputs to domain types

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==2 && isequal(varargin{end},-1))
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
        narginchk(1,2);

        % % Check the source_points input to be length exactly equal to 2
        % fcn_DebugTools_checkInputsToFunctions(...
        %     source_points, '2column_of_numbers',[2 2]);
        % 
        % % Check the associated_points_in_domain input to be length greater than or equal to 2
        % fcn_DebugTools_checkInputsToFunctions(...
        %     associated_points_in_domain, '2column_of_numbers',[2 3]);
    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (2<= nargin)
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

% Initialze the domain structure for output
regression_domain = fcn_geometry_fillEmptyDomainStructure;

% Pull out key variables
best_fit_type            = Hough_domain.best_fit_type;
points_in_domain         = Hough_domain.points_in_domain;
best_fit_parameters      = Hough_domain.best_fit_parameters;

% Check to make sure we have a Hough fit input
switch best_fit_type
    case {'Hough line'}
        regression_domain.best_fit_type = 'Vector regression line fit';
    case {'Hough segment'}
        regression_domain.best_fit_type = 'Vector regression segment fit';
    otherwise
        error('A domain was tested for fitting that was not a Hough fit');
end
regression_domain.points_in_domain = points_in_domain;


% Pull out parameters from the Hough domain fit
associated_points_in_domain = points_in_domain;

%% First, sort all the points in the original direction
% Find the unit vector (in fast mode) connecting base to end point.
unit_tangent_vector_of_domain = best_fit_parameters(1,1:2);
base_point_of_domain          = best_fit_parameters(1,3:4);

% Find domain ranges related to point-to-point fit
projection_vectors_from_base = associated_points_in_domain - base_point_of_domain;
tangent_distances = sum(projection_vectors_from_base.*unit_tangent_vector_of_domain,2);
% orthogonal_distances = sum(projection_vectors_from_base.*unit_orthogonal_vector_of_domain,2);

% Sort the points in direction of point-to-point projection
[~,sorted_indicies] = sort(tangent_distances);
sorted_points_in_domain = associated_points_in_domain(sorted_indicies,:);

% What angle is the projection vector at?
projection_vector_angle = atan2(unit_tangent_vector_of_domain(:,2),unit_tangent_vector_of_domain(:,1));

%% Next, find best-fit line
% In the early form of this function, linear regression was used. The code
% for this is kept below in case we need to use it again sometime later.

flag_use_vector_regression = 1;
if flag_use_vector_regression
    [root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(sorted_points_in_domain,-1);
    base_point_on_line = root_point;
    regression_vector_angle = atan2(unit_vector(1,2),unit_vector(1,1));


else
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(sorted_points_in_domain,-1);

    % Convert slope into projection and orthogonal vectors
    regression_vector_angle = atan2(slope,1);

    % Find new base point that uses regression fit
    y_value_on_line = slope*base_point_of_domain(1,1)+intercept;
    base_point_on_line = [base_point_of_domain(1,1) y_value_on_line];

end

% Fix both angles to 0 to 2pi range
projection_vector_angle = mod(projection_vector_angle,2*pi);
regression_vector_angle = mod(regression_vector_angle,2*pi);


% Check if regression slope angle was calculated in wrong direction. If so,
% add pi onto it and fix
if abs(projection_vector_angle - regression_vector_angle)>pi/2
    if flag_do_debug
        fprintf(1,'Correcting angle...\n')
    end
    regression_vector_angle = regression_vector_angle + pi;
    regression_vector_angle = mod(regression_vector_angle,2*pi);
end

% Check if angles are in alignment?
% Convert angles into positions on unit circle, so that we do not worry
% about wrap-around errors
projection_angle_position = [cos(projection_vector_angle) sin(projection_vector_angle)];
regression_angle_position = [cos(regression_vector_angle) sin(regression_vector_angle)];
distance_squared_in_unit_circle = sum((projection_angle_position - regression_angle_position).^2,2);

threshold_angle_squared = (20*pi/180)^2;
if distance_squared_in_unit_circle>threshold_angle_squared
    warning('on','backtrace');
    warning('Regression angle and fitted angle are signficantly different. ');
    warning('Projection angle is (degrees): %.3f',projection_vector_angle*180/pi);
    warning('Regression angle is (degrees): %.3f',regression_vector_angle*180/pi);
    error('Regression angle and fitted angle are so different that an error likely occurred!');
end

unit_converted_projection_vector = [cos(regression_vector_angle) sin(regression_vector_angle)];
unit_converted_orthogonal_vector = unit_converted_projection_vector*[0 1; -1 0];




% Project the base point onto the line segment using the dot
% product to create a new base point
projection_to_lowest_point = sorted_points_in_domain(1,:) - base_point_on_line;
transverse_distance_to_lowest_point = sum(projection_to_lowest_point.*unit_converted_projection_vector,2);
min_regression_point = base_point_on_line + transverse_distance_to_lowest_point*unit_converted_projection_vector;

projection_to_highest_point = sorted_points_in_domain(end,:) - base_point_on_line;
transverse_distance_to_highest_point = sum(projection_to_highest_point.*unit_converted_projection_vector,2);
max_regression_point = base_point_on_line + transverse_distance_to_highest_point*unit_converted_projection_vector;

% Find the transverse distances and standard deviations of
% transverse distances
projections = sorted_points_in_domain - base_point_on_line;
orthogonal_distances = sum(projections.*unit_converted_orthogonal_vector,2);

% Save point-to-point best-fit line segment
% regression_fit_line_segment = [min_regression_point; max_regression_point];
% [unit_projection_vectors(best_agreement_index,:) points(base_point_index,:) current_station_distances];
% regression_domain.best_fit_parameters = [min_regression_point max_regression_point];
regression_domain.best_fit_parameters = [unit_converted_projection_vector base_point_on_line transverse_distance_to_lowest_point transverse_distance_to_highest_point];

% Calculate the domain boxes
std_dev_orthogonal_distance = std(orthogonal_distances);
sigma_orthogonal_distance = std_dev_orthogonal_distance; % max(abs(orthogonal_distances));

N_sigmas = 1;
domain_box_1_sigma = ...
    [...
    min_regression_point - N_sigmas*sigma_orthogonal_distance*unit_converted_orthogonal_vector;
    max_regression_point - N_sigmas*sigma_orthogonal_distance*unit_converted_orthogonal_vector;
    max_regression_point + N_sigmas*sigma_orthogonal_distance*unit_converted_orthogonal_vector;
    min_regression_point + N_sigmas*sigma_orthogonal_distance*unit_converted_orthogonal_vector;
    ];

N_sigmas = 2;
domain_box_2_sigma = ...
    [...
    min_regression_point - N_sigmas*sigma_orthogonal_distance*unit_converted_orthogonal_vector;
    max_regression_point - N_sigmas*sigma_orthogonal_distance*unit_converted_orthogonal_vector;
    max_regression_point + N_sigmas*sigma_orthogonal_distance*unit_converted_orthogonal_vector;
    min_regression_point + N_sigmas*sigma_orthogonal_distance*unit_converted_orthogonal_vector;
    ];

N_sigmas = 3;
domain_box_3_sigma = ...
    [...
    min_regression_point - N_sigmas*sigma_orthogonal_distance*unit_converted_orthogonal_vector;
    max_regression_point - N_sigmas*sigma_orthogonal_distance*unit_converted_orthogonal_vector;
    max_regression_point + N_sigmas*sigma_orthogonal_distance*unit_converted_orthogonal_vector;
    min_regression_point + N_sigmas*sigma_orthogonal_distance*unit_converted_orthogonal_vector;
    ];

regression_domain.best_fit_1_sigma_box = polyshape(domain_box_1_sigma(:,1),domain_box_1_sigma(:,2),'Simplify',false,'KeepCollinearPoints',true);
regression_domain.best_fit_2_sigma_box = polyshape(domain_box_2_sigma(:,1),domain_box_2_sigma(:,2),'Simplify',false,'KeepCollinearPoints',true);
regression_domain.best_fit_3_sigma_box = polyshape(domain_box_3_sigma(:,1),domain_box_3_sigma(:,2),'Simplify',false,'KeepCollinearPoints',true);
regression_domain.best_fit_domain_box  = regression_domain.best_fit_2_sigma_box;


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

    % Get the color ordering?
    try
        color_ordering = orderedcolors('gem12');
    catch
        color_ordering = colororder;
    end

    N_colors = length(color_ordering(:,1));

    hold on;
    grid on;
    axis equal;

    % Plot the fits    
    ith_domain = 1;
    current_color = color_ordering(mod(ith_domain,N_colors)+1,:); %#ok<NASGU>
      
    % Plot the base point and vector
    plot(base_point_of_domain(1,1),base_point_of_domain(1,2),'g.','MarkerSize',30);
    quiver(base_point_of_domain(:,1),base_point_of_domain(:,2),unit_tangent_vector_of_domain(:,1),unit_tangent_vector_of_domain(:,2),0,'g','Linewidth',5);
        
    % Plot the domain
    fcn_geometry_plotFitDomains(regression_domain, fig_num);

    % Make axis slightly larger? And since this is the first one, save the
    % axis limits.
    if flag_rescale_axis
        temp = axis;
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        new_axis = [temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y];
        axis(new_axis);
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


