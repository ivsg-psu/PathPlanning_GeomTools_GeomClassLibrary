function [best_fit_parameters, std_dev_orthogonal_distance] = fcn_geometry_fitSegmentViaRegression(points_in_domain, varargin)
%% fcn_geometry_fitSegmentViaRegression 
% finds the regression fit vector that is the total least squares best fit
% of a set of points. The points are assumed to be ordered in the direction
% of the fit.
% 
% NOTE: the vector fit is not the same as a linear least-squares linear
% regression, which minimizes sum-of-squares of the VERTICAL errors between
% a line fit and the respective points. The method here is similar to
% total-least-squares. Namely, the vector is found that approximately
% minimizes the sum-of-squares distance between the vector and the
% orthogonal projection to each point. Thus, the fit is minimizing
% ORTHOGONAL errors. Unlike linear regression, this method works for data
% aligned vertically. The method is approximate because the 
% 
% Format: 
% [regression_fit_line_segment, domain_box] = fcn_geometry_fitSegmentViaRegression(Hough_domain, (fig_num))
%
% INPUTS:
%      points_in_domain: an semi-ordered set of points in the domain. Note:
%      the points are re-sorted using a projection vector of the first and
%      last point.
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      best_fit_parameters: a vector of fitted parameters for a segment type, namely:
%
%             [
%              base_point_x, 
%              base_point_y, 
%              heading,
%              s_Length,
%             ]
%
%      See fcn_geometry_fillEmptyDomainStructure for details.
%
%      std_dev: the standard deviation in the point
%      fit, as measured in the transverse direction (orthogonal to the line
%      fit). E.g., this is the total-least-squares standard deviation.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_calcUnitVector
%      fcn_geometry_fitVectorToNPoints
%      fcn_geometry_fitSlopeInterceptNPoints 
%      fcn_geometry_domainBoxByType
%      fcn_geometry_fillColorFromNumberOrName
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_fitSegmentViaRegression
% for a full test suite.
%
% This function was written on 2024_07_21 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2024_07_21 - S Brennan
% -- wrote the code using fcn_geometry_fitLinearRegressionFromHoughFit as a
% starter

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

% % Does user want to specify best_fit_domain_box_projection_distance?
% best_fit_domain_box_projection_distance = [];
% if (2<=nargin)
%     temp = varargin{1};
%     if ~isempty(temp)
%         best_fit_domain_box_projection_distance = temp;
%     end
% end

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


%% First, sort all the points in the approximate direction

% Find the unit vector (in fast mode) connecting base to end point.
approximate_base_point_of_domain          = points_in_domain(1,1:2);
approximate_end_point_of_domain           = points_in_domain(end,1:2);
approximate_vector_of_domain              = approximate_end_point_of_domain - approximate_base_point_of_domain;
approximate_unit_tangent_vector_of_domain = fcn_geometry_calcUnitVector(approximate_vector_of_domain);

% Find domain ranges related to point-to-point fit
approximate_projection_vectors_from_base = points_in_domain - approximate_base_point_of_domain;
approximate_tangent_distances = sum(approximate_projection_vectors_from_base.*approximate_unit_tangent_vector_of_domain,2);
% orthogonal_distances = sum(projection_vectors_from_base.*unit_orthogonal_vector_of_domain,2);

% Re-sort the points in direction of point-to-point projection
[~,sorted_indicies] = sort(approximate_tangent_distances);
resorted_points_in_domain = points_in_domain(sorted_indicies,:);

% What angle is the projection vector at?
approximate_projection_vector_angle = atan2(approximate_unit_tangent_vector_of_domain(:,2),approximate_unit_tangent_vector_of_domain(:,1));

%% Next, find best-fit line direction
% In the early form of this function, linear regression was used. The code
% for this is kept below in case we need to use it again sometime later.

flag_use_vector_regression = 1;
if flag_use_vector_regression
    [fitted_root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(resorted_points_in_domain,-1);
    fitted_base_point_on_line = fitted_root_point;
    regression_vector_angle = atan2(unit_vector(1,2),unit_vector(1,1));

else
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(resorted_points_in_domain,-1);

    % Convert slope into projection and orthogonal vectors
    regression_vector_angle = atan2(slope,1);

    % Find new base point that uses regression fit
    y_value_on_line = slope*approximate_base_point_of_domain(1,1)+intercept;
    fitted_base_point_on_line = [approximate_base_point_of_domain(1,1) y_value_on_line];

end

% Fix both angles to 0 to 2pi range
approximate_projection_vector_angle = mod(approximate_projection_vector_angle,2*pi);

% Make sure angle is between 0 and 2*pi
regression_vector_angle = mod(regression_vector_angle,2*pi);

% Check if regression slope angle was calculated in wrong direction. If so,
% the vector should be approximately pi off from the correct direction.
% So if this is the case, add pi onto it to fix
if abs(approximate_projection_vector_angle - regression_vector_angle)>pi/2
    if flag_do_debug
        fprintf(1,'Correcting angle...\n')
    end
    regression_vector_angle = regression_vector_angle + pi;

    % Make sure angle is between 0 and 2*pi
    regression_vector_angle = mod(regression_vector_angle,2*pi);
end


unit_converted_regression_vector = [cos(regression_vector_angle) sin(regression_vector_angle)];
unit_converted_orthogonal_vector = unit_converted_regression_vector*[0 1; -1 0];

%% Find the new base point for the segment 
% Find transverse distances for each point and keep only the lowest and
% highest. The lowest defines the base point. The difference defines the
% length.
vectors_from_base_point_to_each_point = resorted_points_in_domain - fitted_base_point_on_line;
transverse_distance_to_each_point = sum(vectors_from_base_point_to_each_point.*unit_converted_regression_vector,2);
min_transverse_distance = min(transverse_distance_to_each_point);
max_transverse_distance = max(transverse_distance_to_each_point);

updated_base_point_on_line  = fitted_base_point_on_line + min_transverse_distance*unit_converted_regression_vector;

updated_segment_angle = regression_vector_angle;
updated_segment_length = max_transverse_distance - min_transverse_distance;
best_fit_parameters = [updated_base_point_on_line updated_segment_angle updated_segment_length];


%% Calculate the standard deviation of the fit
% Find the transverse distances and standard deviations of
% transverse distances
orthogonal_distances = sum(vectors_from_base_point_to_each_point.*unit_converted_orthogonal_vector,2);
std_dev_orthogonal_distance = std(orthogonal_distances);


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

    hold on;
    grid on;
    axis equal;

    % Plot the input points
    plot(resorted_points_in_domain(:,1),resorted_points_in_domain(:,2),'k.','MarkerSize',10);

    % Plot the fits    
    ith_domain = 1;
    current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain); 
     
    fcn_geometry_plotGeometry('segment',best_fit_parameters, [],current_color);

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


