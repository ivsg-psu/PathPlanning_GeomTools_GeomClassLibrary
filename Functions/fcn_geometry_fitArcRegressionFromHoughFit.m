function [regression_domain, std_dev_orthogonal_distance] = fcn_geometry_fitArcRegressionFromHoughFit(Hough_domain, varargin)
%% fcn_geometry_fitArcRegressionFromHoughFit
% Given a domain containing a set of points that are matched to an arc via
% a Hough vote, finds the arc regression fit and domain box.
% 
% Format: 
% [regression_domain, std_dev_transverse_distance] = fcn_geometry_fitArcRegressionFromHoughFit(Hough_domain, (best_fit_domain_box_projection_distance), (fig_num))
%
% INPUTS:
%      Hough_domain: a structure that records details of the domain of
%      fitting. See fcn_geometry_fillEmptyDomainStructure for details.
%
%      (OPTIONAL INPUTS)
%
%      best_fit_domain_box_projection_distance: the distance from the curve
%      fit, in the transverse direction, to project in both the positive
%      and negative directions to produce the best_fit_domain_box. If left
%      empty, defaults to 2 standard deviations to thus give a box that is
%      +/- 2 sigma.
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
%      fcn_geometry_fitCircleRegressionFromHoughFit
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_fitArcRegressionFromHoughFit
% for a full test suite.
%
% This function was written on 2024_01_09 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2024_01_09 - S Brennan
% -- wrote the code
% 2024_01_18 - S Brennan
% -- changed to domain inputs and outputs
% 2024_04_02 - S Brennan
% -- added best_fit_domain_box_projection_distance as an input option

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
        narginchk(1,3);

        % % Check the source_points input to be length exactly equal to 3
        % fcn_DebugTools_checkInputsToFunctions(...
        %     source_points, '2column_of_numbers',[3 3]);
        %
        % % Check the associated_points_in_domain input to be length greater
        % % than or equal to 3
        % fcn_DebugTools_checkInputsToFunctions(...
        %     associated_points_in_domain, '2column_of_numbers',[3 4]);
    end
end

% Does user want to specify best_fit_domain_box_projection_distance?
best_fit_domain_box_projection_distance = [];
if (2<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        best_fit_domain_box_projection_distance = temp;
    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if  (0==flag_max_speed) && (3<= nargin)
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
Hough_best_fit_source_indicies = Hough_domain.best_fit_source_indicies;
% flag_arc_is_counterclockwise = Hough_domain.best_fit_parameters(1,7);

% Check to make sure we have a Hough fit input
switch best_fit_type
    case {'Hough circle'}
        regression_domain.best_fit_type = 'Regression circle';
        Hough_flag_this_is_a_circle  = 1;
    case {'Hough arc'}
        regression_domain.best_fit_type = 'Regression arc';
        Hough_flag_this_is_a_circle  = Hough_domain.best_fit_parameters(1,6);
    otherwise
        error('A domain was tested for fitting that was not a Hough fit');
end
regression_domain.points_in_domain = points_in_domain;


% Pull out parameters from the Hough domain fit
associated_points_in_domain = points_in_domain;
source_points = [...
    associated_points_in_domain(Hough_best_fit_source_indicies(1),:); 
    associated_points_in_domain(Hough_best_fit_source_indicies(2),:); 
    associated_points_in_domain(Hough_best_fit_source_indicies(3),:)];

regression_domain.points_in_domain = points_in_domain;

%% The solution approach is simple: fit with a circle and then find the start/end angles

% Fit the circle portion
[regression_fit_circle_center_and_radius, ~, ~, standard_deviation] = ...
    fcn_geometry_fitCircleRegressionFromHoughFit(source_points, associated_points_in_domain);

circleCenter = regression_fit_circle_center_and_radius(1,1:2);
circleRadius = regression_fit_circle_center_and_radius(1,3);

% Calculate the remaining arc details?
if Hough_flag_this_is_a_circle
    regression_domain.best_fit_parameters = [circleCenter, circleRadius];
    start_angle_in_radians = 0;
    end_angle_in_radians = 2*pi;
    degree_step = 1;
else
    % Find the arc angles associated with the regression fit
    % These may be slightly different from Hough fit because the circle center
    % and radius will be different, due to regression fit
    index_source_point = Hough_best_fit_source_indicies(1);
    station_tolerance = circleRadius*3*pi; % Allow the arc to go all the way around, we only do this to force the check below to use all the points

    [~, ~, start_angle_in_radians, end_angle_in_radians] = ...
        fcn_geometry_findArcAgreementIndicies(associated_points_in_domain, circleCenter, circleRadius, index_source_point, station_tolerance, -1);

    % Check the direction. If direction is not aligned with the input regression direction, flip it
    degree_step = min(1,(end_angle_in_radians-start_angle_in_radians)*1/10*180/pi);
    % is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(source_points(1,:), source_points(2,:), source_points(3,:));
    flag_arc_is_counterclockwise = fcn_geometry_arcDirectionFromCircleCenter(associated_points_in_domain, circleCenter, -1);

    if flag_arc_is_counterclockwise~=1
        degree_step = -1*degree_step;
        temp = start_angle_in_radians;
        start_angle_in_radians = end_angle_in_radians;
        end_angle_in_radians = temp;
    end

    % Fill in the results
    regression_fit_arc_center_and_radius_and_angles_and_flag = [circleCenter, circleRadius, start_angle_in_radians, end_angle_in_radians, Hough_flag_this_is_a_circle, flag_arc_is_counterclockwise];
    regression_domain.best_fit_parameters = regression_fit_arc_center_and_radius_and_angles_and_flag;
end

%% Calculate the domain boxes
% Create a domain by doing a range of angles across inner and
% outer arc that spans the test area

N_points = ceil(abs(end_angle_in_radians - start_angle_in_radians)/abs(degree_step*pi/180));
angles = linspace(start_angle_in_radians, end_angle_in_radians,N_points)';

% Calculate the domain boxes
std_dev_orthogonal_distance = standard_deviation;
sigma_orthogonal_distance = std_dev_orthogonal_distance; % max(abs(orthogonal_distances));

domain_box_1_sigma = fcn_geometry_domainBoxByType('arc', circleCenter, circleRadius, angles,  1*sigma_orthogonal_distance,-1);
domain_box_2_sigma = fcn_geometry_domainBoxByType('arc', circleCenter, circleRadius, angles,  2*sigma_orthogonal_distance,-1);
domain_box_3_sigma = fcn_geometry_domainBoxByType('arc', circleCenter, circleRadius, angles,  3*sigma_orthogonal_distance,-1);

regression_domain.best_fit_1_sigma_box = domain_box_1_sigma;
regression_domain.best_fit_2_sigma_box = domain_box_2_sigma;
regression_domain.best_fit_3_sigma_box = domain_box_3_sigma;
if isempty(best_fit_domain_box_projection_distance)
    regression_domain.best_fit_domain_box  = regression_domain.best_fit_2_sigma_box;
else
    additional_arc_radians = best_fit_domain_box_projection_distance/circleRadius;
    if end_angle_in_radians>start_angle_in_radians
        direction_of_angles = 1;
    else
        direction_of_angles = -1;
    end
    angles_padded = [start_angle_in_radians - direction_of_angles*additional_arc_radians; angles; end_angle_in_radians + direction_of_angles*additional_arc_radians];
    regression_domain.best_fit_domain_box  = fcn_geometry_domainBoxByType('arc', circleCenter, circleRadius, angles_padded,  best_fit_domain_box_projection_distance,-1);
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

    % Plot the results in point space
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
    current_color = color_ordering(mod(ith_domain,N_colors)+1,:); 
      
    
    if flag_do_debug
        % Plot the source_points
        plot(source_points(1,1),source_points(1,2),'g.','MarkerSize',30);
        plot(source_points(2,1),source_points(2,2),'b.','MarkerSize',30);
        plot(source_points(3,1),source_points(3,2),'r.','MarkerSize',30);
    end

    % Plot the associated_points_in_domain
    plot(associated_points_in_domain(:,1),associated_points_in_domain(:,2),'.','MarkerSize',5,'Color',current_color);

    % Plot the domains
    fcn_geometry_plotFitDomains(regression_domain, fig_num);

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


