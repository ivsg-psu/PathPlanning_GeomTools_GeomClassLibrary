function [regression_fit_arc_center_and_radius_and_angles, domain_box, radial_errors, standard_deviation]  =  ...
    fcn_geometry_fitArcRegressionFromHoughFit(source_points, associated_points_in_domain, varargin)
% fcn_geometry_fitArcRegressionFromHoughFit
% Given a set of points that are matched to an arc via a Hough vote,
% finds the arc regression fit and domain box. 
% 
% Format: 
% [regression_fit_arc_center_and_radius_and_angles, domain_box] = fcn_geometry_fitArcRegressionFromHoughFit(source_points,associated_points_in_domain, (fig_num))
%
% INPUTS:
%      source_points: a 3x2 matrix of the points used to create the Hough circle fit (used to find direction): 
% 
%      [point1_x  point1_y; point2_x  point2_y; point3_x  point3_y;]
%
%      These points are used to determine the direction of the resulting
%      circle fit.
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
%      regression_fit_arc_center_and_radius_and_angles: a 5x1 matrix of the format:
%
%      [circleCenter_x  circleCenter_y circleRadius arc_start_angle_in_radians arc_end_angle_in_radians]
%
%      domain_box: the box that encloses the 2-standard-deviation interval
%      around the regression circle fit.
%
%      radial_errors: the individual errors in each point, radially, in an
%      [N x 1] matrix
%
%      standard_deviation: the standard deviation in the errors
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

if 0==flag_max_speed
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(2,3);

        % Check the source_points input to be length exactly equal to 3
        fcn_DebugTools_checkInputsToFunctions(...
            source_points, '2column_of_numbers',[3 3]);

        % Check the associated_points_in_domain input to be length greater
        % than or equal to 3
        fcn_DebugTools_checkInputsToFunctions(...
            associated_points_in_domain, '2column_of_numbers',[3 4]);
    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (3<= nargin) && (0==flag_max_speed)
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

% The solution approach is simple: fit with a circle and then find the start/end angles
[regression_fit_circle_center_and_radius, ~, radial_errors, standard_deviation] = ...
    fcn_geometry_fitCircleRegressionFromHoughFit(...
    [associated_points_in_domain(1,:); associated_points_in_domain(2,:); associated_points_in_domain(end,:)], ...
    associated_points_in_domain);


circleCenter = regression_fit_circle_center_and_radius(1,1:2); 
circleRadius = regression_fit_circle_center_and_radius(1,3);
index_source_point = 1;
station_tolerance = circleRadius*3*pi; % Allow the arc to go all the way around

[~, ~, start_angle_in_radians, end_angle_in_radians] = ...
    fcn_geometry_findArcAgreementIndicies(associated_points_in_domain, circleCenter, circleRadius, index_source_point, station_tolerance, -1);

% Check the direction
degree_step = 1;
is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(source_points(1,:), source_points(2,:), source_points(3,:));
if is_counterClockwise~=1
    degree_step = -1*degree_step;
    temp = start_angle_in_radians;
    start_angle_in_radians = end_angle_in_radians;
    end_angle_in_radians = temp;
end

% Fill in the results
regression_fit_arc_center_and_radius_and_angles = [circleCenter, circleRadius, start_angle_in_radians, end_angle_in_radians];

% Fill in domain box
sigma_multiplier = 2;
max_orthogonal_distance = sigma_multiplier*standard_deviation; % max(abs(orthogonal_distances));

% Create a domain by doing a large range of angles across an inner and
% outer arc that spans the test area

angles = (start_angle_in_radians:degree_step*pi/180:end_angle_in_radians)';
inner_radius = max(0,(circleRadius - max_orthogonal_distance));
outer_radius = circleRadius + max_orthogonal_distance;
inner_arc = inner_radius*[cos(angles) sin(angles)] + ones(length(angles(:,1)),1)*circleCenter;
outer_arc = outer_radius*[cos(angles) sin(angles)] + ones(length(angles(:,1)),1)*circleCenter;
domain_box = [inner_arc; flipud(outer_arc)];


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

    hold on;
    grid on;
    title('Regression fit of circle data');
    xlabel('X [meters]');
    ylabel('Y [meters]')
    
      
    % Plot the inputs: the source_points and associated_points_in_domain 
    plot(source_points(1,1),source_points(1,2),'g.','MarkerSize',30);
    plot(source_points(2,1),source_points(2,2),'b.','MarkerSize',30);
    plot(source_points(3,1),source_points(3,2),'r.','MarkerSize',30);
    h_plot = plot(associated_points_in_domain(:,1),associated_points_in_domain(:,2),'.','MarkerSize',10);
    current_color = get(h_plot,'Color');

    % Plot the circle fit
    fcn_geometry_plotCircle(circleCenter, circleRadius, current_color,fig_num)
    plot(circleCenter(1,1),circleCenter(1,2),'b+','MarkerSize',30);

    % Plot the domain
    domainShape = polyshape(domain_box(:,1),domain_box(:,2),'Simplify',false,'KeepCollinearPoints',true);
    plot(domainShape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',1,'EdgeAlpha',0);

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


