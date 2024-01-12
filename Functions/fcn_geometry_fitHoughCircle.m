function [best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies, best_fit_is_a_circle] = ...
    fcn_geometry_fitHoughCircle(points, transverse_tolerance, varargin)
% fcn_geometry_fitHoughCircle
%
% This function takes the input points and tolerance as the input and
% outputs the fitted parameters and agreement indices. 
% 
% [best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies, best_fit_is_a_circle]  = fcn_geometry_fitHoughCircle(points, transverse_tolerance, ...
%         (station_tolerance), (expected_radii_range), (flag_use_permutations), (fig_num))
% 
% INPUTS:
%      points: a Nx2 vector where N is the number of points, but at least 2 rows. 
%
%      transverse_tolerance: the orthogonal distance between the points and
%      the linear curve fit that indicate whether a point "belongs" to the
%      fit (if distance is less than or equal to the tolerance), or is
%      "outside" the fit (if distance is greater than the tolerance).
%
%      (OPTIONAL INPUTS)
%
%      station_tolerance: the projection distance between the points in a
%      curve fit, along the direction of the line, that indicate whether a
%      point "belongs" to the linear fit (if distance is less than or equal
%      to the tolerance), or is "outside" the fit (if distance is greater
%      than the tolerance). If left empty, then a circle fit is performed
%      using only the transverse_tolerance.
% 
%      expected_radii_range: a vector in form of [r_min r_max] indicating
%      expected radius. Any radii outside this range will not be assessed. 
% 
%      flag_use_permutations: specify permutation type. Can be set to:
%
%          0:search for permutations via ordering 1-2-3, then 2-3-4, etc.
%          This assumes best-fit circle will be formed by points in direct
%          sequence, and the points are in order (VERY fast)
%
%          1:search for permutations via nchoosek, which assumes best-fit
%          circle will not be formed by points not in direct sequence, but
%          the points are in order. (somewhat slow, especially in number of
%          points is greater than 100).  
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%      Let K represent the number of permutations possible in N-choose-2
%      standard ordering (see function nchoosek). For each K value, this
%      returns:
%
%      best_fitted_parameters: the best-fitted circle or arc in format of
%      [x_center, y_center, radius, start_angle_in_radians,
%      end_angle_in_radians] for the best 3-point combination of the point
%      permutation. This is returned as a [1 x 5] matrix
%
%      best_fit_source_indicies: the three indicies, in [1x3] vector format,
%      of the points that produced the best fit.
% 
%      best_agreement_indicies: the indicies of the points that are within
%      agreement of the best-fit parameters, given the transverse and
%      station tolerance settings.
%
%      best_fit_is_a_circle: returns 1 if the best fit is a circle, 0 if
%      the best fit was an arc. Note: it will only perform circle-fitting
%      if the station_tolerance is left blank.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_circleCenterFrom3Points
%      fcn_geometry_plotCircle 
%      fcn_geometry_findArcAgreementIndicies
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_fitHoughCircle
% for a full test suite.
%
% This function was written on 2023_12_15 by Aneesh Batchu
% Questions or comments? abb6486@psu.edu or sbrennan@psu.edu

% Revision history:
% 2023_12_15 - A. Batchu
% -- wrote the code
% 2023_12_18 - S. Brennan
% -- fixed broken for-loop
% -- moved circle fitting to external call to existing function
% -- added transverse and station tolerances
% -- added plotting figure for results as varargin
% 2023_12_27 - S. Brennan
% -- changed inequality-based region testing with polygon region tests
% 2023_12_28 - S. Brennan
% -- added fast mode by allowing fig_num set to -1
% 2024_01_03 - S. Brennan
% -- output only the best fits, keeping same format as fitHoughLine code
% 2024_01_05 - S. Brennan
% -- allow user to specify circle radii to search over
% -- allow user to specify search sequence
% -- added station tolerance to enable arc-fits


%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==6 && isequal(varargin{end},-1))
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
    debug_fig_num = 234343;
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

if 0==flag_max_speed
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(2,6);

        % Check the points input to be length greater than or equal to 2
        fcn_DebugTools_checkInputsToFunctions(...
            points, '2column_of_numbers',[2 3]);

        % Check the transverse_tolerance input is a positive single number
        fcn_DebugTools_checkInputsToFunctions(transverse_tolerance, 'positive_1column_of_numbers',1);

    end
end


% Does user want to specify station_tolerance?
station_tolerance = [];
if (4<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        station_tolerance = temp;
    end
end

% Does user want to specify expected_radii_range?
expected_radii_range = [-inf inf];
if (4<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        expected_radii_range = temp;
    end
end

% Does user want to specify flag_use_permutations?
flag_use_permutations = 1;
if (5<=nargin)
    temp = varargin{3};
    if ~isempty(temp)
        flag_use_permutations = temp;
    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (6<=nargin) && (0==flag_max_speed)
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


% Total number of points
N_points = size(points,1);

% Do debugging plots?
if flag_do_debug
    figure(debug_fig_num);
    clf;
    hold on;
    grid on;
    axis equal

    plot(points(:,1),points(:,2),'k.','MarkerSize',20);
    temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;

    debug_axis_limits = [temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y];
    axis(debug_axis_limits);
end


% Find all possible 3-point combinations
% NOTE: to find the permutations, use perms (for example: perms(1:3)). This
% will return all the ordering of the columns used in nchoosek that can be
% tried.
switch flag_use_permutations
    case 0
        combos_paired = [(1:N_points-2)' (2:N_points-1)' (3:N_points)'];
    case 1
        combos_paired = nchoosek(1:N_points,3);
    otherwise
        if isinf(flag_use_permutations)
            error('Cannot set flag_use_permutations to infinity.');
        elseif flag_use_permutations>0
            combos_paired = [(1:N_points-2)' (2:N_points-1)' (3:N_points)'];
            random_combos = rand(N_points-2,1);

            % Throw a random number
            if flag_use_permutations<1
                % check if random number for each combo is less than the
                % value given by user - if so, keep it!
                keeps = random_combos<=flag_use_permutations;
                combos_paired = combos_paired(keeps,:);
            else

                % sort the random numbers, keeping only top N of them
                N_keeps = ceil(flag_use_permutations);
                N_keeps = min([N_keeps length(combos_paired(:,1))]);
                [~,indicies_sorted] = sort(random_combos);
                combos_paired = combos_paired(indicies_sorted(1:N_keeps),:);
            end
        else
            error('Invalid setting for flag_use_permutations')
        end
end

% How many combinations are there?
N_permutations = size(combos_paired,1);

% Pre-allocation of fittedParameters and agreementIndices for saving
% computation time
fitted_parameters = zeros(N_permutations,3);

% Initialize output vectors
best_agreement_indicies_binary_form = zeros(1,N_points);
best_agreement_count = -inf;
best_agreement_index = [];
best_start_angle_in_radians = 0;
best_end_angle_in_radians = 0;
best_fit_is_a_circle = 0;


agreements = zeros(N_permutations,1);

for ith_combo = 1:N_permutations

    % fprintf(1,'Checking %.0d of %.0d\n',ith_combo,N_permutations);
    if 0==flag_max_speed
        if 0==mod(ith_combo,100)
            fprintf(1,'Checking %.0d of %.0d\n',ith_combo,N_permutations);
        end
    end

    % Extract the source points
    test_source_points = points(combos_paired(ith_combo,:),:);

    % Find fitted curve - call the function in "fast" mode
    [circleCenter, circleRadius] = fcn_geometry_circleCenterFrom3Points(test_source_points(1,:),test_source_points(2,:),test_source_points(3,:),-1);
    
    % Store resulting fitted parameters in "fitted_parameters" matrix
    fitted_parameters(ith_combo,:) = [circleCenter, circleRadius];

    if (circleRadius>=expected_radii_range(1)) && (circleRadius<=expected_radii_range(2))
        % Find points in agreement?
        if flag_do_debug
            % Do region-based testing (this is very slow)
            [agreement_indicies, domainPolyShape, start_angle_in_radians, end_angle_in_radians, flag_is_a_circle] = ...
                fcn_INTERNAL_findAgreementIndicies(points, circleCenter, circleRadius, combos_paired(ith_combo,:), transverse_tolerance, station_tolerance, best_agreement_count, 0);
        else
            % Do inequality-based testing (fast)
            [agreement_indicies, domainPolyShape, start_angle_in_radians, end_angle_in_radians, flag_is_a_circle] = ...
                fcn_INTERNAL_findAgreementIndicies(points, circleCenter, circleRadius, combos_paired(ith_combo,:), transverse_tolerance, station_tolerance, best_agreement_count, 1);
        end
    else
        agreement_indicies = [];
        domainPolyShape = [];
        start_angle_in_radians = [];
        end_angle_in_radians = [];
        flag_is_a_circle = 0;
    end

    
    % Save the count of points in agreement
    agreement_count = length(agreement_indicies);
    agreements(ith_combo) = agreement_count;

    % Check if this is the best agreement so far
    if agreement_count>best_agreement_count
        best_agreement_index = ith_combo;

        % Save new "best" agreement total count
        best_agreement_count = agreement_count;
        
        % Save angles
        best_start_angle_in_radians = start_angle_in_radians;
        best_end_angle_in_radians = end_angle_in_radians;

        % Save new "best" indicies
        best_agreement_indicies_binary_form = zeros(1,N_points);
        best_agreement_indicies_binary_form(1,agreement_indicies) = 1;

        % Save if it is a circle
        best_fit_is_a_circle = flag_is_a_circle;
    end



    % Redo debugging plots to see how fit worked?
    if flag_do_debug
        figure(debug_fig_num);
        clf;
        hold on;
        grid on;
        axis equal

        % Plot the points in black        
        plot(points(:,1),points(:,2),'k.','MarkerSize',20);

        % Plot the domain
        plot(domainPolyShape);

        % Plot the points in agreement
        plot(points(agreement_indicies,1),points(agreement_indicies,2),'r.','MarkerSize',20);

        % Plot the test points in blue
        plot(test_source_points(:,1),test_source_points(:,2),'b.','MarkerSize',20);

        % Use the same axis each time, to make plots stay in 'same' place
        axis(debug_axis_limits);
    end

end


% Save results to outputs
best_fitted_parameters = [fitted_parameters(best_agreement_index,:) best_start_angle_in_radians best_end_angle_in_radians];
best_fit_source_indicies = combos_paired(best_agreement_index,:);
best_agreement_indicies = find(best_agreement_indicies_binary_form==1);


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
    title('Points and maximum-vote fit, plotted in point-space');
    xlabel('X [meters]');
    ylabel('Y [meters]')

    [~,sorted_indicies] = sort(agreements,'ascend');


    % Plot the input points
    plot(points(:,1),points(:,2),'k.','MarkerSize',20);

    % Plot the circle fit
    if best_fit_is_a_circle==1
        title('Circle fit');
    else
        title('Arc fit');
    end

    fcn_geometry_plotCircle(best_fitted_parameters(1:2), best_fitted_parameters(3), 'b-',fig_num)
    fcn_geometry_plotCircle(best_fitted_parameters(1:2), best_fitted_parameters(3)-transverse_tolerance, 'r-',fig_num) 
    fcn_geometry_plotCircle(best_fitted_parameters(1:2), best_fitted_parameters(3)+transverse_tolerance, 'r-',fig_num) 
    plot(best_fitted_parameters(1),best_fitted_parameters(2),'b+','MarkerSize',15);

    % Plot the best-fit points
    plot(points(best_agreement_indicies,1),points(best_agreement_indicies,2),'r.','MarkerSize',15);

    % Label the points
    text(points(best_agreement_indicies(end),1),points(best_agreement_indicies(end),2),sprintf('Start angle: %.3f deg', best_start_angle_in_radians*180/pi));
    text(points(best_agreement_indicies(1),1),points(best_agreement_indicies(1),2),sprintf('End angle: %.3f deg', best_end_angle_in_radians*180/pi));

    % Plot the source points
    plot(points(best_fit_source_indicies,1),points(best_fit_source_indicies,2),'bo','MarkerSize',15);

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end

    % Plot the Hough space results, from least to best
    figure(fig_num+1);
    clf;

    subplot(3,1,1);
    hold on;
    grid on;
    xlabel('Rho [radians]');
    ylabel('Distance of circle center from origin [meters]');
    xlim([-pi pi]);
    ylim([0 20]);

    subplot(3,1,2);
    hold on;
    grid on;
    xlabel('Rho [radians]');
    ylabel('Curvature [1/meters]');
    xlim([-pi pi]);
    ylim([0 1]);

    subplot(3,1,3);
    hold on;
    grid on;
    xlabel('Rho [radians]');
    ylabel('Distance of circle center from origin [meters]');
    zlabel('Curvature [1/meters]');
    xlim([-pi pi]);
    ylim([0 20]);
    zlim([0 1]);
    view(3);


    fitted_parameters(ith_combo,:) = [circleCenter, circleRadius];
    rho = atan2(fitted_parameters(:,2),fitted_parameters(:,1));
    distance_circle_center_from_origin = sum((fitted_parameters(:,1).^2+fitted_parameters(:,2).^2),2).^0.5;
    curvature = 1./fitted_parameters(:,3);

    % Plot the results in increasing order of better fit
    N_steps = 20;    
    indicies = linspace(1,N_permutations,N_steps);
    for ith_plot = 1:N_steps-1
        plot_indicies_start = ceil(indicies(ith_plot));
        plot_indicies_end = floor(indicies(ith_plot+1));

        sorted_indicies_to_plot = sorted_indicies(plot_indicies_start:plot_indicies_end);
        plot_color = (N_steps - ith_plot + 1)/N_steps*[1 1 1];
        marker_size = ceil(ith_plot*20/N_steps);

        subplot(3,1,1);
        plot(rho(sorted_indicies_to_plot,1),distance_circle_center_from_origin(sorted_indicies_to_plot,1),'.','MarkerSize',marker_size,'Color',plot_color);

        subplot(3,1,2);
        plot(rho(sorted_indicies_to_plot,1),curvature(sorted_indicies_to_plot,1),'.','MarkerSize',marker_size,'Color',plot_color);

        subplot(3,1,3);
        plot3(rho(sorted_indicies_to_plot,1),distance_circle_center_from_origin(sorted_indicies_to_plot,1),curvature(sorted_indicies_to_plot,1),'.','MarkerSize',marker_size,'Color',plot_color);

    end

    % Plot the best fits
    best_rho = atan2(best_fitted_parameters(:,2),best_fitted_parameters(:,1));
    best_distance_circle_center_from_origin = sum((best_fitted_parameters(:,1).^2+best_fitted_parameters(:,2).^2),2).^0.5;
    best_curvature = 1./best_fitted_parameters(:,3);

    subplot(3,1,1);
    plot(best_rho(1,1),best_distance_circle_center_from_origin(1,1),'.','MarkerSize',30,'Color',[1 0 0]);

    subplot(3,1,2);
    plot(best_rho(1,1),best_curvature(1,1),'.','MarkerSize',30,'Color',[1 0 0]);

    subplot(3,1,3);
    plot3(best_rho(1,1),best_distance_circle_center_from_origin(1,1),best_curvature(1,1),'.','MarkerSize',30,'Color',[1 0 0]);

    
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

%% fcn_INTERNAL_findAgreementIndicies
function [agreement_indicies, domainPolyShape, start_angle_in_radians, end_angle_in_radians, flag_is_a_circle] = ...
    fcn_INTERNAL_findAgreementIndicies(points, circleCenter, circleRadius, test_source_point_indicies, transverse_tolerance, station_tolerance, best_agreement_count, flag_max_speed) %#ok<INUSD>


flag_do_debug = 0;

% Set a flag to use enclosed domains to find points in agreement, rather
% than inequalities. This method is more general and easier to debug, but
% probably is slower.
if 0==flag_max_speed
    % Slow method
    flag_use_domain_method = 1;
else
    % Fast method
    flag_use_domain_method = 0;
end
domainPolyShape = []; 

% Create a plot for debugging?
if flag_do_debug
    figure(38383);
    clf;
    hold on;
    grid on;

    % Plot the input points
    plot(points(:,1),points(:,2),'k.','MarkerSize',20);

    % Plot the source points
    test_source_points = points(test_source_point_indicies,:);
    plot(test_source_points(:,1),test_source_points(:,2),'bo','MarkerSize',15);
    plot(test_source_points(1,1),test_source_points(1,2),'b.','MarkerSize',30);

    % Plot the circle fit
    fcn_geometry_plotCircle(circleCenter, circleRadius, 'b-',38383) 
    fcn_geometry_plotCircle(circleCenter, circleRadius-transverse_tolerance, 'r-',38383) 
    fcn_geometry_plotCircle(circleCenter, circleRadius+transverse_tolerance, 'r-',38383) 
    plot(circleCenter(:,1),circleCenter(:,2),'b+','MarkerSize',15);


    % Make axis slightly larger?
    temp = axis;
    %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

end



% Find agreement domain
if flag_use_domain_method

    % Create a domain by doing a large range of angles across an inner and
    % outer arc that spans the test area
    angles = (0:1:359.9999)'*pi/180;
    inner_radius = max(0,(circleRadius - transverse_tolerance));
    outer_radius = circleRadius + transverse_tolerance;
    inner_arc = inner_radius*[cos(angles) sin(angles)] + ones(length(angles(:,1)),1)*circleCenter;
    outer_arc = outer_radius*[cos(angles) sin(angles)] + ones(length(angles(:,1)),1)*circleCenter;
    best_fit_domain_box = [inner_arc; flipud(outer_arc)];

    % Find the points within the domain using isinterior
    domainPolyShape = polyshape(best_fit_domain_box(:,1),best_fit_domain_box(:,2),'Simplify',false,'KeepCollinearPoints',true);
    indicies_in_transverse_agreement = find(isinterior(domainPolyShape,points) == 1);

else
    % Distance of all the input points from the center
    point_radii = sum((points - circleCenter).^2,2).^0.5;

    % Absolute Error to find the indices in agreement
    absolute_radial_error = abs(point_radii - circleRadius);

    % Indices in transverse agreement
    indicies_in_transverse_agreement_binary_form = ((absolute_radial_error <= transverse_tolerance))' ;
    indicies_in_transverse_agreement = find(indicies_in_transverse_agreement_binary_form==1);
end


if flag_do_debug
    figure(38383);
   
    % Plot the radial agreement points
    plot(points(indicies_in_transverse_agreement,1),points(indicies_in_transverse_agreement,2),'r.','MarkerSize',15);
end

% Find the indicies in station agreement? Only do this if the transverse
% count is larger than current best count (this needs to be debugged)
if 1==1 %length(indicies_in_transverse_agreement)<best_agreement_count 
    if ~isempty(station_tolerance)
        % Grab only the points in radial agreement
        points_in_radial_agreement = points(indicies_in_transverse_agreement,:);

        % Find index of the source point
        index_source_point = find(indicies_in_transverse_agreement == test_source_point_indicies(1),1);

        % Call a function to find which points are in station agreement
        % also
        [indicies_in_station_agreement, flag_is_a_circle, start_angle_in_radians, end_angle_in_radians] = ...
            fcn_geometry_findArcAgreementIndicies(points_in_radial_agreement, circleCenter, circleRadius, index_source_point, station_tolerance, -1);

        % Calculate agreement with both vectors
        % agreement_indicies_binary_form = indicies_in_transverse_agreement.*points_within_angle;
        agreement_indicies = indicies_in_transverse_agreement(indicies_in_station_agreement);

    else
        agreement_indicies = indicies_in_transverse_agreement;
        flag_is_a_circle = 1;
        start_angle_in_radians = 0;
        end_angle_in_radians = 2*pi;
    end
else
    agreement_indicies = indicies_in_transverse_agreement;
    flag_is_a_circle = [];
    start_angle_in_radians = [];
    end_angle_in_radians = [];
end

if flag_do_debug
    figure(38383);
   
    % Plot the radial agreement points
    plot(points(agreement_indicies,1),points(agreement_indicies,2),'g.','MarkerSize',10);
end

end % Ends fcn_INTERNAL_findAgreementIndicies
