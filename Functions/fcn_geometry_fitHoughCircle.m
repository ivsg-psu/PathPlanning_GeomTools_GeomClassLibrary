function domains = fcn_geometry_fitHoughCircle(points, transverse_tolerance, varargin)
% fcn_geometry_fitHoughCircle
%
% This function takes the input points and tolerance as the input and
% outputs the fitted parameters and agreement indices of the top-voted
% circle.
% 
% [best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies, best_fit_is_a_circle]  = fcn_geometry_fitHoughCircle(points, transverse_tolerance, ...
%         (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_use_permutations), (fig_num))
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
%      point "belongs" to the circle fit (if distance is less than or equal
%      to the tolerance), or is "outside" the fit (if distance is greater
%      than the tolerance). If left empty, then a circle fit is performed
%      using only the transverse_tolerance.
%
%      points_required_for_agreement: the number of points required for an
%      agreement to be valid, with minimum value of 3. If left empty, then
%      the best agreement will always be returned. If a value is given,
%      line fitting will continue by clustering data until there are no
%      fits greater than or equal to points_required_for_agreement. The
%      results of the line fit will be saved in a cell array.
% 
%      flag_force_circle_fit: specify that only circle fits are allowed.
%      Can be set to:
%
%          0:search for both arcs and circles (default)
%
%          1:search only for circle fits. If station tolerance is given,
%          then a circle will only be returned if all points meet both the
%          transverse and station tolerances. Note: to force a circle fit
%          irregardless of station tolerance, set station_tolerance to an
%          empty value, e.g. station_tolerance = [];
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
%          1: (default) search for permutations via nchoosek, which assumes best-fit
%          circle will not be formed by points not in direct sequence, but
%          the points are in order. (somewhat slow, especially in number of
%          points is greater than 100).  
% 
%          N: search for permutations by randomly selecting N of them to
%          test. This essentially turns the Hough transform vote into
%          RANSAC.
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      domains: a structure that records details of the domain of fitting.
%      See fcn_geometry_fillEmptyDomainStructure for details.
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
% 2024_01_12 - S. Brennan
% -- added flag to force circle fitting, so allow circle searching even if
% station tolerance is given. Otherwise, both arcs and circles could be
% returned.
% 2024_01_15 - S. Brennan
% -- changed outputs to domain types
% -- added fcn_geometry_plotFitDomains for plotting

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==8 && isequal(varargin{end},-1))
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
        narginchk(2,8);

        % Check the points input to be length greater than or equal to 2
        fcn_DebugTools_checkInputsToFunctions(...
            points, '2column_of_numbers',[2 3]);

        % Check the transverse_tolerance input is a positive single number
        fcn_DebugTools_checkInputsToFunctions(transverse_tolerance, 'positive_1column_of_numbers',1);

    end
end


% Does user want to specify station_tolerance?
station_tolerance = [];
if (3<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        station_tolerance = temp;
    end
end


% Does user specify points_required_for_agreement?
points_required_for_agreement = 1;
if (4<= nargin)
    temp = varargin{2};
    if ~isempty(temp)
        points_required_for_agreement = temp;
        if points_required_for_agreement<3
            error('The input points_required_for_agreement must be greater than or equal to 3.')
        end
    end
end

% Does user want to specify expected_radii_range?
flag_force_circle_fit = 0;
if (5<=nargin)
    temp = varargin{3};
    if ~isempty(temp)
        flag_force_circle_fit = temp;
    end
end

% Does user want to specify expected_radii_range?
expected_radii_range = [-inf inf];
if (6<=nargin)
    temp = varargin{4};
    if ~isempty(temp)
        expected_radii_range = temp;
    end
end

% Does user want to specify flag_use_permutations?
flag_use_permutations = 1;
if (7<=nargin)
    temp = varargin{5};
    if ~isempty(temp)
        flag_use_permutations = temp;
    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (8<=nargin)
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


%% Step 1: find all the agreement counts, save in array "agreements"
agreements = zeros(N_permutations,1);

for ith_combo = 1:N_permutations

    % fprintf(1,'Checking %.0d of %.0d\n',ith_combo,N_permutations);
    if 0==flag_max_speed
        if 0==mod(ith_combo,10000)
            fprintf(1,'Checking %.0d of %.0d\n',ith_combo,N_permutations);
        end
    end

    % Extract the source points
    test_source_points = points(combos_paired(ith_combo,:),:);

    % Find fitted curve - call the function in "fast" mode
    [circleCenter, circleRadius] = fcn_geometry_circleCenterFrom3Points(test_source_points(1,:),test_source_points(2,:),test_source_points(3,:),-1);

    % Store resulting fitted parameters in "fitted_parameters" matrix
    fitted_parameters(ith_combo,:) = [circleCenter, circleRadius];

    % Check results only in cases that can possibly be valid fits
    if (circleRadius>=expected_radii_range(1)) && (circleRadius<=expected_radii_range(2))
        % Find points in agreement?
        agreement_indicies = ...
            fcn_geometry_findAgreementsOfPointsToArc(points,  combos_paired(ith_combo,1), circleCenter, circleRadius, transverse_tolerance, (station_tolerance), (flag_force_circle_fit),(points_required_for_agreement), (-1));

     end

    
    % Save the count of points in agreement
    agreement_count = length(agreement_indicies);
    agreements(ith_combo) = agreement_count;

    % Redo debugging plots to see how fit worked?
    if flag_do_debug
        figure(debug_fig_num);
        clf;
        hold on;
        grid on;
        axis equal

        % Plot the points in black        
        plot(points(:,1),points(:,2),'k.','MarkerSize',20);

        % Plot the points in agreement
        plot(points(agreement_indicies,1),points(agreement_indicies,2),'r.','MarkerSize',20);

        % Plot the test points in blue
        plot(test_source_points(:,1),test_source_points(:,2),'b.','MarkerSize',20);

        % Use the same axis each time, to make plots stay in 'same' place
        axis(debug_axis_limits);
    end

end


%% Step 2: take the top agreements, accumulate them in prep for domains

% Initialize outputs as cell arrays
best_fit_source_indicies          = {};
best_fit_agreement_indicies       = {};
best_fit_flag_is_a_circle         = {};
best_fit_start_angle_in_radians   = {};
best_fit_end_angle_in_radians     = {};


% Find best agreements
flag_find_only_best_agreement = 0;

if flag_find_only_best_agreement ==1
    %% CODE THIS LATER
    % [~, best_agreement_index] = max(agreements);
    % base_point_index = combos_paired(best_agreement_index,1);
    % 
    % % Find the indicies corresponding to the best transverse agreement and station agreement.
    % % For debugging:  figure(2345); clf;
    % best_fit_agreement_indicies{1} = fcn_geometry_findAgreementsOfPointsToLineVector( points, unit_projection_vectors(best_agreement_index,:), base_point_index, transverse_tolerance, station_tolerance, -1);
    % 
    % % Save results to outputs
    % % best_fit_parameters_phi_rho_form{1} = fitted_parameters(best_agreement_index,:);
    % best_fit_source_indicies{1} = combos_paired(best_agreement_index,:);
    % % best_fit_parameters_unitvector_basepoint_distance_form{1} = [unit_projection_vectors(best_agreement_index,:) points(base_point_index,:) station_distances];

else % Find many agreements, up to the number of specified votes
    N_fits = 0;
    remaining_agreements = agreements;

    [best_agreement_count, best_agreement_index] = max(remaining_agreements);

    while best_agreement_count >= points_required_for_agreement

        % Which point index is the base point?
        base_point_index = combos_paired(best_agreement_index,1);

        % Find the indicies in transverse agreement and station agreement
        % with this best fit
        circleCenter  = fitted_parameters(best_agreement_index,1:2);
        circleRadius  = fitted_parameters(best_agreement_index,3);

        [current_agreement_indicies, flag_is_a_circle, start_angle_in_radians, end_angle_in_radians] = ...
            fcn_geometry_findAgreementsOfPointsToArc(points, base_point_index, circleCenter, circleRadius, transverse_tolerance, (station_tolerance), (flag_force_circle_fit),(points_required_for_agreement), (-1));

        % Accumulate results from best fits to prep for domain formation
        N_fits = N_fits+1;
        best_fit_agreement_indicies{N_fits}     = current_agreement_indicies; %#ok<AGROW>
        best_fit_source_indicies{N_fits}        = combos_paired(best_agreement_index,:); %#ok<AGROW>
        best_fit_flag_is_a_circle{N_fits}       = flag_is_a_circle; %#ok<AGROW>
        best_fit_start_angle_in_radians{N_fits} = start_angle_in_radians;  %#ok<AGROW>
        best_fit_end_angle_in_radians{N_fits}   = end_angle_in_radians;  %#ok<AGROW>

        % Remove all combos_paired sums that include this agreement
        for ith_agreement = 1:length(current_agreement_indicies)
            % Grab current agreement
            agreement_to_remove = current_agreement_indicies(ith_agreement);

            % Find which combos contain this agreement
            combos_with_either_index_matching = any(combos_paired==agreement_to_remove,2);

            % Set these combo sums to zero
            remaining_agreements(combos_with_either_index_matching) = 0;
        end

        % Recalculate best fit
        [best_agreement_count, best_agreement_index] = max(remaining_agreements);

    end

end

%% Step 3: form the domains

% Convert outputs to domain type
domains = fcn_INTERNAL_filldomains(points, circleCenter, circleRadius, best_fit_agreement_indicies, best_fit_source_indicies, best_fit_flag_is_a_circle, best_fit_start_angle_in_radians, best_fit_end_angle_in_radians, transverse_tolerance, station_tolerance);


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


    % Produce the sorted list, to create the Hough plot
    [~,sorted_indicies] = sort(agreements,'ascend');

    %% Plot the results in point space

    hold on;
    grid on;
    title('Points and maximum-vote fit, plotted in point-space');
    xlabel('X [meters]');
    ylabel('Y [meters]')

    % Plot the input points
    plot(points(:,1),points(:,2),'k.','MarkerSize',20);

    % Plot the domains
    fcn_geometry_plotFitDomains(domains, fig_num);


    % % Plot the circle fit
    % if best_fit_is_a_circle==1
    %     title('Circle fit');
    % else
    %     title('Arc fit');
    % end
    % 
    % fcn_geometry_plotCircle(best_fitted_parameters(1:2), best_fitted_parameters(3), 'b-',fig_num)
    % fcn_geometry_plotCircle(best_fitted_parameters(1:2), best_fitted_parameters(3)-transverse_tolerance, 'r-',fig_num) 
    % fcn_geometry_plotCircle(best_fitted_parameters(1:2), best_fitted_parameters(3)+transverse_tolerance, 'r-',fig_num) 
    % plot(best_fitted_parameters(1),best_fitted_parameters(2),'b+','MarkerSize',15);
    % 
    % % Plot the best-fit points
    % plot(points(best_agreement_indicies,1),points(best_agreement_indicies,2),'r.','MarkerSize',15);
    % 
    % % Label the points
    % text(points(best_agreement_indicies(end),1),points(best_agreement_indicies(end),2),sprintf('Start angle: %.3f deg', best_start_angle_in_radians*180/pi));
    % text(points(best_agreement_indicies(1),1),points(best_agreement_indicies(1),2),sprintf('End angle: %.3f deg', best_end_angle_in_radians*180/pi));
    % 
    % % Plot the source points
    % plot(points(best_fit_source_indicies,1),points(best_fit_source_indicies,2),'bo','MarkerSize',15);

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end

    %% Plot the Hough space results, from least to best
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


    % TEMPLATE: fitted_parameters(ith_combo,:) = [circleCenter, circleRadius];
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

    % % Plot the best fits
    % best_rho = atan2(best_fitted_parameters(:,2),best_fitted_parameters(:,1));
    % best_distance_circle_center_from_origin = sum((best_fitted_parameters(:,1).^2+best_fitted_parameters(:,2).^2),2).^0.5;
    % best_curvature = 1./best_fitted_parameters(:,3);

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
%% fcn_INTERNAL_filldomains
function domains = fcn_INTERNAL_filldomains(points, circleCenter, circleRadius, best_fit_agreement_indicies, fit_source_indicies, best_fit_flag_is_a_circle, best_fit_start_angle_in_radians, best_fit_end_angle_in_radians, transverse_tolerance, station_tolerance)

% How many domains are there?
N_domains = length(best_fit_agreement_indicies);

% Preallocate the domain structure, to include the unfitted points
domains{N_domains+1} = struct;

unfitted_indicies = ones(length(points(:,1)),1);

% Loop through domains
for ith_domain = 1:N_domains
    
    flag_this_is_a_circle = 0;
    if isempty(station_tolerance) || best_fit_flag_is_a_circle{ith_domain}==1
        flag_this_is_a_circle = 1;
    end

    % Find the points in the domain
    points_in_domain = points(best_fit_agreement_indicies{ith_domain},:);

    % Find the indicies of these points, in the domain. Note that the
    % indicies of the fit are in the original points. Once the points are
    % saved that only belong to this domain, the indicies change. So we
    % need to find what they change to.
    fit_source_index_list = fit_source_indicies{ith_domain};
    domain_specific_indicies = 0*fit_source_index_list; % Initialization
    for ith_index = 1:length(fit_source_index_list(1,:))
        external_index = fit_source_indicies{ith_domain}(ith_index);
        domain_index = find(best_fit_agreement_indicies{ith_domain}==external_index,1);
        domain_specific_indicies(ith_index) = domain_index;
    end
    best_fit_source_indicies    = domain_specific_indicies;


    % Calculate the best-fit domain box
    URHERE
    domain_box = fcn_geometry_domainBoxByType(type_of_domain, varargin)   

    inner_radius = circleRadius
    if flag_this_is_a_circle
        % Find the points that start and end the domain Hough fit
        domain_specific_start_point = points_in_domain(domain_specific_start_index,:);
        domain_specific_end_point = points_in_domain(domain_specific_end_index,:);
        unit_projection_vector = fcn_geometry_calcUnitVector(domain_specific_end_point-domain_specific_start_point);
        unit_orthogonal_vector = unit_projection_vector*[0 1; -1 0];

        base_projection_vectors = points_in_domain - domain_specific_start_point;

        station_distances_of_points_in_domain = ...
            sum(base_projection_vectors.* unit_projection_vector,2);
        max_station = max(station_distances_of_points_in_domain);
        min_station = min(station_distances_of_points_in_domain);
        Hough_fit_start_point = domain_specific_start_point + unit_projection_vector*min_station;
        Hough_fit_end_point   = domain_specific_start_point + unit_projection_vector*max_station;

        best_fit_parameters = [Hough_fit_start_point Hough_fit_end_point];

        boundary_left  = [Hough_fit_start_point; Hough_fit_end_point] + ones(2,1)*unit_orthogonal_vector*transverse_tolerance;
        boundary_right = [Hough_fit_start_point; Hough_fit_end_point] - ones(2,1)*unit_orthogonal_vector*transverse_tolerance;
    else


    end

    centers = [3 4];
    radii = 2;
    start_angle_in_radians = 45 * pi/180;
    end_angle_in_radians = 135 * pi/180;
    degree_step = []; % Default is 1 degree
    format = [];
    fig_num = [];

    arc_points_matrix = fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (degree_step), (format),(fig_num));


    domain_box = [boundary_left; flipud(boundary_right)];
    domainShape = polyshape(domain_box(:,1),domain_box(:,2),'Simplify',false,'KeepCollinearPoints',true);

    % Which points remain unfitted?
    unfitted_indicies(best_fit_agreement_indicies{ith_domain})=0;

    %% Save results into the domain structure
    domains{ith_domain} = fcn_geometry_fillEmptyDomainStructure;

    if flag_this_is_a_circle
        domains{ith_domain}.best_fit_type = 'Hough circle';
    else
        domains{ith_domain}.best_fit_type = 'Hough arc';
    end
    domains{ith_domain}.points_in_domain         = points_in_domain;
    domains{ith_domain}.best_fit_source_indicies = best_fit_source_indicies;
    domains{ith_domain}.best_fit_domain_box      = domainShape;
    domains{ith_domain}.best_fit_parameters      = best_fit_parameters;
end


% Save unfitted points into last structure
domains{end} = fcn_geometry_fillEmptyDomainStructure;
domains{end}.best_fit_type = 'unfitted';
domains{end}.points_in_domain = points(unfitted_indicies==1,:);

end % Ends fcn_INTERNAL_filldomains
