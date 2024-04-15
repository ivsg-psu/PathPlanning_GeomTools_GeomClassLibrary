function domains = fcn_geometry_fitHoughLine(points, transverse_tolerance, station_tolerance, varargin)
% fcn_geometry_fitHoughLine
% Checks all permutations between points to fit a line (N choose 2), then
% calculates the line fit in polar form (rho and phi), and determines which
% of the points are within a tolerance distance of that line fit. This is
% calculated for all possible 2-point permutations, ordered in N-choose-2
% format. The best-fit parameters are returned
%
% Format:
% [best_fitted_parameters, agreement_indicies] = fcn_geometry_fitHoughLine(points,varargin)
%
% INPUTS:
%      points: a Nx2 vector where N is the number of points, but at least 2 rows.
%
%      transverse_tolerance: the orthogonal distance between the points and
%      the linear curve fit that indicate whether a point "belongs" to the
%      fit. A point belongs to the fit if the transverse distance is less
%      than or equal to the transverse_tolerance
%
%      station_tolerance: the projection distance between the points in a
%      curve fit, along the direction of the line, that indicate whether a
%      point "belongs" to the linear fit of a line segment or not. A line
%      segment is considered continous if station interval distance is less
%      than or equal to the station_tolerance. Set station_tolerance to an
%      empty value, [], to avoid line segment fitting and instead find pure
%      line fits.
%
%      (OPTIONAL INPUTS)
%
%      points_required_for_agreement: the number of points required for an
%      agreement to be valid, with minimum value of 3. If left empty, then
%      the best agreement will always be returned. If a value is given,
%      line fitting will continue by clustering data until there are no
%      fits greater than or equal to points_required_for_agreement. The
%      results of the line fit will be saved in a cell array.
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
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_calcUnitVector
%      fcn_geometry_findAgreementsOfPointsToLineVector    
%      fcn_geometry_fillEmptyDomainStructure
%      fcn_geometry_plotFitDomains
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_fitHoughLine
% for a full test suite.

% This function was written on 2023_12_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2023_12_14 - S. Brennan
% -- wrote the code
% 2023_12_29 - S. Brennan
% -- added debug mode that uses environmental variables
% -- added speed mode
% 2023_12_31 - S. Brennan
% -- externalized station checks to fcn_geometry_findPointsInSequence
% -- now allows both line and line segment fitting
% -- modified to only ouput the best fit
% 2024_01_14 - S. Brennan
% -- added points_required_for_agreement option
% -- fixed minor bug where fastmode looking at varargin{1}, not
% varargin{end}
% -- added best_fitted_parameters_unitvector_basepoint_distance_form output
% -- added multi-line fitting option, much faster than doing one-at-a-time
% 2024_01_15 - S. Brennan
% -- functionalized fcn_INTERNAL_findPhisRhos
% -- changed outputs to domain types
% -- added fcn_geometry_plotFitDomains for plotting

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==5 && isequal(varargin{end},-1))
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
        narginchk(3,5);

        % Check the points input to be length greater than or equal to 2
        fcn_DebugTools_checkInputsToFunctions(...
            points, '2column_of_numbers',[2 3]);

        % Check the transverse_tolerance input is a positive single number
        fcn_DebugTools_checkInputsToFunctions(transverse_tolerance, 'positive_1column_of_numbers',1);

        % Check the station_tolerance input is a positive single number
        if ~isempty(station_tolerance)
            fcn_DebugTools_checkInputsToFunctions(station_tolerance, 'positive_1column_of_numbers',1);
        end
    end
end

% Does user specify points_required_for_agreement?
flag_find_only_best_agreement = 1;
if 4<= nargin
    temp = varargin{1};
    if ~isempty(temp)
        points_required_for_agreement = temp;
        flag_find_only_best_agreement = 0;
        if points_required_for_agreement<3
            error('The input points_required_for_agreement must be greater than or equal to 3.')
        end
    end
end


% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if 5<= nargin && 0==flag_max_speed
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end


%% Solve for the line fit
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
N_points = length(points(:,1));

% Find all possible 2-point permutations
combos_paired = nchoosek(1:N_points,2);

% Find all projections from the paired combos to each other
projection_vectors = points(combos_paired(:,2),:) - points(combos_paired(:,1),:);

% Convert the projection vectors into unit vectors in projection and
% orthgonal directions
unit_projection_vectors = fcn_geometry_calcUnitVector(projection_vectors);

% Find phi and rho parameters corresponding to each combo and vector
[phis, rhos] = fcn_INTERNAL_findPhisRhos(points, unit_projection_vectors, combos_paired);

% Save the results for output from the function and for plotting
% fitted_parameters = [phis rhos];

% Find the agreements between points and every single line fit
agreements = ...
    fcn_INTERNAL_findAgreementsOfPointsToFits(...
    points, unit_projection_vectors, combos_paired, transverse_tolerance, station_tolerance, flag_max_speed);


% Initialize outputs as cell arrays
% best_fit_parameters_phi_rho_form    = {};
best_fit_parameters_unitvector_basepoint_distance_form = {};
% best_fit_source_indicies  = {};
best_fit_agreement_indicies   = {};
max_agreement_indicies = {}; % For plotting

% Find best agreement
if flag_find_only_best_agreement ==1
    [~, max_agreement_index] = max(agreements);
    base_point_index = combos_paired(max_agreement_index,1);

    % Find the indicies corresponding to the best transverse agreement and station agreement.
    % For debugging:  figure(2345); clf;
    [best_fit_agreement_indicies{1}, station_distances] = fcn_geometry_findAgreementsOfPointsToLineVector( points, unit_projection_vectors(max_agreement_index,:), base_point_index, transverse_tolerance, station_tolerance, -1);

    % Save results to outputs
    % best_fit_parameters_phi_rho_form{1} = fitted_parameters(best_agreement_index,:);
    % best_fit_source_indicies{1} = combos_paired(max_agreement_index,:);   
    best_fit_parameters_unitvector_basepoint_distance_form{1} = [unit_projection_vectors(max_agreement_index,:) points(base_point_index,:) station_distances];

    max_agreement_indicies{1} = max_agreement_index;
else % Find many agreements, up to the number of specified votes
    N_fits = 0;
    remaining_agreements = agreements;

    [best_agreement_count, max_agreement_index] = max(remaining_agreements);


    while best_agreement_count >= points_required_for_agreement

        % Which point index is the base point?
        base_point_index = combos_paired(max_agreement_index,1);

        % Find the indicies in transverse agreement and station agreement.
        [current_agreement_indicies, current_station_distances] = ...
            fcn_geometry_findAgreementsOfPointsToLineVector( points, unit_projection_vectors(max_agreement_index,:), base_point_index, transverse_tolerance, station_tolerance, -1);

        % Save results from prior best fit to outputs
        N_fits = N_fits+1;
        max_agreement_indicies{N_fits} = max_agreement_index; %#ok<AGROW>
        best_fit_agreement_indicies{N_fits} = current_agreement_indicies; %#ok<AGROW>
        % best_fit_source_indicies{N_fits} = combos_paired(max_agreement_index,:); %#ok<AGROW>
        % best_fit_parameters_phi_rho_form{N_fits} = fitted_parameters(best_agreement_index,:); %#ok<AGROW>
        best_fit_parameters_unitvector_basepoint_distance_form{N_fits} = [unit_projection_vectors(max_agreement_index,:) points(base_point_index,:) current_station_distances]; %#ok<AGROW>
        % best_fit_source_indicies{N_fits} = combos_paired(best_agreement_index,:); %#ok<AGROW>

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
        [best_agreement_count, max_agreement_index] = max(remaining_agreements);

    end
    
end


% Convert outputs to domain type
domains = fcn_INTERNAL_filldomains(points, best_fit_parameters_unitvector_basepoint_distance_form, best_fit_agreement_indicies, transverse_tolerance, station_tolerance);


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

    % Produce the sorted list, to create the Hough plot
    [~,sorted_indicies] = sort(agreements,'ascend');

    % Save the number of permutations
    N_permutations = length(combos_paired(:,1));

    %% Plot the results in point space
    hold on;
    grid on;
    axis equal;
    title('Points and maximum-vote fit, plotted in point-space');
    xlabel('X [meters]');
    ylabel('Y [meters]')

    % Plot the input points
    plot(points(:,1),points(:,2),'k.','MarkerSize',40);

    % Plot the domains
    fcn_geometry_plotFitDomains(domains, fig_num);

    % Plot all line fits? NOTE: this can be confusing as these are the
    % INPUT line segments, not the actual line segments at the end of the
    % fit. NOTE 2: this code was NOT updated and needs to be edited for
    % cell array changes that occurred after it was written.
    if 1==0
        start_points = points(combos_paired(:,1),:);
        start_angles = atan2(start_points(:,2),start_points(:,1));
        start_radii  = rhos ./ cos(start_angles - phis);
        start_points_calculated = start_radii.*[cos(start_angles) sin(start_angles)];

        end_points = points(combos_paired(:,2),:);
        end_angles = atan2(end_points(:,2),end_points(:,1));
        end_radii  = rhos ./ cos(end_angles - phis);
        end_points_calculated = end_radii.*[cos(end_angles) sin(end_angles)];

        N_steps = 10;
        for ith_best_line = min(N_steps,length(agreements)):-1:1
            ith_line = sorted_indicies(end-ith_best_line+1);
            plot_color = (ith_best_line - 1)/N_steps*[1 1 1];
            line_size = ceil((N_steps - ith_best_line + 1)/N_steps * 5);

            plot(...
                [start_points_calculated(ith_line,1) end_points_calculated(ith_line,1)],...
                [start_points_calculated(ith_line,2) end_points_calculated(ith_line,2)], '-','LineWidth',line_size,'Color',plot_color);
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

    %% Plot the Hough space results, from least to best
    figure(fig_num+1);
    clf;
    hold on;
    grid on;

    title('Hough space plot of all fits');
    subtitle('Points are color-weighted wherein higher darkness indicates higher votes');

    % Plot the Hough fits in batches, making the points bigger and darker
    % in each batch to emphasize which points were the best fits.
    N_steps = 20;
    indicies = linspace(1,N_permutations,N_steps);
    for ith_plot = 1:N_steps-1
        plot_indicies_start = ceil(indicies(ith_plot));
        plot_indicies_end = floor(indicies(ith_plot+1));

        sorted_indicies_to_plot = sorted_indicies(plot_indicies_start:plot_indicies_end);
        plot_color = (N_steps - ith_plot + 1)/N_steps*[1 1 1];
        marker_size = ceil(ith_plot*20/N_steps);
        plot(phis(sorted_indicies_to_plot),rhos(sorted_indicies_to_plot),'.','MarkerSize',marker_size,'Color',plot_color);
    end

    % Plot the best-fit phis and rhos with a red dot
    if flag_find_only_best_agreement ==1
        index_to_plot = max_agreement_indicies{1};
        fitted_parameters = [phis(index_to_plot), rhos(index_to_plot)];
        plot(fitted_parameters(1,1),fitted_parameters(1,2),'r.','MarkerSize',30);
    else
        for ith_agreement = 1:length(domains)-1
            % Get current color
            current_color = fcn_geometry_fillColorFromNumberOrName(ith_agreement);

            % Get points to plot
            index_to_plot = max_agreement_indicies{ith_agreement};
            fitted_parameters = [phis(index_to_plot), rhos(index_to_plot)];


            % Update plot
            plot(fitted_parameters(1,1),fitted_parameters(1,2),'.','MarkerSize',50,'Color',current_color);
        end
    end


    xlabel('Phi [radians]');
    ylabel('Rho [meters]');

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

%% fcn_INTERNAL_findPhisRhos
function [phis, rhos] = fcn_INTERNAL_findPhisRhos(points, unit_projection_vectors, combos_paired)
unit_orthogonal_vectors = unit_projection_vectors*[0 1; -1 0];

% Calculate the angles of the unit projection vectors
angles_relative_to_x_axis = atan2(unit_projection_vectors(:,2),unit_projection_vectors(:,1));

% For polar line definitions, the phis are the angle of the vector's
% orthogonal, and so a vector pointed straight up will have an phi angle of
% zero. We can convert by subtracting 90 degrees from the vector angle to
% calculate the phi angles.
phis = angles_relative_to_x_axis - pi/2;

% The radius of the line relative to the origin is just the projection
% distance from the origin to the points on the line, dot producted with
% the negative of the orthogonal vectors
mixed_rhos = -1*sum(points(combos_paired(:,1),:).*unit_orthogonal_vectors,2);

% Sometimes the line has "negative" radius, and if so, we need to flip the
% angle measured and as well make these radii positive
rhos = mixed_rhos;
rhos(mixed_rhos<0) = rhos(mixed_rhos<0)*-1;
phis(mixed_rhos<0) = phis(mixed_rhos<0) + pi;

% Make sure phis wrap around correctly
phis = mod(phis,2*pi);
end % Ends fcn_INTERNAL_findPhisRhos

%% fcn_INTERNAL_filldomains
function domains = fcn_INTERNAL_filldomains(points, best_fit_parameters_unitvector_basepoint_distance_form, best_fit_agreement_indicies, transverse_tolerance, station_tolerance)

% How many domains are there?
N_domains = length(best_fit_agreement_indicies);

% Preallocate the domain structure, to include the unfitted points
domains{N_domains+1} = struct;

unfitted_indicies = ones(length(points(:,1)),1);

% Loop through domains
for ith_domain = 1:N_domains

    % Is this a line or line segment?
    flag_this_is_a_line = 0;
    if isempty(station_tolerance)
        flag_this_is_a_line = 1;
    end

    % Grab the points in the domain
    points_in_domain = points(best_fit_agreement_indicies{ith_domain},:);

    % Calculate the best-fit domain box
    % FORMAT:     best_fit_parameters_unitvector_basepoint_distance_form{1} = [unit_projection_vectors(max_agreement_index,:) points(base_point_index,:) station_distances];

    % Find the points that start and end the domain Hough fit
    fit_parameters              = best_fit_parameters_unitvector_basepoint_distance_form{ith_domain};
    unit_projection_vector      = fit_parameters(1,1:2);
    domain_specific_start_point = fit_parameters(1,3:4);

    % Calculate a box that goes all around the points in the domain
    base_projection_vectors = points_in_domain - domain_specific_start_point;

    % Find the max and min stations among the points
    station_distances_of_points_in_domain = ...
        sum(base_projection_vectors.* unit_projection_vector,2);
    max_station = max(station_distances_of_points_in_domain);
    min_station = min(station_distances_of_points_in_domain);
    
    Hough_fit_start_point = domain_specific_start_point + unit_projection_vector*min_station;
    Hough_fit_end_point   = domain_specific_start_point + unit_projection_vector*max_station;

    unit_orthogonal_vector = unit_projection_vector*[0 1; -1 0];
    boundary_left  = [Hough_fit_start_point; Hough_fit_end_point] + ones(2,1)*unit_orthogonal_vector*transverse_tolerance;
    boundary_right = [Hough_fit_start_point; Hough_fit_end_point] - ones(2,1)*unit_orthogonal_vector*transverse_tolerance;

    domain_box = [boundary_left; flipud(boundary_right)];
    domainShape = polyshape(domain_box(:,1),domain_box(:,2),'Simplify',false,'KeepCollinearPoints',true);

    % Determine the best_fit_parameters
    if flag_this_is_a_line
        best_fit_parameters = fit_parameters(1,1:4);
    else
        best_fit_parameters = fit_parameters(1,1:6);
    end


    % Update which points remain unfitted
    unfitted_indicies(best_fit_agreement_indicies{ith_domain})=0;

    % Save results into the domain structure
    domains{ith_domain} = fcn_geometry_fillEmptyDomainStructure;

    if flag_this_is_a_line
        domains{ith_domain}.best_fit_type = 'Hough line';
    else
        domains{ith_domain}.best_fit_type = 'Hough segment';
    end
    domains{ith_domain}.points_in_domain         = points_in_domain;
    domains{ith_domain}.best_fit_domain_box      = domainShape;
    domains{ith_domain}.best_fit_parameters      = best_fit_parameters;
end


% Save unfitted points into last structure
domains{end} = fcn_geometry_fillEmptyDomainStructure;
domains{end}.best_fit_type = 'unfitted';
domains{end}.points_in_domain = points(unfitted_indicies==1,:);

end % Ends fcn_INTERNAL_filldomains

%% fcn_INTERNAL_findAgreementsOfPointsToFits
function agreement_counts = ...
    fcn_INTERNAL_findAgreementsOfPointsToFits(...
    input_points, ...
    unit_projection_vectors, ...
    combos_paired, ...
    transverse_tolerance, station_tolerance, flag_max_speed)

N_combos = length(unit_projection_vectors(:,1));

% Initialize output vector
agreement_counts = zeros(N_combos,1);

% Loop through all the combos, recording agreement with each combo
if 0==flag_max_speed
    h_waitbar = waitbar(0,'Calculating line and line segment fits...');
end

for ith_vector = 1:N_combos
    % fprintf(1,'Checking %.0d of %.0d\n',ith_vector,N_combos);
    if 0==flag_max_speed
        if 0==mod(ith_vector,1000)
            waitbar(ith_vector/N_combos);
        end
    end

    base_point_index = combos_paired(ith_vector,1);

    % Find the indicies in transverse agreement and station agreement.
    % For debugging:  figure(2345); clf;
    indicies_in_both_lateral_and_station_agreement = fcn_geometry_findAgreementsOfPointsToLineVector( input_points, unit_projection_vectors(ith_vector,:), base_point_index, transverse_tolerance, station_tolerance, -1);

    agreement_count = length(indicies_in_both_lateral_and_station_agreement);
    agreement_counts(ith_vector) = agreement_count;
end
if 0==flag_max_speed
    close(h_waitbar)
end
end % Ends fcn_INTERNAL_findAgreementsOfPointsToFits


