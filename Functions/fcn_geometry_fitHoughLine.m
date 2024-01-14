function [best_fitted_parameters_phi_rho_form, best_fit_source_indicies, best_agreement_indicies, best_fitted_parameters_unitvector_basepoint_distance_form] = ...
    fcn_geometry_fitHoughLine(points, transverse_tolerance, station_tolerance, varargin)
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
%      best_fitted_parameters_phi_rho_form: the Hough space results are
%      shown in polar coordinates as this is the minimum parameter form,
%      requiring only 2 parameters. This is described by the angle, phi, of
%      the segment relative to the x-axis and radius, rho, of the nearest
%      approach of the line projection to the origin, calculated for the
%      line or line segment of the best point-to-point permutation. This is
%      returned as a [1 x 2] matrix ordered as [phi rho]. The line fit is
%      of the form: r = rho / [cos(theta - phi)]
%
%      best_fit_source_indicies: the two indicies, in [1x2] vector format,
%      of the points that produced the best fit.
%
%      best_agreement_indicies: the indicies of the points that are within
%      agreement of the best-fit parameters, given the transverse and
%      station tolerance settings.
%
%      best_fitted_parameters_unitvector_basepoint_distance_form: the Hough
%      voting of the fits are conducted in a form that uses a unit vector,
%      projected from a base point, for a distance range that can be both
%      positive and negative. This result, for the best fits, is returned
%      as a [1 x 6] vector composed as follows:
%      [unitvector_x unitvector_y basepoint_x basepoint_y distance_min
%      distance_max]
%
% DEPENDENCIES:
%      nchoosek
%      fcn_geometry_findPointsInSequence
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

% Save the results for output from the function and for plotting
fitted_parameters = [phis rhos];

% Find the agreements between points and every single line fit
agreements = ...
    fcn_INTERNAL_findAgreementsOfPointsToFits(...
    points, unit_projection_vectors, combos_paired, transverse_tolerance, station_tolerance);

% Find best agreement
if flag_find_only_best_agreement ==1
    [~, best_agreement_index] = max(agreements);
    base_point_index = combos_paired(best_agreement_index,1);

    % Find the indicies in transverse agreement and station agreement.
    % For debugging:  figure(2345); clf;
    best_agreement_indicies = fcn_geometry_findAgreementsOfPointsToLineVector( points, unit_projection_vectors(best_agreement_index,:), base_point_index, transverse_tolerance, station_tolerance, -1);

    % Save results to outputs
    best_fitted_parameters_phi_rho_form = fitted_parameters(best_agreement_index,:);
    best_fit_source_indicies = combos_paired(best_agreement_index,:);

    URHERE
    best_fitted_parameters_unitvector_basepoint_distance_form = [];

else
    % Initialize outputs as a cell array
    best_fitted_parameters_phi_rho_form    = {};
    best_fit_source_indicies  = {};    
    best_agreement_indicies   = {};
    N_fits = 0;

    remaining_agreements = agreements;

    [best_agreement_count, best_agreement_index] = max(remaining_agreements);

    while best_agreement_count >= points_required_for_agreement

        base_point_index = combos_paired(best_agreement_index,1);

        % Find the indicies in transverse agreement and station agreement.
        % For debugging:  figure(2345); clf;
        current_agreement_indicies = ...
            fcn_geometry_findAgreementsOfPointsToLineVector( points, unit_projection_vectors(best_agreement_index,:), base_point_index, transverse_tolerance, station_tolerance, -1);

        % Save results to outputs
        N_fits = N_fits+1;
        best_agreement_indicies{N_fits} = current_agreement_indicies; %#ok<AGROW>
        best_fitted_parameters_phi_rho_form{N_fits} = fitted_parameters(best_agreement_index,:); %#ok<AGROW>
        best_fit_source_indicies{N_fits} = combos_paired(best_agreement_index,:); %#ok<AGROW>

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
    title('Points and maximum-vote fit, plotted in point-space');
    xlabel('X [meters]');
    ylabel('Y [meters]')

    % Produce the sorted list, to create the Hough plot
    [~,sorted_indicies] = sort(agreements,'ascend');

    % Save the number of permutations
    N_permutations = length(combos_paired(:,1));

    % Plot the results in point space

    % Plot the input points
    plot(points(:,1),points(:,2),'k.','MarkerSize',20);

    % Plot the best-fit points
    if flag_find_only_best_agreement ==1
        plot(points(best_agreement_indicies,1),points(best_agreement_indicies,2),'r.','MarkerSize',15);
    else
        for ith_agreement = 1:length(best_agreement_indicies)
            fitted_agreement_indicies = best_agreement_indicies{ith_agreement};
            plot(points(fitted_agreement_indicies,1),points(fitted_agreement_indicies,2),'o','MarkerSize',(15+5*ith_agreement));
        end
    end

    % Plot all line fits NOTE: this can be confusing as these are the INPUT
    % line segments, not the actual line segments at the end of the fit.
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

    % Plot the Hough space results, from least to best
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
            plot(best_fitted_parameters_phi_rho_form(1,1),best_fitted_parameters_phi_rho_form(1,2),'r.','MarkerSize',30);
    else
        for ith_agreement = 1:length(best_agreement_indicies)
            fitted_parameters = best_fitted_parameters_phi_rho_form{ith_agreement};
            plot(fitted_parameters(1,1),fitted_parameters(1,2),'.','MarkerSize',30);
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



%% fcn_INTERNAL_findAgreementsOfPointsToFits
function agreement_counts = ...
    fcn_INTERNAL_findAgreementsOfPointsToFits(...
    input_points, ...
    unit_projection_vectors, ...
    combos_paired, ...
    transverse_tolerance, station_tolerance)

N_combos = length(unit_projection_vectors(:,1));

% Initialize output vector
agreement_counts = zeros(N_combos,1);

% Loop through all the combos, recording agreement with each combo
for ith_vector = 1:N_combos
    base_point_index = combos_paired(ith_vector,1);

    % Find the indicies in transverse agreement and station agreement.
    % For debugging:  figure(2345); clf;
    indicies_in_both_lateral_and_station_agreement = fcn_geometry_findAgreementsOfPointsToLineVector( input_points, unit_projection_vectors(ith_vector,:), base_point_index, transverse_tolerance, station_tolerance, -1);

    agreement_count = length(indicies_in_both_lateral_and_station_agreement);
    agreement_counts(ith_vector) = agreement_count;
end

end % Ends fcn_INTERNAL_findAgreementsOfPointsToFits


