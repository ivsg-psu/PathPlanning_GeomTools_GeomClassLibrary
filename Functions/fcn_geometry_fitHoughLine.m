function [best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies] = ...
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
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      best_fitted_parameters: the angle, phi, and radius, rho, for the
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

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{1},-1))
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
        narginchk(3,4);

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

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if 4<= nargin && 0==flag_max_speed
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
[best_agreement_index, best_agreement_indicies, agreements] = ...
    fcn_INTERNAL_findAgreementsOfPointsToFits(...
    points, unit_projection_vectors, combos_paired, transverse_tolerance, station_tolerance);


% Save results to outputs
best_fitted_parameters = fitted_parameters(best_agreement_index,:);
best_fit_source_indicies = combos_paired(best_agreement_index,:);

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
    plot(points(best_agreement_indicies,1),points(best_agreement_indicies,2),'r.','MarkerSize',15); 


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
    plot(best_fitted_parameters(1,1),best_fitted_parameters(1,2),'r.','MarkerSize',30);

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
function [best_agreement_index, best_agreement_indicies, agreements] = ...
    fcn_INTERNAL_findAgreementsOfPointsToFits(...
    input_points, ...
    unit_projection_vectors, ...
    combos_paired, ...
    transverse_tolerance, station_tolerance)

flag_do_debug = 0;

N_combos = length(unit_projection_vectors(:,1));
N_points = length(input_points(:,1));

% Initialize output vectors
best_agreement_indicies_binary_form = zeros(1,N_points);
best_agreement_count = -inf;
best_agreement_index = [];
agreements = zeros(N_combos,1);

% Loop through all the combos, trying to find best agreement
for ith_vector = 1:N_combos
    base_point_index = combos_paired(ith_vector,1);
    next_point_index = combos_paired(ith_vector,2);
    base_point = input_points(base_point_index,:);
    next_point = input_points(next_point_index,:);

    % Find the indicies in transverse agreement. This is done by finding
    % the lateral distances of the points to the fit line, and then keeping
    % only the points within the lateral distance.
    unit_orthogonal_vector = unit_projection_vectors(ith_vector,:)*[0 1; -1 0];
    base_projection_vectors = input_points - base_point;
    lateral_distances = sum(unit_orthogonal_vector.*base_projection_vectors,2);

    indicies_in_lateral_agreement = find((abs(lateral_distances)<=transverse_tolerance)');
    base_point_index_in_lateral_agreement = find(indicies_in_lateral_agreement == base_point_index);

    % Next find the indicies in station agreement. To do this, find the
    % station coordinates for every point as predicted by their tangent
    % distances. Then sort the points by this distance, keeping track of
    % how indicies move when sorted. Then check differences in station
    % change from point to point, cutting off the fit if and when these
    % differences are larger than the station tolerance.
 
    % Find station distances as projections of tangent vectors    
    tangent_distances_of_points_in_lateral_agreement = ...
        sum(base_projection_vectors(indicies_in_lateral_agreement,:).* unit_projection_vectors(ith_vector,:),2);

    if ~isempty(station_tolerance)
        % Sort the station distances and find those in agreement with station
        % tolerance
        indicies_in_station_agreement = ...
            fcn_geometry_findPointsInSequence(...
            tangent_distances_of_points_in_lateral_agreement, ...
            base_point_index_in_lateral_agreement, ...
            station_tolerance, -1);

        % Find indicies in both lateral and station agreement
        indicies_in_both_lateral_and_station_agreement = indicies_in_lateral_agreement(indicies_in_station_agreement);
    else
        indicies_in_both_lateral_and_station_agreement = indicies_in_lateral_agreement;
    end

    % Check the results? (for debugging)    
    if flag_do_debug
        figure(4747);
        clf;
        hold on;
        grid on;
        axis equal;
        
        plot(base_point(:,1),base_point(:,2),'g.','MarkerSize',40);
        plot(next_point(:,1),next_point(:,2),'r.','MarkerSize',40);
        
        plot(input_points(:,1),input_points(:,2),'k.','MarkerSize',30)
        plot(input_points(indicies_in_lateral_agreement,1),input_points(indicies_in_lateral_agreement,2),'m.','MarkerSize',30)
        
        plot(input_points(indicies_in_both_lateral_and_station_agreement,1),input_points(indicies_in_both_lateral_and_station_agreement,2),'c.-','MarkerSize',20,'LineWidth',3)
    end

    agreement_count = length(indicies_in_both_lateral_and_station_agreement);
    agreements(ith_vector) = agreement_count;

    if agreement_count>best_agreement_count
        best_agreement_index = ith_vector;

        % Save new "best" agreement total count
        best_agreement_count = agreement_count;

        % Save new "best" indicies
        best_agreement_indicies_binary_form = zeros(1,N_points);
        best_agreement_indicies_binary_form(1,indicies_in_both_lateral_and_station_agreement) = 1;
    end

end

best_agreement_indicies = find(best_agreement_indicies_binary_form==1);
end % Ends fcn_INTERNAL_findAgreementsOfPointsToFits

