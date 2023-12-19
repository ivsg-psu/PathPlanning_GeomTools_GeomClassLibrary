function [fitted_parameters, agreement_indicies] = fcn_geometry_fitHoughLine(points, transverse_tolerance, station_tolerance, varargin)
% fcn_geometry_fitHoughLine
% Finds all permutations between points to fit a line (N choose 2),
% calculates the line fit in polar form (rho and phi), and determines which
% of the points are within a tolerance distance of that line fit. This is
% calculated for all possible 2-point permutations, ordered in N-choose-2
% format.
%
% Format: 
% [fitted_parameters, agreement_indicies] = fcn_geometry_fitHoughLine(points,varargin)
%
% INPUTS:
%      points: a Nx2 vector where N is the number of points, but at least 2 rows. 
%
%      transverse_tolerance: the orthogonal distance between the points and
%      the linear curve fit that indicate whether a point "belongs" to the
%      fit (if distance is less than or equal to the tolerance), or is
%      "outside" the fit (if distance is greater than the tolerance).
%
%      station_tolerance: the projection distance between the points in a
%      curve fit, along the direction of the line, that indicate whether a
%      point "belongs" to the linear fit (if distance is less than or equal
%      to the tolerance), or is "outside" the fit (if distance is greater
%      than the tolerance).
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot the results.
%
% OUTPUTS:
%      Let K represent the number of permutations possible in N-choose-2
%      standard ordering (see function nchoosek). For each K value, this
%      returns:
%
%      fitted_parameters: the angle, phi, and radius, rho, for the line
%      fitting two points in the point permutation. This is returned as a
%      [K x 2] matrix ordered as [phi rho]. The line fit is of the form:
%      r = rho / [cos(theta - phi)]
%
%      agreement_indicies: a row of 1 or 0 for all points, with 1
%      indicating that a given point is "within" the line projection, 0 if
%      not. This is returned as a [K x N] matrix.
%
% DEPENDENCIES:
%      nchoosek
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_fitHoughLine
% for a full test suite.

% This function was written on 2023_12_14 by S. Brennan
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
    narginchk(3,4);
    
    % Check the points input to be length greater than or equal to 2
    fcn_DebugTools_checkInputsToFunctions(...
        points, '2column_of_numbers',[2 3]);
    
    % Check the transverse_tolerance input is a positive single number
    fcn_DebugTools_checkInputsToFunctions(transverse_tolerance, 'positive_1column_of_numbers',1);
    
    % Check the station_tolerance input is a positive single number
    fcn_DebugTools_checkInputsToFunctions(station_tolerance, 'positive_1column_of_numbers',1);

end


% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if 4<= nargin
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

% Save the results for output from the function
fitted_parameters = [phis rhos];

%% Find the agreements between points and every single line fit
[agreement_indicies] = fcn_INTERNAL_findAgreementsOfPointsToFits(points,unit_projection_vectors, combos_paired, transverse_tolerance, station_tolerance);


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


    agreements = sum(agreement_indicies,2);
    [~,sorted_indicies] = sort(agreements,'ascend');
    N_permutations = length(combos_paired(:,1));

    % Plot the results in point space

    % Plot the input points
    plot(points(:,1),points(:,2),'k.','MarkerSize',20);

    % Plot the best-fit points
    agreements = sum(agreement_indicies,2);
    [~,best_agreement_index] = max(agreements);
    best_fit_point_indicies = find(agreement_indicies(best_agreement_index,:));
    plot(points(best_fit_point_indicies,1),points(best_fit_point_indicies,2),'r.','MarkerSize',15);


    % Plot all line fits?
    if 1==0
        start_points = points(combos_paired(:,1),:);
        start_angles = atan2(start_points(:,2),start_points(:,1));
        start_radii  = rhos ./ cos(start_angles - phis);
        start_points_calculated = start_radii.*[cos(start_angles) sin(start_angles)];

        end_points = points(combos_paired(:,2),:);
        end_angles = atan2(end_points(:,2),end_points(:,1));
        end_radii  = rhos ./ cos(end_angles - phis);
        end_points_calculated = end_radii.*[cos(end_angles) sin(end_angles)];



        for best_line = 1:min(1000,length(agreements))
            ith_line = sorted_indicies(best_line);
            plot([start_points_calculated(ith_line,1) end_points_calculated(ith_line,1)],[start_points_calculated(ith_line,2) end_points_calculated(ith_line,2)], '-','MarkerSize',10);
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
function [agreement_indicies] = fcn_INTERNAL_findAgreementsOfPointsToFits(input_points, unit_projection_vectors, combos_paired, transverse_tolerance, station_tolerance)

N_combos = length(unit_projection_vectors(:,1));
N_points = length(input_points(:,1));

agreements = zeros(N_combos,1);
agreement_indicies = zeros(N_combos,N_points);

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

    % Next find the indicies in station agreement. To do this, find the
    % station coordinates for every point as predicted by their tangent
    % distances. Then sort the points by this distance, keeping track of
    % how indicies move when sorted. Then check differences in station
    % change from point to point, cutting off the fit if and when these
    % differences are larger than the station tolerance.
 
    % Find station coordinates as tangent distances    
    tangent_distances_of_points_in_lateral_agreement = sum(base_projection_vectors(indicies_in_lateral_agreement,:).* unit_projection_vectors(ith_vector,:),2);

    % Sort the distances in direction of point-to-point projection
    [sorted_distances,sort_order_of_indicies] = sort(tangent_distances_of_points_in_lateral_agreement);
    sorted_indicies = indicies_in_lateral_agreement(sort_order_of_indicies);


    % Find the starting point in the sort associated with the base point
    % index
    starting_point_sort = find(base_point_index == sorted_indicies); 
    distances_after_base_point = sorted_distances(starting_point_sort:end);
    distances_before_base_point = sorted_distances(1:starting_point_sort);
    diff_tangent_distances_after_base_point = diff(distances_after_base_point);
    diff_tangent_distances_before_base_point = diff(distances_before_base_point);
    
    % Which indicies after the base point are in agreement?
    stop_point_after_base_point = find(diff_tangent_distances_after_base_point>station_tolerance,1);
    if ~isempty(stop_point_after_base_point)
        sorted_indicies_after_base_point = (starting_point_sort:(starting_point_sort+stop_point_after_base_point-1));
    else
        sorted_indicies_after_base_point = (starting_point_sort:length(sorted_indicies));
    end

    % Which indicies before the base point are in agreement?
    stop_point_before_base_point = find(diff_tangent_distances_before_base_point>station_tolerance,1,'last');
    if ~isempty(stop_point_before_base_point)
        sorted_indicies_before_base_point = ((stop_point_before_base_point+1):(starting_point_sort-1));
    else
        sorted_indicies_before_base_point = (1:(starting_point_sort-1));
    end

    % Store before and after indicies together to produce indicies in
    % station agreement
    sorted_indicies_to_keep = [sorted_indicies_before_base_point sorted_indicies_after_base_point];

    % Keep only the lateral agreement indicies that meet the station
    % tolerance. These indicies must therefore meet the lateral and station
    % tolerances.
    indicies_in_both_lateral_and_station_agreement = sorted_indicies(sorted_indicies_to_keep);

    % Check the results? (for debugging)
    if 1==0
        figure(4747);
        clf;
        hold on;
        grid on;
        axis equal;
        
        plot(base_point(:,1),base_point(:,2),'g.','MarkerSize',40);
        plot(next_point(:,1),next_point(:,2),'r.','MarkerSize',40);
        
        plot(input_points(:,1),input_points(:,2),'k.','MarkerSize',30)
        plot(input_points(indicies_in_lateral_agreement,1),input_points(indicies_in_lateral_agreement,2),'m.','MarkerSize',30)
        
        plot(input_points(sorted_indicies,1),input_points(sorted_indicies,2),'c.','MarkerSize',20)
        plot(input_points(sorted_indicies,1),input_points(sorted_indicies,2),'c-')

        plot(input_points(indicies_in_both_lateral_and_station_agreement,1),input_points(indicies_in_both_lateral_and_station_agreement,2),'k.','MarkerSize',10)
    end


    agreement_indicies(ith_vector,indicies_in_both_lateral_and_station_agreement) = 1;

    agreements(ith_vector) = sum(indicies_in_lateral_agreement,2);

end
end % Ends fcn_INTERNAL_findAgreementsOfPointsToFits


