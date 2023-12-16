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
 
% Find all possible 2-point permutations
N_points = length(points(:,1));
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

% Find the agreements between points and every single line fit
[agreement_indicies] = fcn_INTERNAL_findAgreementsOfPointsToFits(points,unit_projection_vectors, combos_paired, transverse_tolerance);


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

    % Plot the best fit


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
function [agreement_indicies] = fcn_INTERNAL_findAgreementsOfPointsToFits(input_points,unit_projection_vectors, combos_paired, threshold)

N_combos = length(unit_projection_vectors(:,1));
N_points = length(input_points(:,1));

agreements = zeros(N_combos,1);
agreement_indicies = zeros(N_combos,N_points);

for ith_vector = 1:N_combos
    unit_orthogonal_vector = unit_projection_vectors(ith_vector,:)*[0 1; -1 0];
    base_projection_vectors = input_points - input_points(combos_paired(ith_vector,1),:);
    lateral_distances = sum(unit_orthogonal_vector.*base_projection_vectors,2);

    indicies_in_agreement = (abs(lateral_distances)<threshold)';
    agreement_indicies(ith_vector,:) = indicies_in_agreement;

    agreements(ith_vector) = sum(indicies_in_agreement,2);

end
end % Ends fcn_INTERNAL_findAgreementsOfPointsToFits


