% function [agreement_indices, polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, source_points, fittedParameters, transverse_tolerance, base_point_index, station_tolerance, varargin)
function [agreement_indices, polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly_slow_method(points, source_points, fittedParameters, transverse_tolerance, varargin)
%% fcn_geometry_findAgreementsOfPointsToCubicPoly 
%
% Given a set of XY points, source points, and fitted parameters, finds the
% indicies of the points that are within a transverse_tolerance distance
% away from the cubic polynomial curve, while keeping within
% station_tolerance (if its specified) distance from each point within the
% cluster centered at the base_point_index.
%
% FORMAT:
%
% agreement_indices = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, source_points, fittedParameters, transverse_tolerance, (station_tolerance), (total_points_including_source_points), (fig_num))
%
% INPUTS:
%
%      points: a Nx2 vector where N is the number of points, but at least 2
%      rows.
%
%      source_points: These are the points that are used to fit a cubic
%      polynomial curve using "polyfit"
%
%      fittedParameters: These are fitted parameters of the cubic
%      polynomial curve
%
%      transverse_tolerance: the orthogonal distance between the points and
%      the linear vector fit that indicate whether a point "belongs" to the
%      fit. A point belongs to the fit if the transverse distance is less
%      than or equal to the transverse_tolerance
%
%      (OPTIONAL INPUTS)
%
%      station_tolerance: the projection distance between the points in a
%      curve fit, along the direction of the line, that indicate whether a
%      point "belongs" to the circle fit (if distance is less than or equal
%      to the tolerance), or is "outside" the fit (if distance is greater
%      than the tolerance).
%
%      total_points_including_source_points: To get a better domain, the
%      source points are interpolated. This is the input for total points
%      you want including the source points for getting a better transverse
%      tolerance.
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      agreement_indicies: the indicies of the points that are within
%      agreement of the best-fit parameters, given the transverse and
%      station tolerance settings. 
%
%      polygon_vertices: vertices of the transverse domain.     
% 
% DEPENDENCIES:
%  
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_findAgreementsOfPointsToCubicPoly
% for a full test suite.

% This function was written on 2024_05_10 by Aneesh Batchu
% Questions or comments? abb6486@psu.edu or sbrennan@psu.edu

% Revision history:
% 2024_05_10 - Aneesh Batchu
% -- wrote the code
% 2024_05_13 - Aneesh Batchu
% -- "polygon_vertices" are included as the output to plot the polygon
% domain
% 2024_05_17 - Aneesh Batchu
% -- Added station_tolerance to the code. This code finds the agreement
% indices that are not just in transverse agreement but also in station
% agreement
% 2024_05_29 - Aneesh Batchu
% -- Added an internal function to interpolate the source points,
% generating more domain points and eventually obtaining a better domain
% box around the area of interest.
% -- Added polyfit to the upper and lower boundary points to get a better
% domain
% 2024_05_30 - Aneesh Batchu 
% -- If the slope at any point (source/interpolated) point is zero, the
% corresponding unit tangent vectors is multiplied with transverse
% tolerance (instead of orthogonal vectors to find the domain boundaries)
% -- Added a conditional statement 
% (&& length(indices_in_transverse_agreement)>=2) for finding the points
% that are within the station_tolerance limit.
% -- Removed base_point_index from the inputs
% 2024_06_03 - Aneesh Batchu
% -- Functionalized unit orthogonal vectors and domain
% vertices(fcn_INTERNAL_findUnitOrthogonalVectors and
% fcn_INTERNAL_findDomainVertices) 
% -- Added "total_points_including_source_points" as one of the inputs
% 2024_06_05 - Aneesh Batchu
% -- polygon_vertices are changed when station tolerance is given as the
% input
% -- Functionalized the slope finding method
% "fcn_INTERNAL_findSlopesAtEachPoint"

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==7 && isequal(varargin{end},-1))
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
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(4,7);

        % Check the points input to be length greater than or equal to 2
        fcn_DebugTools_checkInputsToFunctions(...
            points, '2column_of_numbers',[2 3]);

        % Check the source points input to be length greater than or equal to 2
        fcn_DebugTools_checkInputsToFunctions(...
            source_points, '2column_of_numbers',[2 3]);

        % Check the transverse_tolerance input is a positive single number
        fcn_DebugTools_checkInputsToFunctions(transverse_tolerance, 'positive_1column_of_numbers',1);

        % % Check the transverse_tolerance input is a positive single number
        % fcn_DebugTools_checkInputsToFunctions(station_tolerance, 'positive_1column_of_numbers',1);
    end
end

% % Does user want to specify base_point_index?
% base_point_index = [];
% if 5<= nargin
%     temp = varargin{1};
%     if ~isempty(temp)
%         base_point_index = temp;
%     end
% end
% 

% Does user want to specify station_tolerance?
station_tolerance = [];
if 5<= nargin
    temp = varargin{1};
    if ~isempty(temp)
        station_tolerance = temp;
    end
end

% Does user specify total_points_including_source_points?
total_points_including_source_points = 20;
if (6<= nargin)
    temp = varargin{2};
    if ~isempty(temp)
        total_points_including_source_points = temp;
        if total_points_including_source_points<4
            error('The input points_required_for_agreement must be greater than or equal to 4.')
        end
    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (7<= nargin)
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

% Generate N points by interpolating the x-coordinates of source points
x_interpolated_source_points = fcn_INTERNAL_interpolateSourcePoints(source_points, total_points_including_source_points);

% Find the y coordinates of interpolated source points by substituting x
% coordinates of interpolated source points in cubic polynomial using
% "polyval"
[y_interpolated_source_points, slopes_at_each_interpolated_source_point] = fcn_INTERNAL_findSlopesAtEachPoint(x_interpolated_source_points, fittedParameters);

% Find interpolated source points matrix by putting x and y coordinates of
% interpolated source points together. 
interpolated_source_points = [x_interpolated_source_points, y_interpolated_source_points];

% Find the unit orthogonal vectors at each interpolated source points based
% on the slope of the cubic polynomial
unit_orthogonal_vectors = fcn_INTERNAL_findUnitOrthogonalVectors(slopes_at_each_interpolated_source_point);

% Find the polygon vertices of the domain based on given transverse
% tolerance
[polygon_vertices_transverse_tolerance,upper_boundary_points,lower_boundary_points, y_values_upperboundary, y_values_lowerboundary] = fcn_INTERNAL_findDomainVertices(interpolated_source_points, unit_orthogonal_vectors, transverse_tolerance); 

% Determine which points lie inside the polygonal region (domain)
inlier_indices = inpolygon(points(:,1), points(:,2), polygon_vertices_transverse_tolerance(:,1), polygon_vertices_transverse_tolerance(:,2));
% outliers_indices = ~inliers_indices;

% Find the transverse agreement indices based on inlier_indices
indices_in_transverse_agreement = find(inlier_indices==1);
% indices_in_transverse_agreement
% If the station distance is given
if ~isempty(station_tolerance) && length(indices_in_transverse_agreement)>=2

    agreement_indices = fcn_INTERNAL_findIndicesInStationAgreement(points, indices_in_transverse_agreement, station_tolerance);

    if 1 < length(agreement_indices)
        if total_points_including_source_points < length(agreement_indices)
            total_points_including_agreement_points = length(agreement_indices);
        else
            total_points_including_agreement_points = total_points_including_source_points;
        end

        % Interpolate x-coordinates of points within the station tolerance limit
        x_interpolated_agreement_points_station_tolerance = fcn_INTERNAL_interpolateSourcePoints(points(agreement_indices,:), total_points_including_agreement_points);

        % Find the y coordinates and slopes at each interpolated agreement point
        [y_interpolated_agreement_points_station_tolerance, slopes_at_each_interpolated_agreement_point] = fcn_INTERNAL_findSlopesAtEachPoint(x_interpolated_agreement_points_station_tolerance, fittedParameters);

        % Find the unit orthogonal vectors at each interpolated agreement
        % points based on the slope of the cubic polynomial
        unit_orthogonal_vectors_interpolated_agreement_points = fcn_INTERNAL_findUnitOrthogonalVectors(slopes_at_each_interpolated_agreement_point);

        % Find interpolated agreement points matrix by putting x and y coordinates of
        % interpolated agreement points together.
        interpolated_agreement_points = [x_interpolated_agreement_points_station_tolerance, y_interpolated_agreement_points_station_tolerance];

        % Find the polygon vertices of the domain based on given transverse
        % tolerance and station tolerance
        [polygon_vertices_station_tolerance,upper_boundary_points,lower_boundary_points, y_values_upperboundary, y_values_lowerboundary] = fcn_INTERNAL_findDomainVertices(interpolated_agreement_points, unit_orthogonal_vectors_interpolated_agreement_points, transverse_tolerance);

        polygon_vertices = polygon_vertices_station_tolerance;

    else
        
        agreement_indices = [];
        polygon_vertices = [nan, nan]; 
  
    end
    
else
    agreement_indices = indices_in_transverse_agreement;
    polygon_vertices = polygon_vertices_transverse_tolerance; 
end

% sort the agreement indices
agreement_indices = sort(agreement_indices); 

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
    
    figure(fig_num)
    hold on;
    grid on;
    % axis equal
    xlabel('X [m]')
    ylabel('Y [m]')
    
    % plot the test points
    plot(points(:,1), points(:,2), 'c.', 'MarkerSize',30)

    % plot the source points of the cubic polynomial curve
    plot(source_points(:,1), source_points(:,2), 'b.', 'MarkerSize',30)
    
    % Plot the fitted polynomial
    x_fit = linspace(min(source_points(:,1)), max(source_points(:,1)), 100);
    y_fit = polyval(fittedParameters, x_fit);
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

    % % % plot the unit tangent vectors
    % % end_points_tangent = interpolated_source_points + unit_tangent_vectors;
    % % plot(unit_tangent_vectors(:,1), unit_tangent_vectors(:,2), 'r.', 'MarkerSize',30)
    % % end_points_quiver_tangent = end_points_tangent - interpolated_source_points;
    % % quiver(0*interpolated_source_points(:,1),0*interpolated_source_points(:,2),unit_tangent_vectors(:,1),unit_tangent_vectors(:,2),0,'b','Linewidth',2);
    % if ~isempty(station_tolerance) && length(indices_in_transverse_agreement)>=2
    if ~isempty(station_tolerance) && length(indices_in_transverse_agreement)>=2 && 1 < length(agreement_indices)

        % plot the upper boundary points
        plot(upper_boundary_points(:,1), upper_boundary_points(:,2), 'g.', 'MarkerSize',30)

        % Show the unit vectors projecting from curve to upper boundary points
        % in green
        end_points_quiver_upper = upper_boundary_points - interpolated_agreement_points;
        quiver(interpolated_agreement_points(:,1),interpolated_agreement_points(:,2),end_points_quiver_upper(:,1),end_points_quiver_upper(:,2),0,'g','Linewidth',2);

        % plot the fitted upper boundary curve
        plot(interpolated_agreement_points(:,1), y_values_upperboundary, 'g-', 'LineWidth', 2);

        % plot the lower boundary points
        plot(lower_boundary_points(:,1), lower_boundary_points(:,2), 'm.', 'MarkerSize',30)

        % Show the unit vectors projecting from curve to lower boundary points
        % in green
        end_points_quiver_lower = lower_boundary_points - interpolated_agreement_points;
        quiver(interpolated_agreement_points(:,1),interpolated_agreement_points(:,2),end_points_quiver_lower(:,1),end_points_quiver_lower(:,2),0,'m','Linewidth',2);

        % plot the fitted lower boundary curve
        plot(interpolated_agreement_points(:,1), y_values_lowerboundary, 'm-', 'LineWidth', 2);

    else

        % plot the interpolated_source_points
        % plot(interpolated_source_points(:,1), y_interpolated_source_points, 'b.', 'MarkerSize',30)
        plot(interpolated_source_points(:,1), interpolated_source_points(:,2), 'go', 'MarkerSize',30)

        % plot the upper boundary points
        plot(upper_boundary_points(:,1), upper_boundary_points(:,2), 'g.', 'MarkerSize',30)

        % Show the unit vectors projecting from curve to upper boundary points
        % in green
        end_points_quiver_upper = upper_boundary_points - interpolated_source_points;
        quiver(interpolated_source_points(:,1),interpolated_source_points(:,2),end_points_quiver_upper(:,1),end_points_quiver_upper(:,2),0,'g','Linewidth',2);

        % plot the fitted upper boundary curve
        plot(interpolated_source_points(:,1), y_values_upperboundary, 'g-', 'LineWidth', 2);

        % plot the lower boundary points
        plot(lower_boundary_points(:,1), lower_boundary_points(:,2), 'm.', 'MarkerSize',30)

        % Show the unit vectors projecting from curve to lower boundary points
        % in green
        end_points_quiver_lower = lower_boundary_points - interpolated_source_points;
        quiver(interpolated_source_points(:,1),interpolated_source_points(:,2),end_points_quiver_lower(:,1),end_points_quiver_lower(:,2),0,'m','Linewidth',2);

        % plot the fitted lower boundary curve
        plot(interpolated_source_points(:,1), y_values_lowerboundary, 'm-', 'LineWidth', 2);

    end
    % drawpolygon function is used to set the domain
    drawpolygon(gca,"Position",polygon_vertices);

    % plot the points in agreement of the cubic polynomial curve
    plot(points(agreement_indices,1), points(agreement_indices,2), 'r.', 'MarkerSize',40)

    % % plot the points in agreement of the cubic polynomial curve
    % plot(points(agreement_indices,1), points(agreement_indices,2), 'ro', 'MarkerSize',30)

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end

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

function interpolated_source_points = fcn_INTERNAL_interpolateSourcePoints(source_points, n_points)

% Sort the x_coordinates of source points
x_coordinates_source_points = sort(source_points(:,1),1);

% Generate the positions for the 4 x_coordinates of source points within
% the (n_points: 2 to (n_points -1)) 
indices = round(linspace(1, n_points, length(x_coordinates_source_points)));

% Initialize the result vector with NaNs (or any placeholder)
interpolated_source_points = NaN(1, n_points);

% Place the given numbers at the calculated indices
interpolated_source_points(indices) = x_coordinates_source_points;

% Find the indices of the NaNs (to interpolate)
nan_indices = find(isnan(interpolated_source_points));

% Interpolate the NaNs
interpolated_source_points(nan_indices) = interp1(indices, x_coordinates_source_points, nan_indices, 'makima'); 

% The interpolated points are transposed
interpolated_source_points = interpolated_source_points';

end

function [y_interpolated_source_points, slopes_at_each_interpolated_source_point] = fcn_INTERNAL_findSlopesAtEachPoint(x_interpolated_source_points, fittedParameters)
% Find the y coordinates of interpolated source points by substituting x
% coordinates of interpolated source points in cubic polynomial using
% "polyval"
y_interpolated_source_points = polyval(fittedParameters, x_interpolated_source_points);

% Calculate squares and cubes of x_interpolated_source_points for speed
% squares_x_interpolated_source_points = [x_interpolated_source_points, x_interpolated_source_points.^2];

% % Find the slopes (first derivative) of cubic polynomial at each test source point
slopes_at_each_interpolated_source_point = 3*fittedParameters(1,1)*x_interpolated_source_points.^2 + 2*fittedParameters(1,2)*x_interpolated_source_points + fittedParameters(1,3);
% Find the slopes (first derivative) of cubic polynomial at each test source point
% slopes_at_each_interpolated_source_point = 3*fittedParameters(1,1)*squares_x_interpolated_source_points(:,2) + 2*fittedParameters(1,2)*squares_x_interpolated_source_points(:,1) + fittedParameters(1,3);


% round off the slopes to a 4th decimal
slopes_at_each_interpolated_source_point = round(slopes_at_each_interpolated_source_point,4);

end

function unit_orthogonal_vectors = fcn_INTERNAL_findUnitOrthogonalVectors(slopes_at_each_test_source_point)

% Find angle of inclination to find the unit_tangent_vector
theta = atan(slopes_at_each_test_source_point); 
% theta = atan2(slopes_at_each_test_source_point, interpolated_source_points(:,1)); 

% Unit tangent vectors of all test source points
unit_tangent_vectors = [sin(theta), cos(theta)]; 

% Pre-allocate unit_orthogonal_vectors with zeros
unit_orthogonal_vectors = zeros(size(unit_tangent_vectors));

% If the slopes are negative, rotate the unit tangent vector in clockwise
% direction.
indices_slope_is_negative = slopes_at_each_test_source_point < 0;
unit_orthogonal_vectors(indices_slope_is_negative,:) = unit_tangent_vectors(indices_slope_is_negative,:)*[0 -1; 1 0];

% If the slopes are positive, rotate the unit tangent vector in anti
% clockwise direction
indices_slope_is_positive = slopes_at_each_test_source_point > 0;
unit_orthogonal_vectors(indices_slope_is_positive,:) = unit_tangent_vectors(indices_slope_is_positive,:)*[0 1; -1 0];

% If the slopes are zero, use tangent vectors to give transverse tolerance
% to generate the domain box
indices_slope_is_zero = slopes_at_each_test_source_point == 0;
unit_orthogonal_vectors(indices_slope_is_zero,:) = unit_tangent_vectors(indices_slope_is_zero,:);


end

function [polygon_vertices,upper_boundary_points,lower_boundary_points, y_values_upperboundary, y_values_lowerboundary]  = fcn_INTERNAL_findDomainVertices(interpolated_source_points, unit_orthogonal_vectors, transverse_tolerance)

% Find the points on upper and lower boundaries
upper_boundary_points = interpolated_source_points + unit_orthogonal_vectors*transverse_tolerance;
lower_boundary_points = interpolated_source_points - unit_orthogonal_vectors*transverse_tolerance;

% Fit the upper and lower boundary points using polyfit
fit_upper_boundary_parameters = polyfit(upper_boundary_points(:,1), upper_boundary_points(:,2), 3);
fit_lower_boundary_parameters = polyfit(lower_boundary_points(:,1), lower_boundary_points(:,2), 3);

% Plot the upper and lower fitted polynomial to fix the vertices of domain
y_values_upperboundary = polyval(fit_upper_boundary_parameters, interpolated_source_points(:,1));
y_values_lowerboundary = polyval(fit_lower_boundary_parameters, interpolated_source_points(:,1));

% % Define the domain. These are the boundary vertices
% polygon_vertices = [upper_boundary_points; lower_boundary_points(end:-1:1,:)];

% Define the domain. These are the boundary vertices
upper_polygon_vertices = [interpolated_source_points(:,1), y_values_upperboundary];
lower_polygon_vertices = [interpolated_source_points(:,1), y_values_lowerboundary];
polygon_vertices = [upper_polygon_vertices; lower_polygon_vertices(end:-1:1,:)];

end

% function agreement_indices = fcn_INTERNAL_findIndicesInStationAgreement(points, indices_in_transverse_agreement, base_point_index, station_tolerance)
function agreement_indices = fcn_INTERNAL_findIndicesInStationAgreement(points, indices_in_transverse_agreement, station_tolerance)
base_point_index = indices_in_transverse_agreement(1,1);

% Grab only the points in transverse agreement
points_in_transverse_agreement = points(indices_in_transverse_agreement,:);

% Find index of the source point in the rearranged list
index_source_point_in_transverse_agreement = find(indices_in_transverse_agreement == base_point_index,1);

% Find the length of the points in transverse agreement
% N = length(points_in_transverse_agreement(:,1));

% Find the point pairs of the points in transverse agreement to compute
% the station distances between them.
% point_pairs = [1:N-1; 2:N]';

% Find the differences (vectors) of point_pairs(:,2) and
% point_pairs(:,1)
% diff_between_pts_in_points_pair = points_in_transverse_agreement(point_pairs(:,2),:) - points_in_transverse_agreement(point_pairs(:,1),:);
diff_between_pts_in_points_pair = diff(points_in_transverse_agreement); 

% Station distance calculation
station_distances_of_points_in_transverse_agreement = sum(diff_between_pts_in_points_pair.^2,2).^0.5;
% index_source_point_in_transverse_agreement
% length_of_input_distances = length(station_distances_of_points_in_transverse_agreement)

% Sort the station distances and find those in agreement with station
% tolerance
indices_in_station_agreement = ...
    fcn_geometry_findPointsInSequence(...
    station_distances_of_points_in_transverse_agreement, ...
    index_source_point_in_transverse_agreement, ...
    station_tolerance, -1);

indices_in_both_transverse_and_station_agreement = indices_in_transverse_agreement(indices_in_station_agreement);

agreement_indices = indices_in_both_transverse_and_station_agreement;

end
