function [agreement_indices,orothogonal_dist_btw_test_points_and_cubic_curve] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, fittedParameters, tolerance, varargin)
%% fcn_geometry_findAgreementsOfPointsToCubicPoly 
%
% Given a set of XY points, source points, and fitted parameters, finds the
% indicies of the points that are within a transverse_tolerance distance
% away from the cubic polynomial curve, while keeping within
% station_tolerance (if its specified) distance from each point within the
% cluster centered at the base_point_index.
% 
% Note: if station tolerance is specified, the first point is used as the
% "anchor" point wherein the station tolerance is tested from that 1st
% point to the 2nd, then the 2nd to the 3rd, etc.
%
% FORMAT:
%
% agreement_indices = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, fittedParameters, tolerance, (fig_num))
%
% INPUTS:
%
%      points: a Nx2 vector where N is the number of points, but at least 2
%      rows.
%
%      fittedParameters: These are fitted parameters of the cubic
%      polynomial curve. The format is defined in
%      fcn_geometry_fillEmptyDomainStructure 
%
%         Namely:
% 
%             [ 
%              base_point_x, 
%              base_point_y, 
%              flag_is_forward_x,
%              x_Length,
%              coeff_x^0, % cubic term coeffiecient
%              coeff_x^1, % square term coeffiecient
%              coeff_x^2, % linear term coeffiecient
%              coeff_x^3, % constant term coeffiecient
%             ]
%
%      tolerance: the offset, in meters, between the cubic fit and the test
%      points such that points within this tolerance are in agreement. If
%      the offset to a point is larger than this, then the point is not
%      within agreement. If tolerance is entered as a 2x1 or 1x2, then this
%      specifies the tolerance in St coordinates, namely first in the
%      station direction, and then in the transverse direction. For
%      example, an entry of [3 0.02] would have 3 meters tolerance in the
%      station direction, but 0.02 meters tolerance in the transverse
%      direction.
%
%      (OPTIONAL INPUTS)
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
% 2024_06_11 - Aneesh Batchu
% Removed inpolygon method to calculate the inliers. Instead, projection
% distance is used to fidn the points in transverse agreement. This
% improves the speed of the code. 
% 2024_06_12 - S. Brennan
% -- removed weird if statements from
% fcn_INTERNAL_findUnitOrthogonalVectors that were calculating differently
% based on slope thresholds
% -- changed plotting style to make agreement points more clear
% -- changed parameter style to be consistent with line, segment, etc.
% definitions (NOTE: requires changes to other functions!)
% -- changed tolerance style to be consistent with St coordinate definition
% used in other functions
% -- fixed bug in station calculations so that this works

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{end},-1))
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

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 34838; 
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
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(3,4);

        % Check the points input to be length greater than or equal to 2
        fcn_DebugTools_checkInputsToFunctions(...
            points, '2column_of_numbers',[2 3]);

        % % Check the source points input to be length greater than or equal to 2
        % fcn_DebugTools_checkInputsToFunctions(...
        %     source_points, '2column_of_numbers',[2 3]);

        % Check the transverse_tolerance input is a positive single number
        fcn_DebugTools_checkInputsToFunctions(tolerance, 'positive_1or2column_of_numbers',1);

        % % Check the transverse_tolerance input is a positive single number
        % fcn_DebugTools_checkInputsToFunctions(station_tolerance, 'positive_1column_of_numbers',1);
    end
end

% % Does user want to specify base_point_index?
% current_combo = [];
% if 4<= nargin
%     temp = varargin{1};
%     if ~isempty(temp)
%         current_combo = temp;
%     end
% end


% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (4<= nargin)
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

% Extract the cubic fit parameters
cubic_fit = fliplr(fittedParameters(1,5:8));

% Extract the tolerances
if length(tolerance)>1
    station_tolerance = tolerance(1);
    transverse_tolerance = tolerance(2);
else
    station_tolerance = [];
    transverse_tolerance = tolerance;
end

% % Remove the source points from the points
% N_points = length(points(:,1)); 

% Find the points excluding source points
% points_excluding_source_points = setdiff(points, points(current_combo,:), 'rows');

% Find the points in x - domain
% points = find(points >= source_points(1,1) & points <= source_points(end,1));

% Substitute the points to find y values using fitted parameters. The y
% values of the points are different from y coordinates of original points
[y_values_of_points_calculated, slopes_at_each_point] = fcn_INTERNAL_findSlopesAtEachPoint(points(:,1), cubic_fit);

% Calculated points (y values) using fitted parameters. These are not
% original points, but the points that are ON the curve in the form of:
%  [X_original, Y_calculated]
points_calculated = [points(:,1), y_values_of_points_calculated]; 

% Find the unit orthogonal vectors at each calculated point based on the
% slope of the cubic polynomial
unit_orthogonal_vectors = fcn_INTERNAL_findUnitOrthogonalVectors(slopes_at_each_point);

% Check results?
if 1==flag_do_debug
    figure(debug_fig_num)
    clf;
    hold on;
    grid on;

    % plot the input points
    plot(points(:,1), points(:,2), 'k.', 'MarkerSize',30)


    % plot the y_values_of_points_calculated
    plot(points(:,1), y_values_of_points_calculated, 'r.', 'MarkerSize',15)

    % Plot the unit vectors
    quiver(points(:,1),points(:,2),unit_orthogonal_vectors(:,1),unit_orthogonal_vectors(:,2),0,'r','Linewidth',2);

end

% Find vectors from curve to the test points
vectors_from_curve_to_test_points = points - points_calculated; 

% Projection of vectors_calc_original on unit_orthogonal_vectors. Simply,
% approximate perpendicular distance between the test point and cubic
% polynomial.
orothogonal_dist_btw_test_points_and_cubic_curve = sum(vectors_from_curve_to_test_points.*unit_orthogonal_vectors,2); 

% Find the indices in transverse agreement
flags_in_transverse_agreement = abs(orothogonal_dist_btw_test_points_and_cubic_curve) <= transverse_tolerance;
indices_in_transverse_agreement = find(flags_in_transverse_agreement==1); 

% Check results?
if 1==flag_do_debug
    figure(debug_fig_num)

    % plot the agreement points
    plot(points(indices_in_transverse_agreement,1), points(indices_in_transverse_agreement,2), 'go', 'MarkerSize',15, 'LineWidth',3)

end

% If the station distance is given
if ~isempty(station_tolerance)
    station_agreement_indicies = fcn_INTERNAL_findIndicesInStationAgreement(points_calculated(indices_in_transverse_agreement,:), station_tolerance);
    agreement_indices = indices_in_transverse_agreement(station_agreement_indicies);
else
    agreement_indices = indices_in_transverse_agreement;
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
    axis equal
    xlabel('X [m]')
    ylabel('Y [m]')
    
    % plot the test points
    plot(points(:,1), points(:,2), 'k.', 'MarkerSize',30)


    % plot the cubic polynomial
    x_fit = linspace(min(points(:,1)), max(points(:,1)), 100);
    y_fit = polyval(cubic_fit, x_fit);
    plot(x_fit, y_fit, 'k-', 'LineWidth', 2);

    % plot the points in agreement of the cubic polynomial curve
    plot(points(agreement_indices,1), points(agreement_indices,2), 'go', 'MarkerSize',15, 'LineWidth',2)

    quiver(points_calculated(:,1),points_calculated(:,2),vectors_from_curve_to_test_points(:,1),vectors_from_curve_to_test_points(:,2),0,'r','Linewidth',2);

    quiver(points_calculated(:,1),points_calculated(:,2),unit_orthogonal_vectors(:,1),unit_orthogonal_vectors(:,2),0,'m','Linewidth',1);

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


%% fcn_INTERNAL_findSlopesAtEachPoint
function [y_calculated_values, slopes_at_each_calculated_point] = fcn_INTERNAL_findSlopesAtEachPoint(x_coordinates, fittedParameters)
% Find the y coordinates of interpolated source points by substituting x
% coordinates of interpolated source points in cubic polynomial using
% "polyval"
y_calculated_values = polyval(fittedParameters, x_coordinates);

% Calculate squares and cubes of x_interpolated_source_points for speed
% squares_x_interpolated_source_points = [x_interpolated_source_points, x_interpolated_source_points.^2];

% % Find the slopes (first derivative) of cubic polynomial at each test source point
slopes_at_each_calculated_point = 3*fittedParameters(1,1)*x_coordinates.^2 + 2*fittedParameters(1,2)*x_coordinates + fittedParameters(1,3);
% Find the slopes (first derivative) of cubic polynomial at each test source point
% slopes_at_each_interpolated_source_point = 3*fittedParameters(1,1)*squares_x_interpolated_source_points(:,2) + 2*fittedParameters(1,2)*squares_x_interpolated_source_points(:,1) + fittedParameters(1,3);

% round off the slopes to a 4th decimal
slopes_at_each_calculated_point = round(slopes_at_each_calculated_point,4);

end % Ends fcn_INTERNAL_findSlopesAtEachPoint


%% fcn_INTERNAL_findUnitOrthogonalVectors
function unit_orthogonal_vectors = fcn_INTERNAL_findUnitOrthogonalVectors(slopes_at_each_calculated_point)

% Find angle of inclination to find the unit_tangent_vector
theta = atan(slopes_at_each_calculated_point); 
% theta = atan2(slopes_at_each_test_source_point, interpolated_source_points(:,1)); 

% Calculate unit tangent vectors of all test source points
unit_tangent_vectors = [cos(theta), sin(theta)];

% Find the orthogonal vectors
unit_orthogonal_vectors = unit_tangent_vectors*[0 1; -1 0];

end % ends fcn_INTERNAL_findUnitOrthogonalVectors

%% fcn_INTERNAL_findIndicesInStationAgreement
function agreement_indices = fcn_INTERNAL_findIndicesInStationAgreement(points_in_transverse_agreement, station_tolerance)

% Distance calculation from base point to other points
base_point = points_in_transverse_agreement(1,:);
station_distances_of_points_in_transverse_agreement = sum((base_point - points_in_transverse_agreement).^2,2).^0.5;

% Sort the station distances and find those in agreement with station
% tolerance
agreement_indices = ...
    fcn_geometry_findPointsInSequence(...
    station_distances_of_points_in_transverse_agreement, ...
    1, ...
    station_tolerance, -1);

% % indices in both transverse and station agreement
% indices_in_both_transverse_and_station_agreement = indices_in_transverse_agreement(indices_in_station_agreement);
% 
% agreement_indices = indices_in_both_transverse_and_station_agreement;

end % Ends fcn_INTERNAL_findIndicesInStationAgreement
