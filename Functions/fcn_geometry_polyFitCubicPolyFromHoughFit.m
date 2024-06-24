function regression_domain = fcn_geometry_polyFitCubicPolyFromHoughFit(Hough_domain, varargin)
%% fcn_geometry_polyFitCubicPolyFromHoughFit
%
% Given a domain containing a set of points that are matched to an arc via
% a Hough vote, finds the polyfit and domain box.
%
% FORMAT:
%
% regression_domain = fcn_geometry_polyFitCubicPolyFromHoughFit(Hough_domain, (best_fit_domain_box_projection_distance), (fig_num))
%
% INPUTS:
%
%      Hough_domain: a structure th
% at records details of the domain of
%      fitting. See fcn_geometry_fillEmptyDomainStructure for details.
%
%      (OPTIONAL INPUTS)
%
%      best_fit_domain_box_projection_distance: entered as an 1x1, this
%      is the distance from the curve fit, in the transverse direction, to
%      project in both the positive and negative directions to produce the
%      best_fit_domain_box.
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      regression_domain: a structure that records details of the domain of
%      fitting. See fcn_geometry_fillEmptyDomainStructure for details.
%
% DEPENDENCIES:
%
%      fcn_geometry_fillEmptyDomainStructure
%      fcn_geometry_plotFitDomains
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_polyFitCubicPolyFromHoughFit
% for a full test suite.
%
% This function was written on 2024_06_08 by Aneesh Batchu
% Questions or comments? abb6486@psu.edu or sbrennan@psu.edu 

% Revision history:
% 2024_06_08 - Aneesh Batchu
% -- wrote the code originally

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
        narginchk(1,3);
    end
end

% Does user want to specify best_fit_domain_box_projection_distance?
best_fit_domain_box_projection_distance = [];
if (2<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        best_fit_domain_box_projection_distance = temp;
    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if  (0==flag_max_speed) && (3<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

%% Solve for the cubic polynomial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialze the domain structure for output
regression_domain = fcn_geometry_fillEmptyDomainStructure;

% Pull out key variables
best_fit_type            = Hough_domain.best_fit_type;
points_in_domain         = Hough_domain.points_in_domain;
best_fit_parameters      = Hough_domain.best_fit_parameters;
% best_fit_domain_box      = Hough_domain.best_fit_domain_box;

if isequal(best_fit_type, 'Hough cubic polynomial')
    regression_domain.best_fit_type = 'Polyfit cubic polynomial';
end

% Set the regression domain's points_in_domain
regression_domain.points_in_domain = points_in_domain;

% Find fitted curve - use the function "polyfit"
cubicPoly_fittedParameters = polyfit(points_in_domain(:,1), points_in_domain(:,2), 3);

% Regression domain fitted parameters
fittedParameters = [best_fit_parameters(1,1:4) fliplr(cubicPoly_fittedParameters)];  

regression_domain.best_fit_parameters = fittedParameters;

if 20 > length(points_in_domain(:,1))
    n_points = 20;
else
    n_points = length(points_in_domain(:,1));
end

% Generate N points by interpolating the x-coordinates of points in domain
xcoor_interpolated_points_in_domain = fcn_INTERNAL_interpolatePointsInDomain(points_in_domain, n_points);

% Find the y coordinates of interpolated points in domain by substituting x
% coordinates of interpolated points in domain in cubic polynomial using
% "polyval"
[ycoor_interpolated_points_in_domain, slopes_at_each_interpolated_points_in_domain] = fcn_INTERNAL_findSlopesAtEachPoint(xcoor_interpolated_points_in_domain, cubicPoly_fittedParameters);

% Find the unit orthogonal vectors at each interpolated points in domain
% based on the slope of the cubic polynomial
unit_orthogonal_vectors = fcn_INTERNAL_findUnitOrthogonalVectors(slopes_at_each_interpolated_points_in_domain);

% Find interpolated points in domain matrix by putting x and y coordinates
% of interpolated points in domain together.
interpolated_points_in_domain = [xcoor_interpolated_points_in_domain, ycoor_interpolated_points_in_domain];

% Find the polygon vertices of the domain based on given transverse
% tolerance
[best_fit_polygon_vertices,~,~, ~, ~]  = fcn_INTERNAL_findDomainVertices(interpolated_points_in_domain, unit_orthogonal_vectors, best_fit_domain_box_projection_distance);

% Create domain box based on the best fit polygon vertices
domainShape = fcn_geometry_domainBoxByType('cubic polynomial', best_fit_polygon_vertices,-1);

% Best fit domain box
regression_domain.best_fit_domain_box = domainShape;

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
    axis equal;

    % Plot the fits    
    % current_color = fcn_geometry_fillColorFromNumberOrName(1,'points',-1);

    % Plot the associated_points_in_domain
    % plot(associated_points_in_domain(:,1),associated_points_in_domain(:,2),'.','MarkerSize',5,'Color',current_color);

    % Plot the domains
    fcn_geometry_plotFitDomains(regression_domain, fig_num);

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


function xcoor_interpolated_points_in_domain = fcn_INTERNAL_interpolatePointsInDomain(points_in_domain, n_points)

% Sort the x_coordinates of source points
x_coordinates_points_in_domain = sort(points_in_domain(:,1),1);

% Generate the positions for the 4 x_coordinates of source points within
% the (n_points: 2 to (n_points -1)) 
indices = round(linspace(1, n_points, length(x_coordinates_points_in_domain)));

% Initialize the result vector with NaNs (or any placeholder)
xcoor_interpolated_points_in_domain = NaN(1, n_points);

% Place the given numbers at the calculated indices
xcoor_interpolated_points_in_domain(indices) = x_coordinates_points_in_domain;

% Find the indices of the NaNs (to interpolate)
nan_indices = find(isnan(xcoor_interpolated_points_in_domain));

% Interpolate the NaNs
xcoor_interpolated_points_in_domain(nan_indices) = interp1(indices, x_coordinates_points_in_domain, nan_indices, 'makima'); 

% The interpolated points are transposed
xcoor_interpolated_points_in_domain = xcoor_interpolated_points_in_domain';

end

% function xcoor_interpolated_points_in_domain = fcn_INTERNAL_interpolateSourcePoints(points_in_domain, n_points)
% 
% % Sort the x_coordinates of points in domain
% x_coordinates_points_in_domain = sort(points_in_domain(:,1),1);
% 
% % Generate the positions for the x_coordinates of points in domain 
% indices = round(linspace(1, n_points, length(x_coordinates_points_in_domain)));
% 
% % Initialize the result vector with NaNs (or any placeholder)
% xcoor_interpolated_points_in_domain = NaN(1, n_points);
% 
% % Place the given numbers at the calculated indices
% xcoor_interpolated_points_in_domain(indices) = x_coordinates_points_in_domain;
% 
% % Find the indices of the NaNs (to interpolate)
% nan_indices = find(isnan(xcoor_interpolated_points_in_domain));
% 
% % Interpolate the NaNs
% xcoor_interpolated_points_in_domain(nan_indices) = interp1(indices, x_coordinates_points_in_domain, nan_indices, 'makima'); 
% 
% % The interpolated points are transposed
% xcoor_interpolated_points_in_domain = xcoor_interpolated_points_in_domain';
% 
% end

function [ycoor_interpolated_points_in_domain, slopes_at_each_interpolated_points_in_domain] = fcn_INTERNAL_findSlopesAtEachPoint(x_interpolated_points_in_domain, fitted_parameters)
% Find the y coordinates of interpolated source points by substituting x
% coordinates of interpolated source points in cubic polynomial using
% "polyval"
ycoor_interpolated_points_in_domain = polyval(fitted_parameters, x_interpolated_points_in_domain);

% Calculate squares and cubes of x_interpolated_source_points for speed
squares_x_interpolated_points_in_domain = [x_interpolated_points_in_domain, x_interpolated_points_in_domain.^2];

% % Find the slopes (first derivative) of cubic polynomial at each test source point
% slopes_at_each_test_source_point = 3*fittedParameters(1,1)*interpolated_source_points(:,1).^2 + 2*fittedParameters(1,2)*interpolated_source_points(:,1) + fittedParameters(1,3);
% Find the slopes (first derivative) of cubic polynomial at each test source point
slopes_at_each_interpolated_points_in_domain = 3*fitted_parameters(1,1)*squares_x_interpolated_points_in_domain(:,2) + 2*fitted_parameters(1,2)*squares_x_interpolated_points_in_domain(:,1) + fitted_parameters(1,3);

% round off the slopes to a 4th decimal
slopes_at_each_interpolated_points_in_domain = round(slopes_at_each_interpolated_points_in_domain,4);

end

function unit_orthogonal_vectors = fcn_INTERNAL_findUnitOrthogonalVectors(slopes_at_each_interpolated_points_in_domain)

% Find angle of inclination to find the unit_tangent_vector
theta = atan(slopes_at_each_interpolated_points_in_domain); 
% theta = atan2(slopes_at_each_test_source_point, interpolated_source_points(:,1)); 

% Unit tangent vectors of all test source points
unit_tangent_vectors = [cos(theta), sin(theta)]; 

% Find the orthogonal vectors
unit_orthogonal_vectors = unit_tangent_vectors*[0 1; -1 0];

end

function [polygon_vertices,upper_boundary_points,lower_boundary_points, y_values_upperboundary, y_values_lowerboundary]  = fcn_INTERNAL_findDomainVertices(interpolated_points_in_domain, unit_orthogonal_vectors, transverse_tolerance)

% Find the points on upper and lower boundaries
upper_boundary_points = interpolated_points_in_domain + unit_orthogonal_vectors*transverse_tolerance;
lower_boundary_points = interpolated_points_in_domain - unit_orthogonal_vectors*transverse_tolerance;

% Fit the upper and lower boundary points using polyfit
fit_upper_boundary_parameters = polyfit(upper_boundary_points(:,1), upper_boundary_points(:,2), 3);
fit_lower_boundary_parameters = polyfit(lower_boundary_points(:,1), lower_boundary_points(:,2), 3);

% Plot the upper and lower fitted polynomial to fix the vertices of domain
y_values_upperboundary = polyval(fit_upper_boundary_parameters, interpolated_points_in_domain(:,1));
y_values_lowerboundary = polyval(fit_lower_boundary_parameters, interpolated_points_in_domain(:,1));

% % Define the domain. These are the boundary vertices
% polygon_vertices = [upper_boundary_points; lower_boundary_points(end:-1:1,:)];

% Define the domain. These are the boundary vertices
upper_polygon_vertices = [interpolated_points_in_domain(:,1), y_values_upperboundary];
lower_polygon_vertices = [interpolated_points_in_domain(:,1), y_values_lowerboundary];
polygon_vertices = [upper_polygon_vertices; lower_polygon_vertices(end:-1:1,:)];

end


