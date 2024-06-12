function domains = fcn_geometry_fitHoughCubicPolynomial(points, transverse_tolerance, varargin)
%% fcn_geometry_fitHoughCubicPolynomial
%
% This function takes the input points and tolerance as the input and
% outputs the fitted parameters, agreement indices, and best fit type in a
% cell array.
%
% FORMAT: 
%
% domains = fcn_geometry_fitHoughCubicPolynomial(points, transverse_tolerance, (station_tolerance), (points_required_for_agreement), (flag_find_only_best_agreement), (total_points_including_source_points), (fig_num)) 
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
%      points_required_for_agreement: the number of points required for an
%      agreement to be valid, with minimum value of 3. Default is 10. If
%      left empty, then the best agreement will always be returned. If a
%      value is given, line fitting will continue by clustering data until
%      there are no fits greater than or equal to
%      points_required_for_agreement. The results of the line fit will be
%      saved in a cell array.
%
%      flag_find_only_best_agreement: set to 1 if want to only keep best
%      agreement. Otherwise, searches will be continued until none are left
%      that have more than points_required_for_agreement. Default is 0, to
%      find all agreements
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
%      fcn_geometry_findAgreementsOfPointsToCubicPoly
%      fcn_geometry_fillEmptyDomainStructure
%      fcn_geometry_fillColorFromNumberOrName
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_fitHoughCubicPolynomial for a
% full test suite.
%
% This function was written on 2024_05_11 by Aneesh Batchu
% Questions or comments? abb6486@psu.edu or sbrennan@psu.edu


% Revision history:
% 2024_05_11 - Aneesh Batchu
% -- Started the code
% 2024_05_12 - Aneesh Batchu
% -- Functionalized a code that finds agreement indices
% 2024_05_14 - Aneesh Batchu
% -- Prepared parameters for domains and formed the domains
% 2024_05_15 - Aneesh Batchu
% -- Added a case for cubic polynomial in fcn_geometry_plotFitDomains
% 2024_05_16 - Aneesh Batchu
% -- Functionalized this code. 
% 2024_05_17 - Aneesh Batchu
% -- Added station tolerance as one of the inputs
% 2024_05_20 - Aneesh Batchu
% -- Added wait bar
% 2024_06_05 - Aneesh Batchu
% -- Added "total_points_including_source_points" as one of the inputs
% 2024_06_05 - Aneesh Batchu
% -- Modified the function "fcn_geometry_findAgreementsOfPointsToCubicPoly"
% to find the agreement indices faster. More details of the function are
% mentioned in the instructions of that function. 

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
    debug_fig_num = 234343; %#ok<NASGU>
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
        narginchk(2,7);

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
points_required_for_agreement = 10;
if (4<= nargin)
    temp = varargin{2};
    if ~isempty(temp)
        points_required_for_agreement = temp;
        if points_required_for_agreement<4
            error('The input points_required_for_agreement must be greater than or equal to 4.')
        end
    end
end

% Does user want to specify flag_find_only_best_agreement?
flag_find_only_best_agreement = 0;
if (5<=nargin)
    temp = varargin{3};
    if ~isempty(temp)
        flag_find_only_best_agreement = temp;
    end
end

% Does user specify total_points_including_source_points?
total_points_after_interpolating_points_in_domain = 20;
if (6<= nargin)
    temp = varargin{4};
    if ~isempty(temp)
        total_points_after_interpolating_points_in_domain = temp;
        if total_points_after_interpolating_points_in_domain<4
            error('The input points_required_for_agreement must be greater than or equal to 4.')
        end
    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (7<=nargin)
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

slow_method = 0; 

% Total number of points
N_points = size(points,1);

% Find all possible 4-point combinations
combos_paired = nchoosek(1:N_points,4);

% How many combinations are there?
N_permutations = size(combos_paired,1);

% Pre-allocation of fittedParameters and agreementIndices for saving
% computation time
fitted_parameters = zeros(N_permutations,4);

%% Step 1: find all the agreement counts, save in array "agreements"

agreements = zeros(N_permutations,1);

% Loop through all the combos, recording agreement with each combo
if 0==flag_max_speed
    h_waitbar = waitbar(0,'Calculating cubic polynomial fits...');
end

for ith_combo = 1:N_permutations

    % fprintf(1,'Checking %.0d of %.0d\n',ith_combo,N_permutations);
    if 0==flag_max_speed
        if 0==mod(ith_combo,1000)
            waitbar(ith_combo/N_permutations);
        end
    end

    % Extract the source points
    test_source_points = points(combos_paired(ith_combo,:),:);

    % Find fitted curve - use the function "polyfit"
    fittedParameters = polyfit(test_source_points(:,1), test_source_points(:,2), 3);

    % Store resulting fitted parameters in "fitted_parameters" matrix
    fitted_parameters(ith_combo,:) = fittedParameters;
    
    % Find agreement indices
    [agreement_indices,~] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, fittedParameters, transverse_tolerance, (ith_combo), (station_tolerance),(-1));
    if 1 == slow_method
        [agreement_indices,~] = fcn_geometry_findAgreementsOfPointsToCubicPoly_slow_method(points, test_source_points, fittedParameters, transverse_tolerance, (station_tolerance), (total_points_after_interpolating_points_in_domain), (-1));
    end
    % Save the count of points in agreement
    agreement_count = length(agreement_indices);
    agreements(ith_combo) = agreement_count;

    % Break out if all points agree. No need to do more calculations!
    if agreement_count==N_points
        break;
    end

end

if 0==flag_max_speed
    close(h_waitbar)
end

%% Step 2: take the top agreements, accumulate them in prep for domains

% Initialize outputs as cell arrays
best_fit_cubic_poly_coefficients   = {};
best_fit_agreement_indicies      = {};
best_fit_source_indicies         = {};
if 1 == slow_method
    best_fit_polygon_vertices        = {};  % comment this
end

% Find best agreements
if 1 == flag_find_only_best_agreement
    N_fits = 0;
    remaining_agreements = agreements;

    [~, best_agreement_index] = max(remaining_agreements);
    if 1 == slow_method
        % Best fit source points
        only_best_fit_source_points = points(combos_paired(best_agreement_index,:),:);
    end
    % Find the best fit parameters
    only_best_fit_fittedCoefficients  = fitted_parameters(best_agreement_index,:);
    
    % Find the indicies in transverse agreement with this best fit
    if 1 == slow_method
        [only_best_agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly_slow_method(points, only_best_fit_source_points, only_best_fit_fittedCoefficients, transverse_tolerance, (station_tolerance), (total_points_after_interpolating_points_in_domain), (-1));
    end
    [only_best_agreement_indices,~] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, only_best_fit_fittedCoefficients, transverse_tolerance, (combos_paired(best_agreement_index,:)), (station_tolerance), (-1));

    % Accumulate results from best fits to prep for domain formation
    N_fits = N_fits+1;
    best_fit_cubic_poly_coefficients{N_fits}  = only_best_fit_fittedCoefficients;
    best_fit_agreement_indicies{N_fits}     = only_best_agreement_indices;
    best_fit_source_indicies{N_fits}        = combos_paired(best_agreement_index,:);
    if 1 == slow_method
        best_fit_polygon_vertices{N_fits}        = polygon_vertices;  % comment this
    end
else
    N_fits = 0;
    remaining_agreements = agreements;

    [best_agreement_count, best_agreement_index] = max(remaining_agreements);

    while best_agreement_count >= points_required_for_agreement
        if 1 == slow_method
            % test source points
            fit_source_points = points(combos_paired(best_agreement_index,:),:);
        end
        % Find the best fit parameters
        fittedCoefficients  = fitted_parameters(best_agreement_index,:);

        % Find the indicies in transverse agreement with this best fit
        % [current_agreement_indices, polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly_slow_method(points, fit_source_points, fittedCoefficients, transverse_tolerance, station_tolerance, (total_points_after_interpolating_points_in_domain), (-1));
        [current_agreement_indices,dist_btw_points_and_cubic_curve] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, fittedCoefficients, transverse_tolerance, (combos_paired(best_agreement_index,:)), (station_tolerance), (-1));

        % Accumulate results from best fits to prep for domain formation
        N_fits = N_fits+1;
        best_fit_cubic_poly_coefficients{N_fits} = fittedCoefficients;  %#ok<AGROW>
        best_fit_agreement_indicies{N_fits}      = current_agreement_indices;  %#ok<AGROW>
        best_fit_source_indicies{N_fits}         = combos_paired(best_agreement_index,:);  %#ok<AGROW>
        if 1 == slow_method
            best_fit_polygon_vertices{N_fits}        = polygon_vertices; %#ok<AGROW>   % comment this
        end
        % Remove all combos_paired sums that include this agreement
        for ith_agreement = 1:length(current_agreement_indices)
            % Grab current agreement
            agreement_to_remove = current_agreement_indices(ith_agreement);

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

domains = fcn_INTERNAL_filldomains(points, best_fit_cubic_poly_coefficients, best_fit_agreement_indicies, best_fit_source_indicies, total_points_after_interpolating_points_in_domain, transverse_tolerance);
if 1 == slow_method
    domains = fcn_INTERNAL_filldomains(points, best_fit_cubic_poly_coefficients, best_fit_agreement_indicies, best_fit_source_indicies, best_fit_polygon_vertices);
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

    % Plot the results in point space
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end        


    % Produce the sorted list, to create the Hough plot
    % [~,sorted_indicies] = sort(agreements,'ascend');

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

% fcn_INTERNAL_filldomains
% function domains = fcn_INTERNAL_filldomains(points, best_fit_cubic_poly_coefficients, best_fit_agreement_indicies, fit_source_indicies, best_fit_polygon_vertices)
function domains = fcn_INTERNAL_filldomains(points, best_fit_cubic_poly_coefficients, best_fit_agreement_indicies, fit_source_indicies, total_points_after_interpolating_points_in_domain, transverse_tolerance)

slow_method = 0; 
% How many domains are there?
N_domains = length(best_fit_agreement_indicies);

% Preallocate the domain structure, to include the unfitted points
domains{N_domains+1} = struct;

unfitted_indicies = ones(length(points(:,1)),1);

% Loop through domains
for ith_domain = 1:N_domains

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
        if isempty(domain_index)
            domain_index = nan;
        end
        domain_specific_indicies(ith_index) = domain_index;
    end
    best_fit_source_indicies    = domain_specific_indicies;


    % Calculate the best-fit domain box
    cubicPolyCoefficients   = best_fit_cubic_poly_coefficients{ith_domain};
    
    % Save the best-fit parameters
    best_fit_parameters = cubicPolyCoefficients;
    
    if length(points_in_domain(:,1)) > total_points_after_interpolating_points_in_domain                                        % uncomment this
        total_points_after_interpolating_points_in_domain = length(points_in_domain(:,1));
    end

    % Find the polygon vertices to show the domain box                                         % uncomment this
    polygon_vertices = fcn_INTERNAL_findPolygonVerticesForDomainBox(points_in_domain, total_points_after_interpolating_points_in_domain, best_fit_parameters, transverse_tolerance);
    domain_box = polygon_vertices; 
    domainShape = polyshape(domain_box(:,1),domain_box(:,2),'Simplify',false,'KeepCollinearPoints',true);

    % % % % Plot the input points
    % % % plot(points(:,1),points(:,2),'k.','MarkerSize',20);
    % % % 
    % % % % drawpolygon function is used to set the domain
    % % % drawpolygon(gca,"Position",polygon_vertices);
    
    if 1 == slow_method
        domainShape = fcn_geometry_domainBoxByType('cubic polynomial', best_fit_polygon_vertices{ith_domain},-1);   % comment this
    end
    %% Save results into the domain structure
    domains{ith_domain} = fcn_geometry_fillEmptyDomainStructure;

    domains{ith_domain}.best_fit_type = 'Hough cubic polynomial';
    domains{ith_domain}.points_in_domain         = points_in_domain;
    domains{ith_domain}.best_fit_source_indicies = best_fit_source_indicies;
    domains{ith_domain}.best_fit_domain_box      = domainShape;
    domains{ith_domain}.best_fit_parameters      = best_fit_parameters;

    %% Which points remain unfitted?
    unfitted_indicies(best_fit_agreement_indicies{ith_domain})=0;

end


% Save unfitted points into last structure
domains{end} = fcn_geometry_fillEmptyDomainStructure;
domains{end}.best_fit_type = 'unfitted';
domains{end}.points_in_domain = points(unfitted_indicies==1,:);

end % Ends fcn_INTERNAL_filldomains

function polygon_vertices = fcn_INTERNAL_findPolygonVerticesForDomainBox(points_in_domain, total_points_after_interpolating_points_in_domain, best_fit_parameters, transverse_tolerance)

% points_in_domain = points_in_domain(best_fit_source_indicies,:); 

% Interpolate x-coordinates of points within the station tolerance limit
x_interpolated_points_in_domain = fcn_INTERNAL_interpolateSourcePoints(points_in_domain, total_points_after_interpolating_points_in_domain);

% Find the y coordinates and slopes at each interpolated agreement point
[y_interpolated_points_in_domain, slopes_at_each_interpolated_points_in_domain] = fcn_INTERNAL_findSlopesAtEachPoint(x_interpolated_points_in_domain, best_fit_parameters);

% Find the unit orthogonal vectors at each interpolated agreement
% points based on the slope of the cubic polynomial
unit_orthogonal_vectors_interpolated_agreement_points = fcn_INTERNAL_findUnitOrthogonalVectors(slopes_at_each_interpolated_points_in_domain);

% Find interpolated agreement points matrix by putting x and y coordinates of
% interpolated agreement points together.
interpolated_points_in_domain = [x_interpolated_points_in_domain, y_interpolated_points_in_domain];

% Find the polygon vertices of the domain based on given transverse
% tolerance and station tolerance
[polygon_vertices,~,~, ~, ~] = fcn_INTERNAL_findDomainVertices(interpolated_points_in_domain, unit_orthogonal_vectors_interpolated_agreement_points, transverse_tolerance);

end

function interpolated_points_in_domain = fcn_INTERNAL_interpolateSourcePoints(points_in_domain, n_points)

% Sort the x_coordinates of points in domain
x_coordinates_points_in_domain = sort(points_in_domain(:,1),1);

% Generate the positions for the x_coordinates of points in domain 
indices = round(linspace(1, n_points, length(x_coordinates_points_in_domain)));

% Initialize the result vector with NaNs (or any placeholder)
interpolated_points_in_domain = NaN(1, n_points);

% Place the given numbers at the calculated indices
interpolated_points_in_domain(indices) = x_coordinates_points_in_domain;

% Find the indices of the NaNs (to interpolate)
nan_indices = find(isnan(interpolated_points_in_domain));

% Interpolate the NaNs
interpolated_points_in_domain(nan_indices) = interp1(indices, x_coordinates_points_in_domain, nan_indices, 'makima'); 

% The interpolated points are transposed
interpolated_points_in_domain = interpolated_points_in_domain';

end

function [y_interpolated_points_in_domain, slopes_at_each_interpolated_points_in_domain] = fcn_INTERNAL_findSlopesAtEachPoint(x_interpolated_points_in_domain, best_fit_parameters)

% Find the y coordinates of interpolated points in domain by substituting x
% coordinates of interpolated points in domain in cubic polynomial using
% "polyval"
y_interpolated_points_in_domain = polyval(best_fit_parameters, x_interpolated_points_in_domain);

%  Find the slopes (first derivative) of cubic polynomial at each
%  interpolated points in domain
slopes_at_each_interpolated_points_in_domain = 3*best_fit_parameters(1,1)*x_interpolated_points_in_domain.^2 + 2*best_fit_parameters(1,2)*x_interpolated_points_in_domain + best_fit_parameters(1,3);

% round off the slopes to a 4th decimal
slopes_at_each_interpolated_points_in_domain = round(slopes_at_each_interpolated_points_in_domain,4);

end

function unit_orthogonal_vectors = fcn_INTERNAL_findUnitOrthogonalVectors(slopes_at_each_interpolated_points_in_domain)

% Find angle of inclination to find the unit_tangent_vector
theta = atan(slopes_at_each_interpolated_points_in_domain); 
% theta = atan2(slopes_at_each_test_source_point, interpolated_source_points(:,1)); 

% Unit tangent vectors of all test source points
unit_tangent_vectors = [sin(theta), cos(theta)]; 

% Pre-allocate unit_orthogonal_vectors with zeros
unit_orthogonal_vectors = zeros(size(unit_tangent_vectors));

% If the slopes are negative, rotate the unit tangent vector in clockwise
% direction.
indices_slope_is_negative = slopes_at_each_interpolated_points_in_domain < 0;
unit_orthogonal_vectors(indices_slope_is_negative,:) = unit_tangent_vectors(indices_slope_is_negative,:)*[0 -1; 1 0];

% If the slopes are positive, rotate the unit tangent vector in anti
% clockwise direction
indices_slope_is_positive = slopes_at_each_interpolated_points_in_domain > 0;
unit_orthogonal_vectors(indices_slope_is_positive,:) = unit_tangent_vectors(indices_slope_is_positive,:)*[0 1; -1 0];

% If the slopes are zero, use tangent vectors to give transverse tolerance
% to generate the domain box
indices_slope_is_zero = slopes_at_each_interpolated_points_in_domain == 0;
unit_orthogonal_vectors(indices_slope_is_zero,:) = unit_tangent_vectors(indices_slope_is_zero,:);

end

function [polygon_vertices,upper_boundary_points,lower_boundary_points, y_values_upperboundary, y_values_lowerboundary]  = fcn_INTERNAL_findDomainVertices(interpolated_points_in_domain, unit_orthogonal_vectors, transverse_tolerance)

% Find the points on upper and lower boundaries
upper_boundary_points = interpolated_points_in_domain + unit_orthogonal_vectors*(2*transverse_tolerance);
lower_boundary_points = interpolated_points_in_domain - unit_orthogonal_vectors*(2*transverse_tolerance);

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
