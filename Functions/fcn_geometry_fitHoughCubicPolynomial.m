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
total_points_including_source_points = 20;
if (6<= nargin)
    temp = varargin{4};
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
    % ith_combo
    % conbo_pair = combos_paired(ith_combo,:)
    % [agreement_indices,~] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, combos_paired(ith_combo,1), station_tolerance,(-1));
    [agreement_indices,~] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, test_source_points, fittedParameters, transverse_tolerance, (station_tolerance),(total_points_including_source_points),(-1));

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
best_fit_polygon_vertices        = {};

% Find best agreements
if 1 == flag_find_only_best_agreement
    N_fits = 0;
    remaining_agreements = agreements;

    [~, best_agreement_index] = max(remaining_agreements);

    % Best fit source points
    only_best_fit_source_points = points(combos_paired(best_agreement_index,:),:);

    % Find the best fit parameters
    only_best_fit_fittedCoefficients  = fitted_parameters(best_agreement_index,:);
    
    % Find the indicies in transverse agreement with this best fit
    % [only_best_agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, only_best_fit_source_points, only_best_fit_fittedCoefficients, transverse_tolerance, combos_paired(best_agreement_index,1), station_tolerance, (-1));
    [only_best_agreement_indices,polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, only_best_fit_source_points, only_best_fit_fittedCoefficients, transverse_tolerance, (station_tolerance), (total_points_including_source_points), (-1));

    % Accumulate results from best fits to prep for domain formation
    N_fits = N_fits+1;
    best_fit_cubic_poly_coefficients{N_fits}  = only_best_fit_fittedCoefficients;
    best_fit_agreement_indicies{N_fits}     = only_best_agreement_indices;
    best_fit_source_indicies{N_fits}        = combos_paired(best_agreement_index,:);
    best_fit_polygon_vertices{N_fits}        = polygon_vertices;
else
    N_fits = 0;
    remaining_agreements = agreements;

    [best_agreement_count, best_agreement_index] = max(remaining_agreements);

    while best_agreement_count >= points_required_for_agreement

        % test source points
        fit_source_points = points(combos_paired(best_agreement_index,:),:);

        % Find the best fit parameters
        fittedCoefficients  = fitted_parameters(best_agreement_index,:);

        % Find the indicies in transverse agreement with this best fit
        % [current_agreement_indices, polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, fit_source_points, fittedCoefficients, transverse_tolerance, combos_paired(best_agreement_index,1), station_tolerance, (-1));
        [current_agreement_indices, polygon_vertices] = fcn_geometry_findAgreementsOfPointsToCubicPoly(points, fit_source_points, fittedCoefficients, transverse_tolerance, (station_tolerance), (total_points_including_source_points), (-1));

        % Accumulate results from best fits to prep for domain formation
        N_fits = N_fits+1;
        best_fit_cubic_poly_coefficients{N_fits} = fittedCoefficients;  %#ok<AGROW>
        best_fit_agreement_indicies{N_fits}      = current_agreement_indices;  %#ok<AGROW>
        best_fit_source_indicies{N_fits}         = combos_paired(best_agreement_index,:);  %#ok<AGROW>
        best_fit_polygon_vertices{N_fits}        = polygon_vertices; %#ok<AGROW>

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

domains = fcn_INTERNAL_filldomains(points, best_fit_cubic_poly_coefficients, best_fit_agreement_indicies, best_fit_source_indicies, best_fit_polygon_vertices);

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
function domains = fcn_INTERNAL_filldomains(points, best_fit_cubic_poly_coefficients, best_fit_agreement_indicies, fit_source_indicies, best_fit_polygon_vertices)

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
    
    domainShape = fcn_geometry_domainBoxByType('cubic polynomial', best_fit_polygon_vertices{ith_domain},-1);

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
