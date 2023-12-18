function domains = fcn_geometry_HoughSegmentation(points, threshold_max_points, transverse_tolerance, station_tolerance, varargin)
% fcn_geometry_HoughSegmentation
% Given a set of points, attempts to segment the points into line segments,
% circles, arcs, etc. by using a Hough Transform methodology. 
% 
% Step 1: Points are checked against each type of fit using N choose K
% combinations of point indicies associated with the number of points
% needed for the minimum fit. For example, a line fit requires a
% minimum of 2 points, so if there are 5 points to fit, then the N choose
% combinations becomes (5 choose 2) which gives:
%
% Note, the nchoosek funciton in MATLAB produces these combination
% sequencies automatically.
%
% Format: 
% domains = fcn_geometry_HoughSegmentation(input_points, transverse_tolerance, station_tolerance, (fig_num))
%
% INPUTS:
%      points: a Nx2 vector where N is the number of points, but at least 2 rows. 
%
%      threshold_max_points: the number of points below which a match is no
%      longer performed.
%
%      transverse_tolerance: the orthogonal distance between the points and
%      a curve fit that indicate whether a point "belongs" to the fit
%      (if distance is less than or equal to the tolerance), or is
%      "outside" the fit (if distance is greater than the tolerance).
%
%      station_tolerance: the projection distance between the points in a
%      curve fit, along the direction of the fit, that indicate whether a
%      point "belongs" to the fit (if distance is less than or equal to the
%      tolerance), or is "outside" the fit (if distance is greater than the
%      tolerance).
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot the results.
%
% OUTPUTS:
%
%      domains: the fitted domains of the segments, ordered from the fit
%      that hast he most points included, to the domain that has the least.
%
% DEPENDENCIES:
%      nchoosek
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_HoughSegmentation
% for a full test suite.

% This function was written on 2023_12_15 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2023_12_14 
% -- wrote the code


flag_do_debug = 1; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking


if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 34838; %#ok<NASGU>
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

if flag_check_inputs == 1
    % Are there the right number of inputs?
    narginchk(4,5);
    
    % Check the points input to be length greater than or equal to 2
    fcn_DebugTools_checkInputsToFunctions(...
        points, '2column_of_numbers',[2 3]);
    
    % Check the threshold_max_points input is a positive single number
    fcn_DebugTools_checkInputsToFunctions(threshold_max_points, 'positive_1column_of_numbers',1);

    % Check the threshold_max_points input is a positive single number
    fcn_DebugTools_checkInputsToFunctions(transverse_tolerance, 'positive_1column_of_numbers',1);

    % Check the threshold_max_points input is a positive single number
    fcn_DebugTools_checkInputsToFunctions(station_tolerance, 'positive_1column_of_numbers',1);

end


% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if 5<= nargin
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

if flag_do_debug
    fig_debug = 8484;
else
    fig_debug = [];
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
 
if flag_do_debug
    figure(fig_debug);
    clf;
    hold on;
    grid on;

    plot(points(:,1),points(:,2),'k.','MarkerSize',20);
    title('Debugging figure for fcn_geometry_HoughSegmentation','Interpreter','none');
end

% Calculate the best agreements before starting the while loop
[best_agreement_count, best_fit_type, best_fit_parameters, best_fit_source_indicies, best_fit_associated_indicies, best_fit_domain_box] = ...
    fcn_INTERNAL_findBestFit(points, transverse_tolerance, station_tolerance, fig_debug);


% best_original_agreement = best_agreement_count; % For plotting
original_indicies = (1:length(points(:,1)))';
non_zero_indicies = find(original_indicies);

% Find domains based on removing peaks one at a time from agreements list
remaining_points = points;
remaining_indicies = ones(length(original_indicies),1);
domain_count = 0;
domains{1} = struct;

while best_agreement_count>threshold_max_points
    domain_count = domain_count+1;

    % Find the points that are inside the last domain that was found


    % Pull out the maximum agreement points and save them
    % Note: the points used for fitting keep changing (reduced each time by
    % the prior fit). So we have to keep track of the non-zero indicies
    % that are numbered relative to the ORIGINAL points, as the functions
    % return the indicies numbered relative to the input (reduced number)
    % points. 
    original_indicies_in_domain = non_zero_indicies(best_fit_associated_indicies);
    points_in_domain = points(original_indicies_in_domain,:);
    domains{domain_count}.points_in_domain = points_in_domain; %#ok<*AGROW>
    domains{domain_count}.best_fit_type = best_fit_type; 
    domains{domain_count}.best_fit_parameters = best_fit_parameters; 
    domains{domain_count}.best_fit_source_indicies = non_zero_indicies(best_fit_source_indicies); 
    domains{domain_count}.best_fit_domain_box = best_fit_domain_box;
    archive_of_remaining_points{domain_count}=remaining_points;% Used for plotting

    % Prune out the points in agreements up to this point, to leave only
    % ones NOT in agreement
    remaining_indicies(original_indicies_in_domain) = 0;
    non_zero_indicies = find(remaining_indicies);
    remaining_points = points(non_zero_indicies,:);

    % Recalculate the fit with remaining points
    [best_agreement_count, best_fit_type, best_fit_parameters, best_fit_source_indicies, best_fit_associated_indicies, best_fit_domain_box] = ...
        fcn_INTERNAL_findBestFit(remaining_points, transverse_tolerance, station_tolerance, fig_debug);

end

% The last domain is unfitted points
domain_count  = domain_count+1;
unfitted_points = remaining_points;
domains{domain_count}.points_in_domain = unfitted_points;
domains{domain_count}.best_fit_type = 'unfitted';
domains{domain_count}.best_fit_parameters = nan;
domains{domain_count}.best_fit_associated_indicies = nan;
domains{domain_count}.best_fit_domain_box = nan;
archive_of_remaining_points{domain_count}=unfitted_points;


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
    temp_h = figure(fig_num+1);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end        

    color_ordering = orderedcolors('gem12');
    N_colors = length(color_ordering(:,1));

    tiledlayout('flow')


    for ith_domain = 1:length(domains)
        nexttile;

        hold on;
        grid on;

        title(sprintf('Stage %.0d: Best fit type is %s', ith_domain, domains{ith_domain}.best_fit_type),'Interpreter','none');

        % Plot the remaining points
        remaining_points = archive_of_remaining_points{ith_domain};
        plot(points(:,1),points(:,2),'.','MarkerSize',30,'Color',[0.7 0.7 0.7]);
        plot(remaining_points(:,1),remaining_points(:,2),'k.','MarkerSize',20);

        % Plot the fitted points
        points_in_domain = domains{ith_domain}.points_in_domain;
        plot(points_in_domain(:,1),points_in_domain(:,2),'.','MarkerSize',15,'Color',color_ordering(mod(ith_domain,N_colors)+1,:));

        % Plot the fits?
        if ~any(isnan(domains{ith_domain}.best_fit_parameters))
            switch domains{ith_domain}.best_fit_type
                case 'line'
                    % Plot the best-fit line segment
                    plot(domains{ith_domain}.best_fit_parameters(:,1),domains{ith_domain}.best_fit_parameters(:,2),'.-','Linewidth',1,'MarkerSize',15,'Color',color_ordering(mod(ith_domain,N_colors)+1,:));

                case 'circle'

                case 'arc'

                otherwise
                    error('Unknown fit type detected - unable to continue!');
            end
        end

        % Make axis slightly larger?
        if flag_rescale_axis
            if ith_domain==1
                temp = axis;
                axis_range_x = temp(2)-temp(1);
                axis_range_y = temp(4)-temp(3);
                percent_larger = 0.3;
                new_axis = [temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y];
            end
            axis(new_axis);
        end


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

%% fcn_INTERNAL_findBestFit
function [best_agreement_count, best_fit_type, best_fit_parameters, best_fit_source_indicies, best_fit_associated_indicies, best_fit_domain_box] = ...
    fcn_INTERNAL_findBestFit(input_points,transverse_tolerance, station_tolerance, fig_num)

% Calculate the best agreements between different fit types, returning the
% best fitting type.

fit_types = {'line','circle','arc'};
best_agreement_counts = zeros(length(fit_types),1);

% Initialize variables
fit_parameters{length(fit_types)} = 0;
fit_associated_point_indicies{length(fit_types)} = 0;
fit_nchoosek_index = zeros(length(fit_types),1);
fit_source_index_list{length(fit_types)} = 0;

% Search through all the fitting types to find the best fits.
% Each fit must do several things:
% 
% Step 1: Given a minimum model order - how many points are needed for the
% minimum fit? - finds all possible model fits.
%
% Step 2: Given all possible model fits, finds which points are associated
% with each fit given tolerance fators. These point totals represent the
% "votes" for that particular fit.
%
%
% The above process is usually repeated among many different types of fits
% to the points: lines, arcs, spirals, etc.
%
% (IN SEPARATE CODE)
% Step 3: Find the best permutation of a fit, namely the one that has the
% highest votes. Using these, find the regression fit coefficients
% and domain of the best fit. 
% 
% Step 4: Given the domain, find the points that are members of that
% best-fit permutation.
%

% Steps 1 and 2 are done within this for-loop
for fit_index = 1:length(fit_types)
    fit_type = fit_types{fit_index};

    switch fit_type
        case 'line'
            % Check line fitting - minimum model order is 2 points
            [fitted_parameters, agreement_indicies] = fcn_geometry_fitHoughLine(input_points, transverse_tolerance, station_tolerance);
            indicies = nchoosek(1:length(input_points(:,1)),2); % Line fit needs 2 points

        case 'circle'
            % Check circle fitting - minimum model order is 3 points
            % [fitted_parameters, agreement_indicies] = fcn_geometry_fitHoughCircle(input_points, tolerance, fig_num);
            fitted_parameters  = [0 0 0];
            agreement_indicies = 1;
            indicies = nchoosek(1:length(input_points(:,1)),3); % Circle fit needs 3 points

        case 'arc'
            % Check arc fitting - minimum model order is 3 points
            % [fitted_parameters, agreement_indicies] = fcn_geometry_fitHoughArc(input_points, tolerance, fig_num);
            fitted_parameters  = [0 0 0 0 0];
            agreement_indicies = 1;
            indicies = nchoosek(1:length(input_points(:,1)),3); % Arc fit needs 3 points
        otherwise
            error('Unknown fit type detected - unable to continue!');
    end
    agreements = sum(agreement_indicies,2);
    [best_agreement_counts(fit_index), index_of_best_fit] = max(agreements);
    fit_parameters{fit_index}     = fitted_parameters(index_of_best_fit,:);
    fit_associated_point_indicies{fit_index} = agreement_indicies(index_of_best_fit,:);
    fit_nchoosek_index(fit_index) = index_of_best_fit;
    fit_source_index_list{fit_index} = indicies(index_of_best_fit,:);
end

% Save the best fit among all the types
[best_agreement_count, fitting_type_number] = max(best_agreement_counts);
best_fit_type = fit_types{fitting_type_number};
best_fit_parameters = fit_parameters{fitting_type_number};
best_fit_source_indicies = fit_source_index_list{fitting_type_number};
best_fit_associated_indicies = find(fit_associated_point_indicies{fitting_type_number});

% Steps 3 and 4 follow within this switch case
% Find regression fit using best agreement and domain
best_fit_domain_box  = nan;
source_points = input_points(best_fit_source_indicies,:);
associated_points_in_domain = input_points(best_fit_associated_indicies,:);
switch best_fit_type
    case 'line'
        % Check line fitting
        [best_fit_parameters, best_fit_domain_box] = fcn_geometry_calcLinearRegressionFromHoughFit(source_points,associated_points_in_domain, fig_num);

        % The following 3 lines might be able to be moved outside the switch case
        % as they are likely going to be the same for all cases.
        domainPolyShape = polyshape(best_fit_domain_box(:,1),best_fit_domain_box(:,2),'Simplify',false,'KeepCollinearPoints',true);
        IndiciesOfPointsInDomain = isinterior(domainPolyShape,input_points);
        best_fit_associated_indicies = find(IndiciesOfPointsInDomain);
    case 'circle'
        % Fill this in
    case 'arc'
        % Fill this in
    otherwise
        error('Unknown fit type detected - unable to continue!');
end


end % Ends fcn_INTERNAL_findBestFit



