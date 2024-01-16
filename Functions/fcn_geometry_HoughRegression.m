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
% domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, (fig_num))
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
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      domains: the fitted domains of the segments, ordered from the fit
%      that hast he most points included, to the domain that has the least.
%
% DEPENDENCIES:
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_fitLinearRegressionFromHoughFit 
%      fcn_geometry_fitCircleRegressionFromHoughFit
%      fcn_geometry_fitCircleRegressionFromHoughFit
%      fcn_geometry_plotCircle
%      fcn_geometry_plotArc
%      fcn_geometry_fitHoughLine
%      fcn_geometry_fitHoughCircle
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_HoughSegmentation
% for a full test suite.
%
% This function was written on 2023_12_15 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2023_12_14 
% -- wrote the code
% 2024_01_01 
% -- modified to use updated versions of fitHoughLine code set
% 2024_01_02 
% -- line segment versions now working
% 2024_01_03 - S. Brennan
% -- added fast mode option
% -- added environmental variable options
% 2024_01_09- S. Brennan
% -- added circle fitting
% 2024_01_10- S. Brennan
% -- added arc fitting
% 2024_01_15 - S. Brennan
% -- typo fix on domain field

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
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (5<= nargin) && (0==flag_max_speed)
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






% Search through all the fitting types to find the best fits. The fits are
% ordered from simplest (and fastest) first, then to more complex. This is
% necessary not only for speed, but also because of degeneracy in the
% fitting process itself.
%
% For example, all lines are circles with infinite radius, all circles are
% arcs that extend up to 2*pi, all arcs are spirals with a radial slope
% equal to zero, etc. In other words, to allow fits of more advanced and
% complex geometric representations before trying simpler ones first, the
% simpler forms would not be used.%
%
% Each fit must do several things:
% 
% Step 1: Given a minimum model order - how many points are needed for the
% minimum fit? Given these points, find all possible index permutations
% that can form a complete model. For example, a line fit requires 2 points
% minimum. If there are N points, then there are N choose 2 permutations
% possible assuming the line is not ordered.
%
% Step 2: Given all possible model fits, find the fit and find also which
% points are associated with each fit given tolerance fators. These point
% totals represent the "votes" for that particular fit.
%
% Step 3: Find the best permutation of a fit, namely the one that has the
% highest votes among all the fits of that type - the best line fit among
% all line fits, for example. Using these, find the regression fit
% coefficients and domain of the best fit.
% 
% Step 4: Given the domain, find the points that are members of that
% best-fit permutation.
%
% The above process is repeated until no more fits are possible for a given
% fit type, starting with the simplest fit type first (lines). After all
% lines are fitted, this proceeds further among many different types of
% fits to the points: lines, arcs, spirals, etc.

% best_original_agreement = best_agreement_count; % For plotting
original_indicies = (1:length(points(:,1)))';
non_zero_indicies = find(original_indicies);

% Find domains based on removing peaks one at a time from agreements list
remaining_points = points;
remaining_indicies = ones(length(original_indicies),1);
domain_count = 0;
domains{1} = struct;

% Define the fit types to search
fit_types = {'line'}; % ,'circle','arc'};
archive_of_remaining_points{1:length(fit_types)} = [];

fprintf(1,'\n\n');

for fit_index = 1:length(fit_types)
    fit_type_to_search = fit_types{fit_index};
    fprintf(1,'Attempting fits of type: %s\n',fit_type_to_search);

    % Calculate the best agreements of this type before starting the while loop
    if length(remaining_points(:,1))>threshold_max_points
        % Perform Hough fit
        Hough_domains = fcn_INTERNAL_findBestFit(remaining_points, transverse_tolerance, station_tolerance, fit_type_to_search, threshold_max_points);
        N_Hough_domains = length(Hough_domains)-1;
        fprintf(1,'Found: %.0d \n',length(Hough_domains)-1);
        
        if N_Hough_domains > 0

            % Perform regression fit on best-fit Hough values
            regression_domains = cell(N_Hough_domains,1);
            for ith_domain = 1:N_Hough_domains
                 regression_domain = fcn_INTERNAL_regressionFitByType(Hough_domains{ith_domain}, fit_type_to_search); 
                 regression_domains{ith_domain} = regression_domain;
            end
        end

        % Save the results
        for ith_domain = 1:N_Hough_domains
            domain_count = domain_count+1;
            domains{domain_count} = Hough_domains{ith_domain};            %#ok<AGROW>
        end
        for ith_domain = 1:N_Hough_domains
            domain_count = domain_count+1;
            domains{domain_count} = regression_domains{ith_domain}; %#ok<AGROW>
        end

        % Update the unfitted points
        remaining_points = Hough_domains{end}.points_in_domain;
        archive_of_remaining_points{fit_index} = remaining_points;
    end
        
end % Ends looping through fitting types

% The last domain is always going to be unfitted points
domain_count  = domain_count+1;
unfitted_points = remaining_points;
domains{domain_count}.points_in_domain = unfitted_points;
domains{domain_count}.best_fit_type = 'unfitted';
domains{domain_count}.best_fit_parameters = nan;
domains{domain_count}.best_fit_source_indicies = nan;
domains{domain_count}.best_fit_domain_box = nan;
archive_of_remaining_points{fit_index}=unfitted_points;


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

    try
        color_ordering = orderedcolors('gem12');
    catch
        color_ordering = colororder;
    end
    
    N_colors = length(color_ordering(:,1));

    tiledlayout('flow')

    % Plot the starting points in the first plot, save the axis
    nexttile;

    hold on;
    grid on;
    axis equal;

    title(sprintf('Starting points'),'Interpreter','none');
    plot(points(:,1),points(:,2),'k.','MarkerSize',20);

    % Make axis slightly larger? And since this is the first one, save the
    % axis limits.
    if flag_rescale_axis
        temp = axis;
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        new_axis = [temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y];
        axis(new_axis);
    end


    % Plot the fits
    for ith_domain = 1:length(domains)
        nexttile;

        hold on;
        grid on;
        axis equal;

        title(sprintf('Stage %.0d: Best fit type is %s', ith_domain, domains{ith_domain}.best_fit_type),'Interpreter','none');

        % Plot all the points in very light grey
        plot(points(:,1),points(:,2),'.','MarkerSize',30,'Color',[0.8 0.8 0.8]);

        fcn_geometry_plotFitDomains(domains{ith_domain}, fig_num);


        % plot(remaining_points(:,1),remaining_points(:,2),'k.','MarkerSize',20);
        % 
        % % Plot the fitted points
        % current_color = color_ordering(mod(ith_domain,N_colors)+1,:);
        % points_in_domain = domains{ith_domain}.points_in_domain;
        % plot(points_in_domain(:,1),points_in_domain(:,2),'.','MarkerSize',15,'Color',current_color);
        % 
        % % Plot the fits?
        % if ~any(isnan(domains{ith_domain}.best_fit_parameters))
        %     switch domains{ith_domain}.best_fit_type
        %         case 'line'
        %             % Plot the best-fit line segment
        %             plot(domains{ith_domain}.best_fit_parameters(:,1),domains{ith_domain}.best_fit_parameters(:,2),'.-','Linewidth',1,'MarkerSize',15,'Color',current_color);
        % 
        %             % Plot the domain box
        %             domain_box = domains{ith_domain}.best_fit_domain_box;
        %             domainShape = polyshape(domain_box(:,1),domain_box(:,2),'Simplify',false,'KeepCollinearPoints',true);
        %             plot(domainShape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',1,'EdgeAlpha',0);
        % 
        %         case 'circle'
        %             % Plot the best-fit circle  
        %             circleCenter = domains{ith_domain}.best_fit_parameters(1,1:2);
        %             circleRadius = domains{ith_domain}.best_fit_parameters(1,3);
        %             fcn_geometry_plotCircle(circleCenter, circleRadius, current_color,fig_num)
        %             plot(circleCenter(1,1),circleCenter(1,2),'+','MarkerSize',30,'Color',current_color);
        % 
        %             % Plot the domain box
        %             domain_box = domains{ith_domain}.best_fit_domain_box;
        %             domainShape = polyshape(domain_box(:,1),domain_box(:,2),'Simplify',false,'KeepCollinearPoints',true);
        %             plot(domainShape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',1,'EdgeAlpha',0);
        % 
        %         case 'arc'
        %             % Plot the best-fit arc  
        %             circleCenter = domains{ith_domain}.best_fit_parameters(1,1:2);
        %             circleRadius = domains{ith_domain}.best_fit_parameters(1,3);
        %             start_angle_in_radians = domains{ith_domain}.best_fit_parameters(1,4);
        %             end_angle_in_radians = domains{ith_domain}.best_fit_parameters(1,5);
        %             fcn_geometry_plotArc(circleCenter, circleRadius, start_angle_in_radians, end_angle_in_radians, current_color);
        %             plot(circleCenter(1,1),circleCenter(1,2),'+','MarkerSize',30,'Color',current_color);
        % 
        %             % Plot the domain box
        %             domain_box = domains{ith_domain}.best_fit_domain_box;
        %             domainShape = polyshape(domain_box(:,1),domain_box(:,2),'Simplify',false,'KeepCollinearPoints',true);
        %             plot(domainShape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',1,'EdgeAlpha',0);
        % 
        %         otherwise
        %             error('Unknown fit type detected - unable to continue!');
        %     end
        % end

        % Make axis slightly larger?
        if flag_rescale_axis
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
function Hough_domains = ...
    fcn_INTERNAL_findBestFit(input_points,transverse_tolerance, station_tolerance, fit_type_to_search, threshold_max_points)

% Calculate the best agreements between different fit types, returning the
% best fitting type.

% %% Start of changes moving fit type out of this function
% fit_types = {'line','circle'}; %,'circle','arc'};
% best_agreement_counts = zeros(length(fit_types),1);
% 
% % Initialize variables
% fit_parameters{length(fit_types)} = 0;
% fit_associated_point_indicies{length(fit_types)} = 0;
% fit_source_index_list{length(fit_types)} = 0;
% 
% % Search through all the fitting types to find the best fits.
% % Each fit must do several things:
% % 
% % Step 1: Given a minimum model order - how many points are needed for the
% % minimum fit? - finds all possible model fits.
% %
% % Step 2: Given all possible model fits, finds which points are associated
% % with each fit given tolerance fators. These point totals represent the
% % "votes" for that particular fit.
% %
% %
% % The above process is usually repeated among many different types of fits
% % to the points: lines, arcs, spirals, etc.
% %
% % (IN SEPARATE CODE)
% % Step 3: Find the best permutation of a fit, namely the one that has the
% % highest votes. Using these, find the regression fit coefficients
% % and domain of the best fit. 
% % 
% % Step 4: Given the domain, find the points that are members of that
% % best-fit permutation.
% %
% 
% % Steps 1 and 2 are done within this for-loop
% for fit_index = 1:length(fit_types)
%     fit_type = fit_types{fit_index};
% % End of changes moving fit type out of this function

switch fit_type_to_search
    case 'line'
        % Check line fitting - minimum model order is 2 points
        Hough_domains = ...
            fcn_geometry_fitHoughLine(input_points, transverse_tolerance, station_tolerance, threshold_max_points, -1);

    case 'circle'
        flag_force_circle_fit = 1;
        expected_radii_range = [2 8];
        flag_use_permutations = [];

        % Check circle fitting - minimum model order is 3 points
        [best_fit_parameters, best_fit_source_indicies, best_fit_associated_indicies] = ...
            fcn_geometry_fitHoughCircle(input_points, transverse_tolerance, (station_tolerance), (flag_force_circle_fit), (expected_radii_range), (flag_use_permutations), (-1));

    case 'arc'
        % If the station_tolerance is empty, then arc fitting should not do
        % anything as this will force circle fits instead. We do not want
        % to fit arcs to possible circles.
        if ~isempty(station_tolerance)
            flag_force_circle_fit = 0;
            expected_radii_range = [2 8];
            flag_use_permutations = [];

            % Check circle fitting - minimum model order is 3 points
            [best_fit_parameters, best_fit_source_indicies, best_fit_associated_indicies] = ...
                fcn_geometry_fitHoughCircle(input_points, transverse_tolerance, (station_tolerance), (flag_force_circle_fit), (expected_radii_range), (flag_use_permutations), (-1));
        else
            best_fit_parameters  = [0 0 0];
            best_fit_source_indicies = 1;
            best_fit_associated_indicies = 1;
        end
    otherwise
        error('Unknown fit type detected - unable to continue!');
end
% [best_agreement_count,...
%     best_fit_parameters, ...
%     best_fit_source_indicies, ...
%     best_fit_associated_indicies] 
% 
% % Count how many points agreed with the fit
% best_agreement_count = length(best_fit_associated_indicies);

% COMMENTED OUT WHEN MOVED RESULTS TO OUTER LOOP
% % Save results to an array for each fit type
% best_agreement_counts(fit_index) = length(best_agreement_indicies);
% fit_parameters{fit_index}     = fitted_parameters;
% fit_associated_point_indicies{fit_index} = best_agreement_indicies;
% fit_source_index_list{fit_index} = best_fit_source_indicies;
% 
% % Save the best fit among all the fitting types
% [best_agreement_count, fitting_type_number] = max(best_agreement_counts);
% best_fit_type = fit_types{fitting_type_number};
% best_fit_parameters = fit_parameters{fitting_type_number};
% best_fit_source_indicies = fit_source_index_list{fitting_type_number};
% best_fit_associated_indicies = fit_associated_point_indicies{fitting_type_number};

end % Ends fcn_INTERNAL_findBestFit

function regression_domain = fcn_INTERNAL_regressionFitByType(Hough_domain, best_fit_type)

% Find regression fit using best agreement and domain
switch best_fit_type
    case 'line'
        % Check line fitting
        regression_domain = fcn_geometry_fitLinearRegressionFromHoughFit(Hough_domain, -1);

    case 'circle'
        % Check circle fitting
        [best_fit_parameters, best_fit_domain_box]  = ...
            fcn_geometry_fitCircleRegressionFromHoughFit(source_points, associated_points_in_domain, -1);

    case 'arc'
        % Check arc fitting
        [best_fit_parameters, best_fit_domain_box]  = ...
            fcn_geometry_fitArcRegressionFromHoughFit(source_points, associated_points_in_domain, -1);

    otherwise
        error('Unknown fit type detected - unable to continue!');
end

end % Ends fcn_INTERNAL_regressionFitByType





