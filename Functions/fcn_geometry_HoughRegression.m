function regression_domains = fcn_geometry_HoughRegression(Hough_domains, varargin)
% fcn_geometry_HoughRegression
% Given a set of Hough fits of points in regions, performs regression
% fitting on each domain.
%
% Format: 
% domains = fcn_geometry_HoughRegression(Hough_domains, (transverse_tolerance), (fig_num))
%
% INPUTS:
%      Hough_domain: a structure that records details of the domain of
%      fitting. See fcn_geometry_fillEmptyDomainStructure for details.
%
%      (OPTIONAL INPUTS)
% 
%      transverse_tolerance: the distance from the curve fit, in the
%      transverse direction, to project in both the positive and negative
%      directions to produce the best_fit_domain_box. If left empty,
%      defaults to 2 standard deviations to thus give a box that is +/- 2
%      sigma.
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
% See the script: script_test_fcn_geometry_HoughRegression
% for a full test suite.
%
% This function was written on 2024_01_18 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2024_01_18
% -- wrote the code
% 2024_05_17 - Aneesh Batchu
% -- Added a case for 'Hough cubic polynomial'

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

        % % Check the points input to be length greater than or equal to 2
        % fcn_DebugTools_checkInputsToFunctions(...
        %     points, '2column_of_numbers',[2 3]);
        %
        % % Check the threshold_max_points input is a positive single number
        % fcn_DebugTools_checkInputsToFunctions(threshold_max_points, 'positive_1column_of_numbers',1);
        %
        % % Check the threshold_max_points input is a positive single number
        % fcn_DebugTools_checkInputsToFunctions(transverse_tolerance, 'positive_1column_of_numbers',1);
        %
        % % Check the threshold_max_points input is a positive single number
        % fcn_DebugTools_checkInputsToFunctions(station_tolerance, 'positive_1column_of_numbers',1);

    end
end

% Does user want to specify best_fit_domain_box_projection_distance?
transverse_tolerance = [];
if (2<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        transverse_tolerance = temp;
    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (3<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

if flag_do_debug
    fig_debug = 8484; %#ok<NASGU>
else
    fig_debug = []; %#ok<NASGU>
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
N_Hough_domains = length(Hough_domains)-1;

if N_Hough_domains > 0
    % Perform regression fit on best-fit Hough values
    regression_domains = cell(N_Hough_domains+1,1);
    for ith_domain = 1:N_Hough_domains
        regression_domain = fcn_INTERNAL_regressionFitByType(Hough_domains{ith_domain}, transverse_tolerance);
        regression_domains{ith_domain} = regression_domain;
    end
end


% The last domain is always going to be unfitted points
domain_count  = N_Hough_domains+1;
regression_domains{domain_count} = Hough_domains{end};


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

    tiledlayout('flow')

    % Plot the starting points in the first plot, save the axis
    nexttile;

    hold on;
    grid on;
    axis equal;
    title(sprintf('Starting points'),'Interpreter','none');

    % Gather the points
    points = [];    
    for ith_domain = 1:length(regression_domains)
        points = [points; regression_domains{ith_domain}.points_in_domain]; %#ok<AGROW>
    end
    
    % Plot the points
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
    for ith_domain = 1:length(regression_domains)
        nexttile;

        hold on;
        grid on;
        axis equal;

        title(sprintf('Stage %.0d: Best fit type is %s', ith_domain, regression_domains{ith_domain}.best_fit_type),'Interpreter','none');

        % Plot all the points in very light grey
        plot(points(:,1),points(:,2),'.','MarkerSize',30,'Color',[0.8 0.8 0.8]);

        fcn_geometry_plotFitDomains(regression_domains{ith_domain}, fig_num);


       
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

function regression_domain = fcn_INTERNAL_regressionFitByType(Hough_domain, transverse_tolerance)

best_fit_type = Hough_domain.best_fit_type;

% Find regression fit using best agreement and domain
switch best_fit_type
    case {'Hough line','Hough segment'}
        % Check line fitting
        regression_domain = fcn_geometry_fitLinearRegressionFromHoughFit(Hough_domain, transverse_tolerance, -1);

    case {'Hough circle','Hough arc'}
        % Check circle fitting
        regression_domain = ...
            fcn_geometry_fitArcRegressionFromHoughFit(Hough_domain, transverse_tolerance, -1);
    case {'Hough cubic polynomial'}
        % Check cubic polynomial fitting
        % regression_domain = fcn_INTERNAL_fitCubicPolynomial(Hough_domain, transverse_tolerance);
        regression_domain = fcn_geometry_polyFitCubicPolyFromHoughFit(Hough_domain, transverse_tolerance, -1);

    otherwise
        error('Unknown fit type detected - unable to continue!: %s', best_fit_type);
end

end % Ends fcn_INTERNAL_regressionFitByType
