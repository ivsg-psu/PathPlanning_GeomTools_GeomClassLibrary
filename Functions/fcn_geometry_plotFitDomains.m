function fcn_geometry_plotFitDomains(domains, varargin)
% fcn_geometry_plotFitDomains
% Given a set of domains obtained by Hough or regression fitting, plots
% these. An optional figure number can be given.
%
% Format:
% fcn_geometry_plotFitDomains(domains, (fig_num))
%
% INPUTS:
%      domains: a structure that includes the domains to be plotted. See
%      fcn_geometry_HoughSegmentation for example outputs
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      (none)
%
% DEPENDENCIES:
%      (none)
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_plotFitDomains
% for a full test suite.
%
% This function was written on 2024_01_12 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_01_12 - S. Brennan
% -- wrote the code
% 2024_01_16
% -- fixed bug with line/segment plotting


%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==2 && isequal(varargin{end},-1))
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
        narginchk(1,2);
    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 1;
if (2<= nargin) && (0==flag_max_speed)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

if isempty(fig_num)
    fig_num = gcf; % Default is to grab the current figure
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

% Nothing to do in here!

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

    % Get the color ordering?
    try
        color_ordering = orderedcolors('gem12');
    catch
        color_ordering = colororder;
    end

    N_colors = length(color_ordering(:,1));

    hold on;
    grid on;
    axis equal;

    %% Plot all the points and domains
    for ith_domain = 1:length(domains)
        % Get current color
        current_color = color_ordering(mod(ith_domain,N_colors)+1,:);

        % Convert domain to one structure, to make things easy
        if iscell(domains)
            domain_to_plot = domains{ith_domain};
        else
            domain_to_plot = domains;
        end

        if ~any(isnan(domain_to_plot.best_fit_parameters))
            % Plot the points
            domainPoints = domain_to_plot.points_in_domain;
            plot(domainPoints(:,1), domainPoints(:,2), '.','MarkerSize',20, 'Color',current_color);

            % Plot the domain shape
            domainShape = domain_to_plot.best_fit_domain_box;
            plot(domainShape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',1,'EdgeAlpha',0);

        else % These are 'unfitted' models, typically

            % Plot the unfitted points in grey
            domainPoints = domain_to_plot.points_in_domain;
            plot(domainPoints(:,1), domainPoints(:,2), '.','MarkerSize',20, 'Color',[0.5 0.5 0.5]);
        end
    end

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

    % Get the maximum plotting distance - we'll need this on some of the
    % fits (lines)
    current_axis = axis;
    axis_range_x = current_axis(2)-current_axis(1);
    axis_range_y = current_axis(4)-current_axis(3);
    max_distance = sum([axis_range_x axis_range_y].^2,2).^0.5;

    % Plot the fits
    for ith_domain = 1:length(domains)
        % Get current color
        current_color = color_ordering(mod(ith_domain,N_colors)+1,:);

        % Convert domain to one structure, to make things easy
        if iscell(domains)
            domain_to_plot = domains{ith_domain};
        else
            domain_to_plot = domains;
        end

        if ~any(isnan(domain_to_plot.best_fit_parameters))

            % Plot the model fit
            model_fit  = domain_to_plot.best_fit_parameters;
            switch domain_to_plot.best_fit_type

                case {'Hough line','Vector regression line fit'}
                    % Plot the best-fit line 
                    %             [unit_projection_vector_x,
                    %              unit_projection_vector_y,
                    %              base_point_x,
                    %              base_point_y,
                    %             ]
                    unit_projection_vector = model_fit(1,1:2);
                    base_point = model_fit(1,3:4);
                    low_station  = -max_distance;
                    high_station = max_distance;
                    line_segment = [base_point + unit_projection_vector*low_station; base_point + unit_projection_vector*high_station];
                    plot(line_segment(:,1),line_segment(:,2),'.-','Linewidth',3,'MarkerSize',15,'Color',current_color);

                case {'Hough segment','Vector regression segment fit'}
                    % Plot the best-fit segment                    
                    %             [unit_projection_vector_x,
                    %              unit_projection_vector_y,
                    %              base_point_x,
                    %              base_point_y,
                    %              station_distance_min,
                    %              station_distance_max,
                    %             ]
                    unit_projection_vector = model_fit(1,1:2);
                    base_point = model_fit(1,3:4);
                    low_station  = model_fit(1,5);
                    high_station = model_fit(1,6);
                    line_segment = [base_point + unit_projection_vector*low_station; base_point + unit_projection_vector*high_station];
                    plot(line_segment(:,1),line_segment(:,2),'.-','Linewidth',3,'MarkerSize',15,'Color',current_color);

                case 'circle'
                    % Plot the best-fit circle
                    circleCenter = domain_to_plot.best_fit_parameters(1,1:2);
                    circleRadius = domain_to_plot.best_fit_parameters(1,3);
                    fcn_geometry_plotCircle(circleCenter, circleRadius, current_color,fig_num)
                    plot(circleCenter(1,1),circleCenter(1,2),'+','MarkerSize',30,'Color',current_color);

                    % Plot the domain box
                    domain_box = domain_to_plot.best_fit_domain_box;
                    domainShape = polyshape(domain_box(:,1),domain_box(:,2),'Simplify',false,'KeepCollinearPoints',true);
                    plot(domainShape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',1,'EdgeAlpha',0);

                case 'arc'
                    % Plot the best-fit arc
                    circleCenter = domain_to_plot.best_fit_parameters(1,1:2);
                    circleRadius = domain_to_plot.best_fit_parameters(1,3);
                    start_angle_in_radians = domain_to_plot.best_fit_parameters(1,4);
                    end_angle_in_radians = domain_to_plot.best_fit_parameters(1,5);
                    fcn_geometry_plotArc(circleCenter, circleRadius, start_angle_in_radians, end_angle_in_radians, current_color,fig_num)  

                    plot(circleCenter(1,1),circleCenter(1,2),'+','MarkerSize',30,'Color',current_color);

                    % Plot the domain box
                    domain_box = domain_to_plot.best_fit_domain_box;
                    domainShape = polyshape(domain_box(:,1),domain_box(:,2),'Simplify',false,'KeepCollinearPoints',true);
                    plot(domainShape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',1,'EdgeAlpha',0);

                otherwise
                    error('Unknown fit type detected: %s - unable to continue!', domain_to_plot.best_fit_type);
            end
        end

        % Force axis back to what it was before
        axis(current_axis);

        if strcmp(domain_to_plot.best_fit_type,'unfitted')
            % Plot the points in grey
            domainPoints = domain_to_plot.points_in_domain;
            plot(domainPoints(:,1), domainPoints(:,2), '.','MarkerSize',20, 'Color',[0.5 0.5 0.5]);
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