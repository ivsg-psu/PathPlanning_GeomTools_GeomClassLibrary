function [fitSequence_points, fitSequence_shapes, fitSequence_endIndicies, fitSequence_parameters, fitSequence_bestFitType] = fcn_geometry_fitSequentialArcsC2(points_to_fit, varargin)
%% fcn_geometry_fitSequentialArcsC2
% This function is a repeat of fcn_geometry_fitSequentialArcs, but forcing
% C2 continuity between arcs in sequence. See the original function for
% more details.
% 
% Format: 
% [fitSequence_points, fitSequence_shapes, fitSequence_endIndicies, fitSequence_parameters, fitSequence_bestFitType] = ...
% fcn_geometry_fitSequentialArcsC2(points_to_fit, (fitting_tolerance), (flag_fit_backwards), (fig_num))
%
% INPUTS:
%      points_to_fit: an [Nx2] matrix of N different [x y] points assumed to
%      be in sequence. Note: the function may break, particularly in the
%      calculation of the arcLength, if the points are not in sequence
%
%      (OPTIONAL INPUTS)
% 
%      fitting_tolerance: the distance allowable from a point to an arc fit
%      wherein the point is considered "fitted" to the arc. Default is 0.1
%      meters. If this is entered as a 2x1 or 1x2, then this specifies the
%      threshold in St form, first in the station direction, and then in
%      the transverse direction. For example, an entry of [3 0.2] would
%      have 0.2 meters threshold in the transverse direction, but 3 meters
%      threshold in the station direction.
%
%      flag_fit_backwards: a flag that, if set to 1, causes the fitting
%      process to proceed "backwards", e.g. from the end of the data set to
%      the beginning. Default is 0.
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
%      NOTE: this code is often used with animations and supports a special
%      fig_num input. Namely, the user can pass an array that includes the
%      figure number and a set of figure handles that, if set, will update
%      an animation of the fitting progress into 4 subplots, so that the
%      user can view results. The array is as follows:
%
%           fig_num(1): the figure number of the plot to use
%           fig_num(2): a handle to the h_plotPoints, that animates the
%           points being tested
%           fig_num(3): h_plotPercentage, to animate the percentage
%           agreement
%           fig_num(4): h_plotFitShape to animate the domain fit shape
%
% OUTPUTS:
%
%      fitSequence_points: a cell array saving the XY points that were found in
%      each domain
%
%      fitSequence_shapes: a polyshape object that defines the boundary of the
%      domain of fit for each domain
%
%      fitSequence_endIndicies: the indicies that indicate the end of each
%      domain, of length N+1 where N is the number of domains (the +1 is
%      because the indicies include the start and end indicies)
%
%      fitSequence_parameters: the parameters for each of the arc fits
%
%      fitSequence_bestFitType: the label of the best fit type (all are arcs)
%
%
% DEPENDENCIES:
%
%      fcn_geometry_fitArcRegressionFromHoughFit
%      fcn_geometry_fitArcRegressionFromHoughFitC2
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_fitSequentialArcsC2
% for a full test suite.
%
% This function was written on 2024_04_03 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2024_04_03 - S. Brennan
% -- wrote the code
% 2024_04_17 - S. Brennan
% -- fixed the animation subfigure issues to be consistent with fig_num
% 2024_04_19 - S. Brennan
% -- added ability to provide tolerances in station and transverse
% directions separately
% 2024_06_21 - Sean Brennan
% -- changed segment parameter format to new standard:
%             [
%              base_point_x, 
%              base_point_y, 
%              heading,
%              s_Length,
%             ]
% 2024_07_01 - Sean Brennan
% -- changed fitting_tolerance to St form, with station first then
% transverse (it was the reverse, previously - confusing)

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
        narginchk(1,4);

    end
end


% Does user want to specify fitting_tolerance?
fitting_tolerance = 0.1;
if (2<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        fitting_tolerance = temp;
    end
end

% Does user want to specify flag_fit_backwards?
flag_fit_backwards = 0;
if (3<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        flag_fit_backwards = temp;
    end
end



% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if  (0==flag_max_speed) && (4<= nargin)
    temp = varargin{end};
    if ~isempty(temp)        
        fig_num = temp;
        flag_do_plots = 1;

        % Does user want to specify animation_figure_handles?
        flag_plot_subfigs = 0;

        if length(temp)>1
            fig_num           = temp(1);
            h_plotPoints      = temp(2);
            h_plotPercentage  = temp(3);
            h_plotFitShape    = temp(4);
            flag_plot_subfigs = 1;
        end
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
NtestPoints = length(points_to_fit(:,1));
minNpoints = 10; % The minimum number of points allowed in a fit

% Get transverse tolerance
if length(fitting_tolerance(1,:))>=2
    transverse_tolerance = fitting_tolerance(1,2);
else
    transverse_tolerance = fitting_tolerance(1,1);
end

% Set the fit direction
if flag_fit_backwards
    current_segment_start_index = NtestPoints;
    direction_of_fit = -1;
    current_point_index = NtestPoints - minNpoints;
    absolute_start_index = NtestPoints;
    absolute_end_index   = 1;
    percentage_point_color = [1 0 0];
    regression_point_color = [1 0 0];

else
    current_segment_start_index = 1;
    direction_of_fit = 1;
    current_point_index = minNpoints;
    absolute_start_index = 1;
    absolute_end_index   = NtestPoints;
    percentage_point_color = [0 1 0];
    regression_point_color = [0 1 0];
end



% Initialize arrays and structures
percentage_of_fits = nan(NtestPoints,1);
Hough_domain.best_fit_type    = 'Hough arc';
Hough_domain.best_fit_parameters  = [nan nan nan nan nan 0]; % The zero indicates this is an arc
empty_data = nan*points_to_fit;
Ndomains = 1;
fitSequence_endIndicies{Ndomains} = current_segment_start_index;

% Initialize subfigure plots?
if flag_do_plots==1 && 1==flag_plot_subfigs
    figure(fig_num);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the input points
    subplot(2,2,1);
    hold on;
    grid on;
    axis equal;
    xlabel('X [meters]');
    ylabel('Y [meters]');
    title('Input points');

    % Plot the input points
    current_color = fcn_geometry_fillColorFromNumberOrName(1,'points',[],-1);
    plot(points_to_fit(:,1),points_to_fit(:,2),'.','Color',current_color,'MarkerSize',2);


    % Grab the axis
    original_axis = axis + [-10 10 -10 10];
    axis(original_axis);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the percentage fit per index
    subplot(2,2,2);

    title('Regression fit');
    hold on;
    grid on;
    axis equal
    xlabel('X [meters]');
    ylabel('Y [meters]');

    % Plot the original data
    current_color = fcn_geometry_fillColorFromNumberOrName(1,'points',[],-1);
    plot(points_to_fit(:,1),points_to_fit(:,2),'.','Color',current_color,'MarkerSize',2);

    % Plot the fit shape
    Hough_domain.points_in_domain = points_to_fit(:,1:2);
    Hough_domain.best_fit_source_indicies = [1 2 NtestPoints];
    regression_fit  =  fcn_geometry_fitArcRegressionFromHoughFit(Hough_domain, transverse_tolerance, -1);

    fitShape = regression_fit.best_fit_domain_box;
    current_color = fcn_geometry_fillColorFromNumberOrName(1,[],[],-1);
    h_plotFitShape = plot(fitShape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',1,'EdgeAlpha',0);


    % Plot the "hit" points within the fit
    empty_data = nan*points_to_fit;
    h_plotPoints = plot(empty_data(:,1),empty_data(:,2),'.','Color',[0 1 0],'MarkerSize',10);
    axis(original_axis);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make a plot of percentage of fits
    subplot(2,2,3);
    hold on;
    percentage_of_fits = nan(NtestPoints,1);
    xlabel('Number of points');
    ylabel('Percentage inside');
    title('Percentage of points fit for fits of N points')

    % Plot a bar going across at 100%
    plot((1:NtestPoints)',ones(NtestPoints,1),'k-','LineWidth',5);

    % Create placeholder points to show progress
    h_plotPercentage = plot((1:NtestPoints)',percentage_of_fits,'.');
    axis([0 NtestPoints -0.1 1.1]);
    grid on;




end


% Perform the fitting
% On the first time through, we do not know orientation, so use a
% general Hough fitting method

flag_keep_going = 1;
plotting_increment_interval = 10;

while 1==flag_keep_going
    current_point_index = current_point_index + 1*direction_of_fit;
    fprintf(1,'Checking current point: %.0d of %.0d\n',current_point_index,NtestPoints)

    flag_update_plots = (0==mod(current_point_index,plotting_increment_interval));
    
    % Grab the points in current domain
    current_points_in_domain = points_to_fit(current_segment_start_index:direction_of_fit:current_point_index,1:2);
    test_points_for_domain = points_to_fit(current_segment_start_index:direction_of_fit:absolute_end_index,:);

    % Perform the regression fit of the arc
    Hough_domain.points_in_domain = current_points_in_domain;
    Hough_domain.best_fit_source_indicies = [1 2 length(current_points_in_domain(:,1))];

    % Do regression
    if 1==Ndomains
        regression_domain  =  ...
            fcn_geometry_fitArcRegressionFromHoughFit(Hough_domain, fitting_tolerance, -1);
    else
        regression_domain  =  ...
            fcn_geometry_fitArcRegressionFromHoughFitC2(Hough_domain, previous_parameters, fitting_tolerance, -1);
    end

    % Use isinterior to check which points belong to the fit
    buffered_box = regression_domain.best_fit_domain_box; % polybuffer(regression_domain.best_fit_domain_box,transverse_tolerance);
    indicies_inside_fit = isinterior(buffered_box, test_points_for_domain(:,1),test_points_for_domain(:,2));
    points_in_fit = empty_data;
    points_in_fit(indicies_inside_fit,:) = test_points_for_domain(indicies_inside_fit,:);
    NpointsInCurrentFit = length(current_points_in_domain(:,1));
    percentage_of_fits(current_point_index,1) = min(sum(indicies_inside_fit)/NpointsInCurrentFit,1);

    % Perform plot updates?
    if flag_plot_subfigs && flag_update_plots
        % Update the fit region
        set(h_plotFitShape,'Shape',regression_domain.best_fit_domain_box);
        
        % Update the points inside the fit
        set(h_plotPoints,'XData',points_in_fit(:,1),'YData',points_in_fit(:,2),'Color',regression_point_color);

        % Update the percentage plot
        set(h_plotPercentage,'YData',percentage_of_fits+0.01*direction_of_fit,'Color',percentage_point_color);
    end

    if percentage_of_fits(current_point_index,1)<1
        
        % Save results from previous segment
        fitSequence_endIndicies{Ndomains+1} = current_point_index; 
        fitSequence_points{Ndomains} = current_points_in_domain; %#ok<AGROW>
        fitSequence_shapes{Ndomains} = regression_domain.best_fit_domain_box; %#ok<AGROW>
        fitSequence_parameters{Ndomains} = regression_domain.best_fit_parameters; %#ok<AGROW>
        fitSequence_bestFitType{Ndomains} = regression_domain.best_fit_type; %#ok<AGROW>
        previous_parameters = regression_domain.best_fit_parameters;

        % Set up for next loop
        Ndomains = Ndomains + 1; % Increment the number of domains

        % Move the index "backwards" to restart the fitting
        desired_new_start_point = current_point_index-(minNpoints-1)*direction_of_fit;
        current_segment_start_index = min(max(1,desired_new_start_point),NtestPoints);

        % Perform plot updates?
        if flag_plot_subfigs 
            % Set up plots for next round
            figure(get(h_plotFitShape.Parent.Parent, 'Number'));
            subplot(2,2,2)
            domainShape = regression_domain.best_fit_domain_box;
            current_color = fcn_geometry_fillColorFromNumberOrName(Ndomains,regression_domain.best_fit_type,[],-1);
            h_plotFitShape = plot(domainShape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',1,'EdgeAlpha',0);

            % Add vertical lines to indicate where the segments are changing
            subplot(2,2,3);
            plot([current_point_index current_point_index],[-0.1 1.1],'-','Color',current_color);
        end
    end

    % Check if we are done
    if flag_fit_backwards
        flag_keep_going =  current_point_index > 1;
    else
        flag_keep_going =  current_point_index < length(points_to_fit(:,1));
    end

    pause(0.01)
end

% Save last results
fitSequence_endIndicies{Ndomains+1} = current_point_index;
fitSequence_points{Ndomains} = current_points_in_domain;
fitSequence_shapes{Ndomains} = regression_domain.best_fit_domain_box;
fitSequence_parameters{Ndomains} = regression_domain.best_fit_parameters; 
fitSequence_bestFitType{Ndomains} = regression_domain.best_fit_type; 

if flag_plot_subfigs
    % Add vertical lines to first indicate where the segments are changing
    figure(get(h_plotFitShape.Parent.Parent, 'Number'));
    subplot(2,2,3);
    current_color = fcn_geometry_fillColorFromNumberOrName(Ndomains,fitSequence_bestFitType{Ndomains},[],-1);
    plot([absolute_start_index absolute_start_index],[-0.1 1.1],'-','Color',current_color);
end

% % Rearrange outputs if fitting in backwards order so that the outputs
% % correspond to the first points in the first indicies and last points in
% % the last indicies. As well, flip the ordering of the parameters
% % front/back.
% if flag_fit_backwards
%     % The indicies have Ndomains+1 in length
%     temp_fitSequence_endIndicies = fitSequence_endIndicies;
%     Nindicies = Ndomains+1;
%     for ith_index = 1:Nindicies
%         fitSequence_endIndicies{ith_index} = temp_fitSequence_endIndicies{Nindicies+1-ith_index};
%     end
% 
% 
%     % These inputs all have Ndomains in length
%     temp_fitSequence_points      = fitSequence_points;
%     temp_fitSequence_shapes      = fitSequence_shapes;
%     temp_fitSequence_parameters  = fitSequence_parameters;
%     temp_fitSequence_bestFitType = fitSequence_bestFitType;
%     for ith_domain = 1:Ndomains
%         fitSequence_points{ith_domain}      = temp_fitSequence_points{Ndomains+1-ith_domain};
%         fitSequence_shapes{ith_domain}      = temp_fitSequence_shapes{Ndomains+1-ith_domain};
%         flipped_order_parameters = fcn_geometry_flipGeom('arc',  temp_fitSequence_parameters{Ndomains+1-ith_domain});
%         fitSequence_parameters{ith_domain}  = flipped_order_parameters;
%         fitSequence_bestFitType{ith_domain} = temp_fitSequence_bestFitType{Ndomains+1-ith_domain};
%     end
% 
% end

%% Check for arcs that are really lines
%            'Regression arc' - 
%
%               [circleCenter_x.
%                circleCenter_y,
%                radius,
%                start_angle_in_radians, 
%                end_angle_in_radians,
%                flag_this_is_a_circle
%               ] 

if 1==0
    for ith_domain = 1:length(fitSequence_parameters)
        if strcmp(fitSequence_bestFitType{ith_domain},'Regression arc')
            % Find the arc's height. See diagram here, for example:
            % https://mathcentral.uregina.ca/QQ/database/QQ.09.07/s/bruce1.html
            angle_sweep_radians = diff(fitSequence_parameters{ith_domain}(4:5));
            half_angle = abs(angle_sweep_radians)/2;
            fit_radius = fitSequence_parameters{ith_domain}(3);
            arc_height = fit_radius*(1-cos(half_angle));

            % It's probably a line if the arc almost fits within the box
            % created by the fitting tolerance. The height of that box is
            % 2*tolerance. However, the first and last points may be very
            % slightly outside this, and so we make the height 3 times the
            % tolerance.
            if arc_height < 3*transverse_tolerance(1)
                % This is a line - redo the fit with a line
                Hough_domain.points_in_domain = fitSequence_points{ith_domain};
                Hough_domain.best_fit_source_indicies = [1 length(fitSequence_points{ith_domain}(:,1))];
                Hough_domain.best_fit_type = 'Hough segment';
                first_point = fitSequence_points{ith_domain}(1,:);
                last_point = fitSequence_points{ith_domain}(end,:);
                vector_from_first_to_last = last_point - first_point;
                angle_of_segment = atan2(vector_from_first_to_last(2),vector_from_first_to_last(1));
                length_of_domain = sum(vector_from_first_to_last.^2,2).^0.5;
                Hough_domain.best_fit_parameters = [first_point angle_of_segment length_of_domain];
                regression_domain = fcn_geometry_fitLinearRegressionFromHoughFit(Hough_domain, transverse_tolerance, -1);

                % Update the data with the regression fit
                fitSequence_shapes{ith_domain} = regression_domain.best_fit_domain_box;
                fitSequence_parameters{ith_domain} = regression_domain.best_fit_parameters;
                fitSequence_bestFitType{ith_domain} = regression_domain.best_fit_type;

            end

        end
    end

end

%% Fix ordering
% Rearrange outputs if fitting in backwards order so that the outputs
% correspond to the first points in the first indicies and last points in
% the last indicies. As well, flip the ordering of the parameters
% front/back.
if flag_fit_backwards
    % The indicies have Ndomains+1 in length
    temp_fitSequence_endIndicies = fitSequence_endIndicies;
    Nindicies = Ndomains+1;
    for ith_index = 1:Nindicies
        fitSequence_endIndicies{ith_index} = temp_fitSequence_endIndicies{Nindicies+1-ith_index};
    end


    % These inputs all have Ndomains in length
    temp_fitSequence_points      = fitSequence_points;
    temp_fitSequence_shapes      = fitSequence_shapes;
    temp_fitSequence_parameters  = fitSequence_parameters;
    temp_fitSequence_bestFitType = fitSequence_bestFitType;
    for ith_domain = 1:Ndomains
        fitSequence_points{ith_domain}      = temp_fitSequence_points{Ndomains+1-ith_domain};
        fitSequence_shapes{ith_domain}      = temp_fitSequence_shapes{Ndomains+1-ith_domain};
        if contains(temp_fitSequence_bestFitType{Ndomains+1-ith_domain},'arc')
            flipped_order_parameters = fcn_geometry_flipGeom('arc',  temp_fitSequence_parameters{Ndomains+1-ith_domain});
        else
            flipped_order_parameters = fcn_geometry_flipGeom('segment',  temp_fitSequence_parameters{Ndomains+1-ith_domain});
        end
        fitSequence_parameters{ith_domain}  = flipped_order_parameters;
        fitSequence_bestFitType{ith_domain} = temp_fitSequence_bestFitType{Ndomains+1-ith_domain};
    end

end

%% Clean the parameters
% Make sure the parameters are good

for ith_domain = 1:length(fitSequence_parameters)
    fitSequence_parameters{ith_domain} = fcn_geometry_cleanGeom(fitSequence_bestFitType{ith_domain},  fitSequence_parameters{ith_domain}, -1); 
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

    if flag_plot_subfigs

        figure(get(h_plotFitShape.Parent.Parent, 'Number'));
        
        % Plot the results in the subplot
        flag_rescale_axis = 0;

        % Match subplot 4 axis with that from subplot 1
        subplot(2,2,1);
        original_axis = axis;

    else
        % Plot the results in the given figure number
        temp_h = figure(fig_num);
        flag_rescale_axis = 0;
        if isempty(get(temp_h,'Children'))
            flag_rescale_axis = 1;
        end
    end

    if flag_plot_subfigs == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Final fit
        figure(fig_num);
        subplot(2,2,4);
        hold on;
        grid on;
        axis equal;
        xlabel('X [meters]');
        ylabel('Y [meters]');

        % Plot the fit results 
        for ith_domain = 1:length(fitSequence_points)
            % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain,fitSequence_bestFitType{ith_domain},[],-1);
            current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain,[],[],-1);
            current_fitSequence_points = fitSequence_points{ith_domain};
            current_fitSequence_shape  = fitSequence_shapes{ith_domain};
            plot(current_fitSequence_points(:,1),current_fitSequence_points(:,2),'.','Color',current_color*0.8,'MarkerSize',10);
            plot(current_fitSequence_shape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',1,'EdgeAlpha',0);
        end

        % Plot the domain fits
        fcn_geometry_plotFitSequences(fitSequence_bestFitType, fitSequence_parameters,(fig_num));

        axis(original_axis);


    else
        figure(fig_num);
        hold on;
        grid on;
        axis equal;
        xlabel('X [meters]');
        ylabel('Y [meters]');

        % Plot the input points very large 
        current_color = fcn_geometry_fillColorFromNumberOrName(1,'points',[],-1);
        plot(points_to_fit(:,1),points_to_fit(:,2),'.','Color',current_color,'MarkerSize',10);



        % Plot the domain points
        for ith_domain = 1:length(fitSequence_points)
            % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain,fitSequence_bestFitType{ith_domain},[],-1);
            current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain,[],[],-1);
            current_fitSequence_points = fitSequence_points{ith_domain};
            current_fitSequence_shape  = fitSequence_shapes{ith_domain};
            plot(current_fitSequence_points(:,1),current_fitSequence_points(:,2),'.','Color',current_color*0.8,'MarkerSize',15);
            plot(current_fitSequence_shape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',3,'EdgeAlpha',0);
        end

        % Plot the domain fits
        fcn_geometry_plotFitSequences(fitSequence_bestFitType, fitSequence_parameters,(fig_num));

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