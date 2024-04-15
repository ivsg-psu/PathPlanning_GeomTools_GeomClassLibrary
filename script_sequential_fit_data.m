%% script_sequential_fit_data
% Exercises the function: fcn_geometry_fitArcRegressionFromHoughFit

% 2024_04_02 - S. Brennan
% -- wrote the code
% 2024_04_14 - S. Brennan
% -- separated out functions - URHERE

close all;

%% Use fillArcSequence to create some test data
fig_num = 1;
figure(fig_num);
clf;

rng(1); % Fix the random number, for debugging

% arc_pattern has [1/R and L] for each segment as a row
arc_pattern = [...
    1/20, 15; 
    0 20;
    -1/5 10; 
    0 10;
    1/15 40; 
    0 15
    -1/10 20];

% arc_pattern = [...
%     1/20, 15; 
%     0 20];

M = 10;
sigma = 0.02;

[test_points, ~, ~, trueArcStartIndicies, trueNamedCurveTypes, trueParameters] = fcn_geometry_fillArcSequenceTestPoints(arc_pattern, M, sigma, -1);

% Add noise?
if 1==0
    % Corrupt the results
    probability_of_corruption = 1;
    magnitude_of_corruption = 0.03;

    test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
        (probability_of_corruption), (magnitude_of_corruption), (-1));
end

% Add outliers?
if 1==0
    % Corrupt the results
    probability_of_corruption = 0.1;
    magnitude_of_corruption = 1;

    test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
        (probability_of_corruption), (magnitude_of_corruption), (-1));
end


% Initialize the subplots
subplot_fig_num = fig_num*100;
animation_figure_handles = fcn_INTERNAL_setupSubplots(test_points, trueArcStartIndicies, trueNamedCurveTypes, subplot_fig_num);

% Perform the fit forwards
fitting_tolerance = 0.1; % Units are meters
flag_fit_backwards = 0;
[fitSequence_points_forward, fitSequence_shapes_forward, fitSequence_endIndicies_forward, fitSequence_parameters_forward, fitSequence_bestFitType_forward] = fcn_geometry_fitSequentialArcs(test_points, fitting_tolerance, flag_fit_backwards, animation_figure_handles, fig_num);

% Perform the fit backwards
fitting_tolerance = 0.1; % Units are meters
flag_fit_backwards = 1;
[fitSequence_points_backward, fitSequence_shapes_backward, fitSequence_endIndicies_backward, fitSequence_parameters_backward, fitSequence_bestFitType_backward] = fcn_geometry_fitSequentialArcs(test_points, fitting_tolerance, flag_fit_backwards, animation_figure_handles, fig_num);

% Compare lengths and parameters

NfitsInSequence = length(fitSequence_points_forward);

% First, make absolutely sure that the number of fits found in the forward
% direction match the same number of fits in the backward direction
if length(fitSequence_points_backward)~=NfitsInSequence
    warning('on','backtrace');
    warning('An error will be thrown at this code location as the fits were directionally different.');
    error('Found different numbers of fits in the fit sequence when comparing forward/backward directions. Code is not able to handle this yet!');
end

for ith_fit = 1:NfitsInSequence
    if ~strcmp(fitSequence_bestFitType_forward{ith_fit},fitSequence_bestFitType_backward{ith_fit})
        warning('on','backtrace');
        warning('An error will be thrown at this code location as the fits were found to be geometrically different.');
        error('Found different geometries as the best fit for the same regions when comparing forward/backward directions. Code is not able to handle this yet!');
    end
end


% Find the probable fit
fitSequence_indicies_matrix_forward = cell2mat(fitSequence_endIndicies_forward)';
fitSequence_indicies_matrix_backward = cell2mat(fitSequence_endIndicies_backward)';
probable_arc_boundary_indicies = round(mean([fitSequence_indicies_matrix_forward fitSequence_indicies_matrix_backward],2));
% probable_arc_boundary_indicies = probable_arc_boundary_indicies(1:end-1,:);

% Print the results

fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('Fit number:'),20));
for ith_fit = 1:NfitsInSequence
    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%.0d',ith_fit),10));
end
fprintf(1,'\n');

fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('True start index:'),20));
for ith_fit = 1:NfitsInSequence
    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%.0d',trueArcStartIndicies(ith_fit)),10));
end
fprintf(1,'\n');

fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('Fit start index:'),20));
for ith_fit = 1:NfitsInSequence
    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%.0d',probable_arc_boundary_indicies(ith_fit)),10));
end
fprintf(1,'\n');

% Add vertical lines to indicate where the segments were identified as
% changing
figure(subplot_fig_num);
subplot(2,2,3);
for ith_start = 1:NfitsInSequence
    
    current_color = fcn_geometry_fillColorFromNumberOrName(ith_start,fitSequence_bestFitType_forward{ith_start},-1);

    plot([probable_arc_boundary_indicies(ith_start) probable_arc_boundary_indicies(ith_start)],[-0.1 1.1],'-','Color',current_color);
end

figure(subplot_fig_num);
subplot(2,2,4);

% Plot the fitted groups of points. If any of the points are mis-labeled,
% there will be one color incorrectly on top of another, for example a red
% point on top of a blue underlying point.
for ith_plot = 1:NfitsInSequence
    current_color = fcn_geometry_fillColorFromNumberOrName(ith_plot,fitSequence_bestFitType_forward{ith_plot},-1);
    index_range = probable_arc_boundary_indicies(ith_plot):probable_arc_boundary_indicies(ith_plot+1);
    plot(test_points(index_range,1),test_points(index_range,2),'.','Color',current_color,'MarkerSize',10);
end

%% Connect the fits so that the lines perfectly align with the arcs

fig_num = 23456;

revised_fitSequence_parameters_forward  = fcn_geometry_alignGeometriesInSequence(fitSequence_bestFitType_forward, fitSequence_parameters_forward, fitting_tolerance*2, fig_num);
revised_fitSequence_parameters_backward = fcn_geometry_alignGeometriesInSequence(fitSequence_bestFitType_backward,fitSequence_parameters_backward, fitting_tolerance*2, fig_num);

fcn_geometry_plotFitSequences(fitSequence_bestFitType_forward, revised_fitSequence_parameters_forward,(fig_num));
fcn_geometry_plotFitSequences(fitSequence_bestFitType_backward, revised_fitSequence_parameters_backward,(fig_num));


%% Print the results


fprintf(1,'\n\nPARAMETER FIT COMPARISON:\n');
for ith_fit = 1:NfitsInSequence
    fprintf(1,'\n\nFit Sequence Number: %.0d\n', ith_fit); 

    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('   '),10));
    fcn_INTERNAL_printFitDetails(trueNamedCurveTypes{ith_fit},trueParameters{ith_fit},1)

    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('TRUE '),10));
    fcn_INTERNAL_printFitDetails(trueNamedCurveTypes{ith_fit},trueParameters{ith_fit},0)

    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('FORWARD '),10));
    fcn_INTERNAL_printFitDetails(fitSequence_bestFitType_forward{ith_fit},revised_fitSequence_parameters_forward{ith_fit},0)

    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('REVERSE '),10));
    fcn_INTERNAL_printFitDetails(fitSequence_bestFitType_backward{ith_fit},revised_fitSequence_parameters_backward{ith_fit},0)
end
fprintf(1,'\n');



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

%% fcn_INTERNAL_printFitDetails
function fcn_INTERNAL_printFitDetails(fit_type,fit_parameters, flag_print_header)

if contains(fit_type,{'arc','circle'})
    print_type = 'arc';            
    header_strings = {'centerX','centerY','radius','startAngle','endAngle','isCircle','turnsLeft'};
elseif contains(fit_type,{'segment','line'})
    print_type = 'line';
    header_strings = {'vect_x','vect_y','start_x','start_y','start_station','end_station'};
end

NumColumnChars = 15;

% Print header?
if flag_print_header
    % Print values
    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf(' '),20));

    parameter_string = '';
    for ith_parameter = 1:length(fit_parameters)
        number_string = fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%s ',header_strings{ith_parameter}),NumColumnChars+1);
        parameter_string = cat(2,parameter_string,number_string);
    end
    fprintf(1,'%s\n',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%s',parameter_string),7*NumColumnChars));
else
    % Print values
    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('Fit type: %s ',print_type),20));

    % Print parameters
    parameter_string = '';
    for ith_parameter = 1:length(fit_parameters)
        number_string = fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%.4f ',fit_parameters(ith_parameter)),NumColumnChars);
        if strcmp(number_string(1),'-')
            parameter_string = cat(2,parameter_string,number_string, ' ');
        else
            parameter_string = cat(2,parameter_string,' ', number_string);
        end
    end
    fprintf(1,'%s\n',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%s',parameter_string),7*NumColumnChars));
end

end % Ends fcn_INTERNAL_printFitDetails

%% fcn_geometry_alignGeometriesInSequence
function revised_fitSequence_parameters = fcn_geometry_alignGeometriesInSequence(fitSequence_bestFitType, fitSequence_parameters, threshold, fig_num)

NfitsInSequence = length(fitSequence_bestFitType);

revised_fitSequence_parameters = fitSequence_parameters;

% Loop through fits, connecting them together
for ith_fit = 1:NfitsInSequence-1
    current_fit_type = fitSequence_bestFitType{ith_fit};
    next_fit_type    = fitSequence_bestFitType{ith_fit+1};

    if strcmp(current_fit_type,'Regression arc') && strcmp(next_fit_type,'Regression arc')
        error('Joining arcs to arcs is not yet implemented');
    elseif strcmp(current_fit_type,'Vector regression segment fit') && strcmp(next_fit_type,'Vector regression segment fit')
        error('Joining line segements to line segments is not yet implemented');
    elseif (strcmp(current_fit_type,'Regression arc') && strcmp(next_fit_type,'Vector regression segment fit')) || (strcmp(current_fit_type,'Vector regression segment fit') && strcmp(next_fit_type,'Regression arc'))
        % This is a connection that is an arc to line or line to arc

        % find which index is which
        if strcmp(current_fit_type,'Regression arc')
            line_index = ith_fit+1;
            arc_index  = ith_fit;
            flag_arc_is_first = 1;
        else
            line_index = ith_fit;
            arc_index  = ith_fit+1;
            flag_arc_is_first = 0;
        end

        % Fix connections of lines to arcs and arcs to lines
        [revised_fitSequence_parameters{line_index}, revised_fitSequence_parameters{arc_index}] = fcn_geometry_joinLineToArc(fitSequence_parameters{line_index}, fitSequence_parameters{arc_index}, flag_arc_is_first, (threshold),(fig_num));

        % Need to update the 2nd one for the next iteration
        if flag_arc_is_first
            fitSequence_parameters{line_index} = revised_fitSequence_parameters{line_index};
        else
            fitSequence_parameters{arc_index} = revised_fitSequence_parameters{arc_index};
        end

    else
        error('unrecognized pattern encountered! Exiting.')
    end
end % Ends looping through fits
end % Ends fcn_geometry_alignGeometriesInSequence



%% fcn_INTERNAL_setupSubplots
function figure_handles = fcn_INTERNAL_setupSubplots(test_points, arcStartIndicies, namedCurveTypes, subplot_fig_num)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE TEST POINTS
figure(subplot_fig_num);
clf;
subplot(2,2,1);
hold on;
grid on;
axis equal;
xlabel('X [meters]');
ylabel('Y [meters]');

% Plot the groups of points
modifiedArcStartIndicies = [arcStartIndicies; length(test_points(:,1))];
for ith_plot = 1:length(arcStartIndicies(:,1))
    current_color = fcn_geometry_fillColorFromNumberOrName(ith_plot,namedCurveTypes{ith_plot},-1);
    index_range = modifiedArcStartIndicies(ith_plot):modifiedArcStartIndicies(ith_plot+1);
    plot(test_points(index_range,1),test_points(index_range,2),'.','Color',current_color,'MarkerSize',10);
end

% Grab the axis
original_axis = axis + [-10 10 -10 10];
axis(original_axis);

% Label the plot
figure(subplot_fig_num);
subplot(2,2,1);
title('Input points');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot that shows sequential segment fitting
figure(subplot_fig_num);
subplot(2,2,2);


NtestPoints = length(test_points(:,1));
Hough_fit.best_fit_type    = 'Hough arc';
Hough_fit.best_fit_parameters  = [nan nan nan nan nan 0]; % The zero indicates this is an arc


figure(subplot_fig_num);
subplot(2,2,2);

title('Regression fit');
hold on;
grid on;
axis equal
xlabel('X [meters]');
ylabel('Y [meters]');

% Plot the original data
plot(test_points(:,1),test_points(:,2),'.','Color',[0 0 0],'MarkerSize',5);

% Plot the fit shape
Hough_fit.points_in_domain = test_points(:,1:2);
Hough_fit.best_fit_source_indicies = [1 2 NtestPoints];
regression_fit  =  ...
    fcn_geometry_fitArcRegressionFromHoughFit(Hough_fit, 0.1, -1);
fitShape = regression_fit.best_fit_domain_box;
current_color = fcn_geometry_fillColorFromNumberOrName(1,[],-1);
h_plotFitShape = plot(fitShape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',1,'EdgeAlpha',0);

% Plot the "hit" points within the fit
empty_data = nan*test_points;
h_plotPoints = plot(empty_data(:,1),empty_data(:,2),'.','Color',[0 1 0],'MarkerSize',10);
axis(original_axis);

% Make a plot of percentage of fits
figure(subplot_fig_num);
subplot(2,2,3);
hold on;
percentage_of_fits = nan(NtestPoints,1);
xlabel('Number of points');
ylabel('Percentage inside');

% Plot a bar going across at 100%
plot((1:NtestPoints)',ones(NtestPoints,1),'k-','LineWidth',5);

% Create placeholder points to show progress
h_plotPercentage = plot((1:NtestPoints)',percentage_of_fits,'.');
axis([0 NtestPoints -0.1 1.1]);
grid on;

% Add vertical lines to indicate where the segments are TRUELY changing
figure(subplot_fig_num);
subplot(2,2,3);
for ith_start = 1:length(arcStartIndicies)
    plot([arcStartIndicies(ith_start) arcStartIndicies(ith_start)],[-0.1 1.1],'k-','LineWidth',5);
end

figure_handles(1) = h_plotPoints;
figure_handles(2) = h_plotPercentage;
figure_handles(3) = h_plotFitShape;

figure(subplot_fig_num);
subplot(2,2,4);
hold on;
grid on;
axis equal;
xlabel('X [meters]');
ylabel('Y [meters]');

% Plot the groups of points
modifiedArcStartIndicies = [arcStartIndicies; length(test_points(:,1))];
for ith_plot = 1:length(arcStartIndicies(:,1))
    current_color = fcn_geometry_fillColorFromNumberOrName(ith_plot,namedCurveTypes{ith_plot},-1);
    index_range = modifiedArcStartIndicies(ith_plot):modifiedArcStartIndicies(ith_plot+1);
    plot(test_points(index_range,1),test_points(index_range,2),'.','Color',current_color,'MarkerSize',30);
end


end % ends fcn_INTERNAL_setupSubplots

function [fitSequence_points, fitSequence_shapes, fitSequence_endIndicies, fitSequence_parameters, fitSequence_bestFitType] = fcn_geometry_fitSequentialArcs(points_to_fit, varargin)
%% fcn_geometry_fitSequentialArcs
% Given a set of XY data, attempts to fit the data in sequential order with
% an arc until the points in the fit fall outside of a fitting tolerance.
% Note that line segment fitting also occurs as line segments are simply
% arcs with infinite radius and this function can work with infinite radius
% arc segments. 
% 
% The function proceeds point-by-point from one end of data to the other
% attempting a fit of points in sequence until the fitted data no longer
% fall within the fitting tolerance. Once the tolerance is violated, a new
% arc fit is started to create a new domain of fitting. When the end of the
% point sequence is reached, the function returns information about each of
% the domains where an arc fit was found. The results are useful to quickly
% approximate XY data as a sequence of joined arcs and lines.
%
% A flag option allows the search to proceed either from first point to end
% point (default) or from last point to first point. If proceeding from the
% last point as the starting point, the outputs are rearranged so that the
% domains and respective indicies are ordered from first point to last
% point, allowing easy comparison to results obtained from fitting in the
% typical first-point as starting piont.
% 
% Format: 
% [fitSequence_points, fitSequence_shapes, fitSequence_endIndicies, fitSequence_parameters, fitSequence_bestFitType] = fcn_geometry_fitSequentialArcs(points_to_fit, (fitting_tolerance), (flag_fit_backwards), (animation_figure_handles),(fig_num))
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
%      meters.
%
%      flag_fit_backwards: a flag that, if set to 1, causes the fitting
%      process to proceed "backwards", e.g. from the end of the data set to
%      the beginning. Default is 0.
%
%      animation_figure_handles: a set of figure handles that, if set, will
%      update an animation of the fitting progress into 4 subplots, so that
%      the user can view results. The handles are stored in an array:
%
%           h_plotPoints      = animation_figure_handles(1);
%           h_plotPercentage  = animation_figure_handles(2);
%           h_plotFitShape = animation_figure_handles(3);
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
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
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_fitSequentialArcs
% for a full test suite.
%
% This function was written on 2024_04_03 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2024_04_03 - S Brennan
% -- wrote the code

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
        narginchk(1,5);

    end
end


% Does user want to specify flag_fit_backwards?
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

% Does user want to specify animation_figure_handles?
flag_do_animations = 0;
if (4<=nargin)
    temp = varargin{3};
    if ~isempty(temp)
        animation_figure_handles = temp;
        flag_do_animations = 1;
    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if  (0==flag_max_speed) && (5<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
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
Ndomains = 1;

% Prep for animations?
if flag_do_animations
    h_plotPoints      = animation_figure_handles(1);
    h_plotPercentage  = animation_figure_handles(2);
    h_plotFitShape = animation_figure_handles(3);       
end

% Set the fit direction
if flag_fit_backwards
    current_segment_start_index = NtestPoints;
    direction_of_fit = -1;
    current_point_index = NtestPoints - 3;
    absolute_start_index = NtestPoints;
    absolute_end_index   = 1;
    percentage_point_color = [1 0 0];
    regression_point_color = [1 0 0];

else
    current_segment_start_index = 1;
    direction_of_fit = 1;
    current_point_index = 3;
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
fitSequence_endIndicies{Ndomains} = current_segment_start_index;


% Perform the fitting
% On the first time through, we do not know orientation, so use a
% general Hough fitting method

flag_keep_going = 1;
plotting_increment_interval = 10;

while 1==flag_keep_going
    current_point_index = current_point_index + 1*direction_of_fit;

    flag_update_plots = (0==mod(current_point_index,plotting_increment_interval));
    
    % Grab the points in current domain
    current_points_in_domain = points_to_fit(current_segment_start_index:direction_of_fit:current_point_index,1:2);
    test_points_for_domain = points_to_fit(current_segment_start_index:direction_of_fit:absolute_end_index,:);

    % Perform the regression fit of the arc
    Hough_domain.points_in_domain = current_points_in_domain;
    Hough_domain.best_fit_source_indicies = [1 2 length(current_points_in_domain(:,1))];
    regression_domain  =  ...
        fcn_geometry_fitArcRegressionFromHoughFit(Hough_domain, fitting_tolerance, -1);

    % Use isinterior to check which points belong to the fit
    buffered_box = regression_domain.best_fit_domain_box; % polybuffer(regression_domain.best_fit_domain_box,fitting_tolerance);
    indicies_inside_fit = isinterior(buffered_box, test_points_for_domain(:,1),test_points_for_domain(:,2));
    points_in_fit = empty_data;
    points_in_fit(indicies_inside_fit,:) = test_points_for_domain(indicies_inside_fit,:);
    NpointsInCurrentFit = length(current_points_in_domain(:,1));
    percentage_of_fits(current_point_index,1) = min(sum(indicies_inside_fit)/NpointsInCurrentFit,1);

    % Perform plot updates?
    if flag_do_animations && flag_update_plots
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

        % Set up for next loop
        Ndomains = Ndomains + 1;        
        current_segment_start_index = min(max(1,current_point_index-2*direction_of_fit),NtestPoints);

        % Perform plot updates?
        if flag_do_animations 
            % Set up plots for next round
            figure(get(h_plotFitShape.Parent.Parent, 'Number'));
            subplot(2,2,2)
            domainShape = regression_domain.best_fit_domain_box;
            current_color = fcn_geometry_fillColorFromNumberOrName(Ndomains,regression_domain.best_fit_type,-1);
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

if flag_do_animations
    % Add vertical lines to first indicate where the segments are changing
    figure(get(h_plotFitShape.Parent.Parent, 'Number'));
    subplot(2,2,3);
    current_color = fcn_geometry_fillColorFromNumberOrName(Ndomains,fitSequence_bestFitType{Ndomains},-1);
    plot([absolute_start_index absolute_start_index],[-0.1 1.1],'-','Color',current_color);
end

% Rearrange outputs if fitting in backwards order so that the outputs
% correspond to the first points in the first indicies and last points in
% the last indicies.
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
        fitSequence_parameters{ith_domain}  = temp_fitSequence_parameters{Ndomains+1-ith_domain};
        fitSequence_bestFitType{ith_domain} = temp_fitSequence_bestFitType{Ndomains+1-ith_domain};
    end

end

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
        if arc_height < 3*fitting_tolerance
            % This is a line - redo the fit with a line
            Hough_domain.points_in_domain = fitSequence_points{ith_domain};
            Hough_domain.best_fit_source_indicies = [1 length(fitSequence_points{ith_domain}(:,1))];
            Hough_domain.best_fit_type = 'Hough segment';
            regression_domain = fcn_geometry_fitLinearRegressionFromHoughFit(Hough_domain, fitting_tolerance, -1);
            
            % Update the data with the regression fit
            fitSequence_shapes{ith_domain} = regression_domain.best_fit_domain_box; 
            fitSequence_parameters{ith_domain} = regression_domain.best_fit_parameters; 
            fitSequence_bestFitType{ith_domain} = regression_domain.best_fit_type; 

        end
        
    end
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

    if flag_do_animations

        figure(get(h_plotFitShape.Parent.Parent, 'Number'));
        
        % Plot the results in the subplot
        flag_rescale_axis = 0;

        % Match subplot 4 axis with that from subplot 1
        subplot(2,2,1);
        original_axis = axis;
        subplot(2,2,4);
        axis(original_axis);

    else
        % Plot the results in the given figure number
        temp_h = figure(fig_num);
        flag_rescale_axis = 0;
        if isempty(get(temp_h,'Children'))
            flag_rescale_axis = 1;
        end
    end

    hold on;
    grid on;
    axis equal;

    title('Domains that were found');
    xlabel('X [meters]');
    ylabel('Y [meters]');

    % Plot the original data
    plot(points_to_fit(:,1),points_to_fit(:,2),'.','Color',[0 0 0],'MarkerSize',5);

    % Plot the domain points
    for ith_domain = 1:length(fitSequence_points)        
        current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain,fitSequence_bestFitType{ith_domain},-1);
        current_fitSequence_points = fitSequence_points{ith_domain};
        current_fitSequence_shape  = fitSequence_shapes{ith_domain};
        plot(current_fitSequence_points(:,1),current_fitSequence_points(:,2),'.','Color',current_color*0.8,'MarkerSize',10);
        plot(current_fitSequence_shape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',1,'EdgeAlpha',0);
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


