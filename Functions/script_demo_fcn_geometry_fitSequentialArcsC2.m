%% script_test_fcn_geometry_fitSequentialArcsC2
% Exercises the function: fcn_geometry_fitSequentialArcsC2

% 2024_04_14 - S. Brennan
% -- wrote the code
% 2024_07_06 - S. Brennan
% -- typo fix in function calls

close all;

%% Artificial data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                _   _  __ _      _       _    _____        _
%     /\        | | (_)/ _(_)    (_)     | |  |  __ \      | |
%    /  \   _ __| |_ _| |_ _  ___ _  __ _| |  | |  | | __ _| |_ __ _
%   / /\ \ | '__| __| |  _| |/ __| |/ _` | |  | |  | |/ _` | __/ _` |
%  / ____ \| |  | |_| | | | | (__| | (_| | |  | |__| | (_| | || (_| |
% /_/    \_\_|   \__|_|_| |_|\___|_|\__,_|_|  |_____/ \__,_|\__\__,_|
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Artificial%20%20Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the random number generator, for repeatability
rng(1);

%%% Use fillArcSequence to create some test data
fig_num = 1;
figure(fig_num);
clf;

rng(1); % Fix the random number, for debugging

% arc_pattern has [1/R and L] for each segment as a row
% arc_pattern = [...
%     1/20, 15; 
%     0 20;
%     -1/5 10; 
%     0 10;
%     1/15 40; 
%     0 15
%     -1/10 20];

arc_pattern = [...
    1/20, 15; 
    -1/40 20];

M = 10; % How many points per meter
sigma = 0.02; % The standard deviation in the points relative to the perfect function fit, in meters

[test_points, ~, ~, trueArcStartIndicies, trueNamedCurveTypes, trueParameters] = fcn_geometry_fillArcSequenceTestPoints(arc_pattern, M, sigma, -1);

% Add more noise?
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
fig_num_array(1) = fig_num;
fig_num_array(2) = 0;
fig_num_array(3) = 0;
fig_num_array(4) = 0;
figure(fig_num); clf;

% Add starter points (truth) onto subplot 2,1,1
subplot(2,2,1);
hold on;
grid on;
axis equal;
xlabel('X [meters]');
ylabel('Y [meters]');

% Plot the groups of true points
modifiedArcStartIndicies = [trueArcStartIndicies; length(test_points(:,1))];
for ith_plot = 1:length(trueArcStartIndicies(:,1))
    if ~isempty(trueNamedCurveTypes)
        current_color = fcn_geometry_fillColorFromNumberOrName(ith_plot,trueNamedCurveTypes{ith_plot},[],-1);
    else
        current_color = [0 0 0];
    end
    index_range = modifiedArcStartIndicies(ith_plot):modifiedArcStartIndicies(ith_plot+1);
    plot(test_points(index_range,1),test_points(index_range,2),'.','Color',current_color,'MarkerSize',10);
end


%%% Perform the fit forwards
fitting_tolerance = [1 0.2]; % Station is 1 meter, transverse is 0.1 meter 
flag_fit_backwards = 0;
[fitSequence_points_forward, fitSequence_shapes_forward, fitSequence_endIndicies_forward, fitSequence_parameters_forward, fitSequence_bestFitType_forward] = ...
    fcn_geometry_fitSequentialArcsC2(test_points, fitting_tolerance, flag_fit_backwards, fig_num_array);

% Plot the true results
subplot(2,2,4);
fcn_geometry_plotFitSequences(trueNamedCurveTypes, trueParameters,(fig_num_array(1)));


%%% Perform the fit backwards
fitting_tolerance = [1 0.2]; % Station is 1 meter, transverse is 0.1 meter 
flag_fit_backwards = 1;
[fitSequence_points_backward, fitSequence_shapes_backward, fitSequence_endIndicies_backward, fitSequence_parameters_backward, fitSequence_bestFitType_backward] = ...
    fcn_geometry_fitSequentialArcsC2(test_points, fitting_tolerance, flag_fit_backwards, fig_num_array);

% Plot the true results
subplot(2,2,4);
fcn_geometry_plotFitSequences(trueNamedCurveTypes, trueParameters,(fig_num_array(1)));

% Print the results
fcn_geometry_printFitSequences(fitSequence_bestFitType_forward, fitSequence_parameters_forward, (1), ('Forward fit'), (1))
fcn_geometry_printFitSequences(fitSequence_bestFitType_backward, fitSequence_parameters_backward, (1), ('Backward fit'), (1))


%%% Merge the forward and backward fits?

%%% Compare lengths and parameters
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
figure(fig_num_array(1));
subplot(2,2,3);
for ith_start = 1:NfitsInSequence
    
    current_color = fcn_geometry_fillColorFromNumberOrName(ith_start,fitSequence_bestFitType_forward{ith_start},[],-1);

    plot([probable_arc_boundary_indicies(ith_start) probable_arc_boundary_indicies(ith_start)],[-0.1 1.1],'-','Color',current_color);
end

figure(fig_num_array(1));
subplot(2,2,4);

% Plot the fitted groups of points. If any of the points are mis-labeled,
% there will be one color incorrectly on top of another, for example a red
% point on top of a blue underlying point.
for ith_plot = 1:NfitsInSequence
    current_color = fcn_geometry_fillColorFromNumberOrName(ith_plot,fitSequence_bestFitType_forward{ith_plot},[],-1);
    index_range = probable_arc_boundary_indicies(ith_plot):probable_arc_boundary_indicies(ith_plot+1);
    plot(test_points(index_range,1),test_points(index_range,2),'.','Color',current_color,'MarkerSize',10);
end


% %% Connect the fits so that the lines perfectly align with the arcs?
% %
% % Check revised_fitSequence_parameters_backward: The backward parameter inputs are not actually backward
% % paramters. Check the inputs. 
% 
% % fig_num = 23456;
% % figure(fig_num);clf;
% % 
% % fitting_tolerance = 5;
% % continuity_level = 2;
% % revised_fitSequence_parameters_forward  = fcn_geometry_alignGeometriesInSequence(fitSequence_bestFitType_forward, fitSequence_parameters_forward,  fitting_tolerance, continuity_level, fig_num);
% % revised_fitSequence_parameters_backward = fcn_geometry_alignGeometriesInSequence(fitSequence_bestFitType_backward,fitSequence_parameters_backward, fitting_tolerance, continuity_level, fig_num);
% % 
% % fcn_geometry_plotFitSequences(fitSequence_bestFitType_forward, revised_fitSequence_parameters_forward,(fig_num));
% % fcn_geometry_plotFitSequences(fitSequence_bestFitType_backward, revised_fitSequence_parameters_backward,(fig_num));
% % 
% % subplot(1,2,1);
% % good_axis_limits = axis;

% Plot the results
% Which results to plot? Comment out one set or another
parameters_forward = fitSequence_parameters_forward;
parameters_backward = fitSequence_parameters_backward;

% parameters_forward = revised_fitSequence_parameters_forward;
% parameters_backward = revised_fitSequence_parameters_backward;



% Plot the original data
fig_num = 23456;
figure(fig_num);clf;

fcn_geometry_plotFitSequences(fitSequence_bestFitType_forward, parameters_forward,(fig_num));
fcn_geometry_plotFitSequences(fitSequence_bestFitType_backward, parameters_backward,(fig_num));
good_axis_limits = axis;


%%% Find avarage parameters
for ith_fit = 1:NfitsInSequence
    fitSequence_parameters_averaged{ith_fit} = (fitSequence_parameters_forward{ith_fit} + fitSequence_parameters_backward{ith_fit})/2; %#ok<SAGROW>
end
fitSequence_bestFitType_averaged = fitSequence_bestFitType_forward;

%%% Plot the errors between true curves and fitted curves
comparison_fig_num = 2828;
figure(comparison_fig_num); clf;
hold on;

threshold = 0.15;

max_forward_error = -inf;
max_backward_error = -inf;
max_averaged_error = -inf;


for ith_fit = 1:NfitsInSequence
    subplot(1,3,1);
    curve_test_segment_length = [];    
    [flag_is_similar, points_XY_on_test_curve, minimum_distance_to_each_point, mean_error, max_error, std_dev_error] = ...
    fcn_geometry_compareCurves(...
    trueNamedCurveTypes{ith_fit}, trueParameters{ith_fit}, fitSequence_bestFitType_forward{ith_fit}, parameters_forward{ith_fit},...
    (threshold), (curve_test_segment_length), (comparison_fig_num));

    max_forward_error = max(max_forward_error,max_error);
    title('Forward fitting');
    axis(good_axis_limits)

    subplot(1,3,2);
    curve_test_segment_length = [];    
    [flag_is_similar, points_XY_on_test_curve, minimum_distance_to_each_point, mean_error, max_error, std_dev_error] = ...
    fcn_geometry_compareCurves(...
    trueNamedCurveTypes{ith_fit}, trueParameters{ith_fit}, fitSequence_bestFitType_backward{ith_fit}, parameters_backward{ith_fit},...
    (threshold), (curve_test_segment_length), (comparison_fig_num));
    max_backward_error = max(max_backward_error,max_error);
    title('Reverse fitting');
    axis(good_axis_limits)

    subplot(1,3,3);
    curve_test_segment_length = [];    
    [flag_is_similar, points_XY_on_test_curve, minimum_distance_to_each_point, mean_error, max_error, std_dev_error] = ...
    fcn_geometry_compareCurves(...
    trueNamedCurveTypes{ith_fit}, trueParameters{ith_fit}, fitSequence_bestFitType_averaged{ith_fit}, fitSequence_parameters_averaged{ith_fit},...
    (threshold), (curve_test_segment_length), (comparison_fig_num));
    max_averaged_error = max(max_averaged_error,max_error);
    title('Averaged fitting');
    axis(good_axis_limits)
end

sgtitle({'Fit quality comparing true curves to fitted curves',sprintf('Red is %.2fm error, blue is 0m error', threshold)});


fprintf(1,'\n');

% Print the results
fcn_geometry_printFitSequences(trueNamedCurveTypes, trueParameters, (1), ('True fit'), (1))
fcn_geometry_printFitSequences(fitSequence_bestFitType_forward, fitSequence_parameters_forward, (1), ('Forward fit'), (1))
fcn_geometry_printFitSequences(fitSequence_bestFitType_backward, fitSequence_parameters_backward, (1), ('Backward fit'), (1))
fcn_geometry_printFitSequences(fitSequence_bestFitType_backward, fitSequence_parameters_averaged, (1), ('Average fit'), (1))


fprintf(1,'Max forward fitting error:  %.3f meters\n',max_forward_error);
fprintf(1,'Max backward fitting error: %.3f meters\n',max_backward_error);
fprintf(1,'Max averaged fitting error: %.3f meters\n',max_averaged_error);

%% Plot the errors between original points and fitted curves
comparison_fig_num = 73434;
figure(comparison_fig_num); clf;
hold on;

threshold = 0.15;
curve_test_segment_length = [];

subplot(1,3,1);

[flag_is_similar_forward, ...
    minimum_distance_to_each_point,...
    indicies_of_nearest_reference_points, ...
    mean_error_forward, max_error_forward, std_dev_error_forward] = ...
    fcn_geometry_comparePointsToCurve(...
    fitSequence_bestFitType_forward, fitSequence_parameters_forward, ...
    test_points, ...
    (threshold), (curve_test_segment_length), (comparison_fig_num));

title('Forward fitting');
axis(good_axis_limits)

subplot(1,3,2);
[flag_is_similar_backward, ...
    minimum_distance_to_each_point,...
    indicies_of_nearest_reference_points, ...
    mean_error_backward, max_error_backward, std_dev_error_backward] = ...
    fcn_geometry_comparePointsToCurve(...
    fitSequence_bestFitType_backward, fitSequence_parameters_backward, ...
    test_points, ...
    (threshold), (curve_test_segment_length), (comparison_fig_num));

title('Backward fitting');
axis(good_axis_limits)

subplot(1,3,3);
[flag_is_similar_averaged, ...
    minimum_distance_to_each_point,...
    indicies_of_nearest_reference_points, ...
    mean_error_averaged, max_error_averaged, std_dev_error_averaged] = ...
    fcn_geometry_comparePointsToCurve(...
    fitSequence_bestFitType_averaged, fitSequence_parameters_averaged, ...
    test_points, ...
    (threshold), (curve_test_segment_length), (comparison_fig_num));

title('Averaged fitting');
axis(good_axis_limits)


sgtitle({'Fit quality comparing test points to fitted curves',sprintf('Red is %.2fm error, blue is 0m error', threshold)});


% Print the results

fprintf(1,'\n');

% Print the results
fcn_geometry_printFitSequences(trueNamedCurveTypes, trueParameters, (1), ('True fit'), (1))
fcn_geometry_printFitSequences(fitSequence_bestFitType_forward, fitSequence_parameters_forward, (1), ('Forward fit'), (1))
fcn_geometry_printFitSequences(fitSequence_bestFitType_backward, fitSequence_parameters_backward, (1), ('Backward fit'), (1))
fcn_geometry_printFitSequences(fitSequence_bestFitType_backward, fitSequence_parameters_averaged, (1), ('Average fit'), (1))


fprintf(1,'IsSimilar forward fitting:  %.3f meters\n',flag_is_similar_forward);
fprintf(1,'IsSimilar backward fitting:  %.3f meters\n',flag_is_similar_backward);
fprintf(1,'IsSimilar averaged fitting:  %.3f meters\n',flag_is_similar_averaged);
fprintf(1,'Max forward fitting error:  %.3f meters\n',max_error_forward);
fprintf(1,'Max backward fitting error: %.3f meters\n',max_error_backward);
fprintf(1,'Max averaged fitting error: %.3f meters\n',max_error_averaged);
fprintf(1,'Mean forward fitting error:  %.3f meters\n',mean_error_forward);
fprintf(1,'Mean backward fitting error: %.3f meters\n',mean_error_backward);
fprintf(1,'Mean averaged fitting error: %.3f meters\n',mean_error_averaged);
fprintf(1,'Std forward fitting error:  %.3f meters\n',std_dev_error_forward);
fprintf(1,'Std backward fitting error: %.3f meters\n',std_dev_error_backward);
fprintf(1,'Std averaged fitting error: %.3f meters\n',std_dev_error_averaged);


%% Test track data
% Initialize the subplots
fig_num = 9999;

fig_num_array(1) = fig_num;
fig_num_array(2) = 0;
fig_num_array(3) = 0;
fig_num_array(4) = 0;
figure(fig_num); clf;


% Check to see if XY data for the centerline of the original track lane was loaded earlier
mat_filename = fullfile(cd,'Data','Centerline_OriginalTrackLane_InnerMarkerClusterCenterOfDoubleYellow.mat');
if exist(mat_filename,'file')
    load(mat_filename,'XY_data');
end

% Since the XY data is very dense, keep only 1 of every "keep_every" points
keep_every = 20;
indicies = (1:length(XY_data(:,1)))';
small_XY_data_indicies = find(0==mod(indicies,keep_every));
small_XY_data = XY_data(small_XY_data_indicies,:);

% Perform the fit forwards
fitting_tolerance = [1 10]; % Units are meters
flag_fit_backwards = 0;
[fitSequence_points_forward, fitSequence_shapes_forward, fitSequence_endIndicies_forward, fitSequence_parameters_forward, fitSequence_bestFitType_forward] = ...
    fcn_geometry_fitSequentialArcsC2(small_XY_data, fitting_tolerance, flag_fit_backwards, fig_num_array);


% Perform the fit backwards
flag_fit_backwards = 1;
[fitSequence_points_backward, fitSequence_shapes_backward, fitSequence_endIndicies_backward, fitSequence_parameters_backward, fitSequence_bestFitType_backward] = ...
    fcn_geometry_fitSequentialArcsC2(small_XY_data, fitting_tolerance, flag_fit_backwards, fig_num_array);

%% Plot results
fig_num = 234343;
figure(fig_num);
clf;

subplot(1,2,1);
hold on;
grid on;
axis equal;
xlabel('X [meters]');
ylabel('Y [meters]');
title('Forward fit');


% Plot the input points very large
current_color = fcn_geometry_fillColorFromNumberOrName(1,'points',[],-1);
plot(small_XY_data(:,1),small_XY_data(:,2),'.','Color',current_color,'MarkerSize',10);

% Plot the domain points
for ith_domain = 1:length(fitSequence_points_forward)
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain,fitSequence_bestFitType{ith_domain},[],-1);
    current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain,[],[],-1);
    current_fitSequence_points = fitSequence_points_forward{ith_domain};
    current_fitSequence_shape  = fitSequence_shapes_forward{ith_domain};
    plot(current_fitSequence_points(:,1),current_fitSequence_points(:,2),'.','Color',current_color*0.8,'MarkerSize',15);
    plot(current_fitSequence_shape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',3,'EdgeAlpha',0);
end

% Plot the domain fits
fcn_geometry_plotFitSequences(fitSequence_bestFitType_forward, fitSequence_parameters_forward,(fig_num));

% Make axis slightly larger?
temp = axis;
%     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
axis_range_x = temp(2)-temp(1);
axis_range_y = temp(4)-temp(3);
percent_larger = 0.3;
axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);


subplot(1,2,2);
hold on;
grid on;
axis equal;
xlabel('X [meters]');
ylabel('Y [meters]');
title('Backward fit');

% Plot the input points very large
current_color = fcn_geometry_fillColorFromNumberOrName(1,'points',[],-1);
plot(small_XY_data(:,1),small_XY_data(:,2),'.','Color',current_color,'MarkerSize',10);

% Plot the domain points
for ith_domain = 1:length(fitSequence_points_backward)
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain,fitSequence_bestFitType{ith_domain},[],-1);
    current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain,[],[],-1);
    current_fitSequence_points = fitSequence_points_backward{ith_domain};
    current_fitSequence_shape  = fitSequence_shapes_backward{ith_domain};
    plot(current_fitSequence_points(:,1),current_fitSequence_points(:,2),'.','Color',current_color*0.8,'MarkerSize',15);
    plot(current_fitSequence_shape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',3,'EdgeAlpha',0);
end

% Plot the domain fits
fcn_geometry_plotFitSequences(fitSequence_bestFitType_backward, fitSequence_parameters_backward,(fig_num));

% Make axis slightly larger?
temp = axis;
%     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
axis_range_x = temp(2)-temp(1);
axis_range_y = temp(4)-temp(3);
percent_larger = 0.3;
axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);



%% Print the results

NleadCharacters = 20;
NfitsInSequence = length(fitSequence_bestFitType_forward);

fprintf(1,'\n\nPARAMETER FIT COMPARISON:\n');
for ith_fit = 1:NfitsInSequence
    fprintf(1,'\n\nFit Sequence Number: %.0d\n', ith_fit); 

    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('   '),NleadCharacters));
    fcn_INTERNAL_printFitDetails(fitSequence_bestFitType_forward{ith_fit}, fitSequence_parameters_forward{ith_fit},1)

    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('FORWARD'),NleadCharacters));
    fcn_INTERNAL_printFitDetails(fitSequence_bestFitType_forward{ith_fit}, fitSequence_parameters_forward{ith_fit},0)

    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('REVERSE'),NleadCharacters));
    fcn_INTERNAL_printFitDetails(fitSequence_bestFitType_backward{ith_fit}, fitSequence_parameters_backward{ith_fit},0)

    fitSequence_parameters_averaged = (fitSequence_parameters_forward{ith_fit} + fitSequence_parameters_backward{ith_fit})/2;

    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('AVERAGED'),NleadCharacters));
    fcn_INTERNAL_printFitDetails(fitSequence_bestFitType_forward{ith_fit},fitSequence_parameters_averaged,0)


    if exist('revised_fitSequence_parameters_forward','var') && ~isempty(revised_fitSequence_parameters_forward)
        fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('FORWARD REV'),NleadCharacters));
        fcn_INTERNAL_printFitDetails(fitSequence_bestFitType_forward{ith_fit},revised_fitSequence_parameters_forward{ith_fit},0)

        fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('REVERSE REV'),NleadCharacters));
        fcn_INTERNAL_printFitDetails(fitSequence_bestFitType_backward{ith_fit},revised_fitSequence_parameters_backward{ith_fit},0)

        averaged_parameters = (revised_fitSequence_parameters_backward{ith_fit} + revised_fitSequence_parameters_forward{ith_fit})/2;

        fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('AVERAGED'),NleadCharacters));
        fcn_INTERNAL_printFitDetails(fitSequence_bestFitType_forward{ith_fit},averaged_parameters,0)
    end
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
    header_strings = {'centerX','centerY','radius','startAngle(deg)','endAngle(deg)','isCircle','turnsLeft'};
elseif contains(fit_type,{'segment','line'})
    print_type = 'line';
    header_strings = {'startX','startY','theta(deg)','Slength'};
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

        % For arcs, the 4th and 5th parameter are angles that must be
        % converted into degrees. Same is true for the line type, 3rd
        % parameter.
        if strcmp(print_type,'arc') && ((ith_parameter==4)||(ith_parameter==5))
            fit_parameters(ith_parameter) = mod(fit_parameters(ith_parameter),2*pi); % Convert to positive angles only
            number_string = fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%.2f ',fit_parameters(ith_parameter)*180/pi),NumColumnChars);
        elseif strcmp(print_type,'line') && (ith_parameter==3)
            fit_parameters(ith_parameter) = mod(fit_parameters(ith_parameter),2*pi); % Convert to positive angles only
            number_string = fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%.2f ',fit_parameters(ith_parameter)*180/pi),NumColumnChars);
        else
            number_string = fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%.4f ',fit_parameters(ith_parameter)),NumColumnChars);
        end

        % If start of the number has a minus sign, end-pad it. Otherwise,
        % front-pad it. 
        if strcmp(number_string(1),'-')
            parameter_string = cat(2,parameter_string,number_string, ' ');
        else
            parameter_string = cat(2,parameter_string,' ', number_string);
        end
    end
    fprintf(1,'%s\n',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%s',parameter_string),7*NumColumnChars));
end

end % Ends fcn_INTERNAL_printFitDetails









