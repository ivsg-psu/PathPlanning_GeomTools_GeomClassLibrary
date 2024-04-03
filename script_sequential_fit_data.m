%% script_sequential_fit_data
% Exercises the function: fcn_geometry_fitArcRegressionFromHoughFit

% 2024_04_02 - S. Brennan
% -- wrote the code

close all;

%% Use fillArcSequence to create some test data
fig_num = 1;
figure(fig_num);
clf;

rng(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE TEST POINTS
figure(fig_num);
subplot(2,2,1);


arc_pattern = [...
    1/20, 15; 
    0 20];

M = 10;
sigma = 0.02;

[test_points, circleCenters, trueStartPointsOfArcs, arcStartIndicies] = fcn_geometry_fillArcSequenceTestPoints(arc_pattern, M, sigma, fig_num);

% Add noise?
if 1==0
    % Corrupt the results
    probability_of_corruption = 1;
    magnitude_of_corruption = 0.03;

    test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
        (probability_of_corruption), (magnitude_of_corruption), (fig_num));
end

% Add outliers?
if 1==0
    % Corrupt the results
    probability_of_corruption = 0.1;
    magnitude_of_corruption = 1;

    test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
        (probability_of_corruption), (magnitude_of_corruption), (fig_num));
end

% Grab the axis
original_axis = axis + [-10 10 -10 10];
axis(original_axis);
axis equal;

% Label the plot
figure(fig_num);
subplot(2,2,1);
title('Input points');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform sequential segment fitting
figure(fig_num);
subplot(2,2,2);


NpointsInFit = 350;
Hough_domain.best_fit_type    = 'Hough arc';
Hough_domain.best_fit_parameters  = [nan nan nan nan nan 0]; % The zero indicates this is an arc

best_fit_domain_box_projection_distance = 0.1;
figure(fig_num);
subplot(2,2,2);
axis(original_axis);
title('Regression fit');
hold on;
grid on;
axis equal
xlabel('X [meters]');
ylabel('Y [meters]');

% Plot the original data
plot(test_points(:,1),test_points(:,2),'.','Color',[0 0 0],'MarkerSize',5);

% Plot the domain shape
Hough_domain.points_in_domain = test_points(1:NpointsInFit,1:2);
Hough_domain.best_fit_source_indicies = [1 2 NpointsInFit];
regression_domain  =  ...
    fcn_geometry_fitArcRegressionFromHoughFit(Hough_domain, best_fit_domain_box_projection_distance, -1);
domainShape = regression_domain.best_fit_domain_box;
current_color = [1 0 0];
h_domainShape = plot(domainShape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',1,'EdgeAlpha',0);

% Plot the "hit" points within the fit
empty_data = nan*test_points;
h_plot = plot(empty_data(:,1),empty_data(:,2),'.','Color',[0 1 0],'MarkerSize',10);


% Make a plot of percentage of fits
subplot(2,2,3);
hold on;
NpointsTotal = length(points_in_fit(:,1));
percentage_of_fits = nan(NpointsTotal,1);
xlabel('Number of points');
ylabel('Percentage inside');
h_percentage = plot((1:NpointsTotal)',percentage_of_fits,'k.-');
axis([0 NpointsTotal -0.1 1.1]);
grid on;

% Add vertical lines to indicate where the segments are changing
for ith_start = 1:length(arcStartIndicies)
    plot([arcStartIndicies(ith_start) arcStartIndicies(ith_start)],[-0.1 1.1],'r-');
end

Ndomains = 1;
current_segment_start_index = 1;

% Perform the fit
current_point_index = 3;
while current_point_index < length(test_points(:,1))
    current_point_index = current_point_index + 1;

    % Grab the points in current domain
    current_points_in_domain = test_points(current_segment_start_index:current_point_index,1:2);
    test_points_for_domain = test_points(current_segment_start_index:end,:);

    Hough_domain.points_in_domain = current_points_in_domain;
    Hough_domain.best_fit_source_indicies = [1 2 length(current_points_in_domain(:,1))];

    regression_domain  =  ...
        fcn_geometry_fitArcRegressionFromHoughFit(Hough_domain, best_fit_domain_box_projection_distance, -1);
    
    % Update the fit region
    set(h_domainShape,'Shape',regression_domain.best_fit_domain_box);

    % Update the points inside the fit
    buffered_box = polybuffer(regression_domain.best_fit_domain_box,0.1);
    indicies_inside_fit = isinterior(buffered_box, test_points_for_domain(:,1),test_points_for_domain(:,2));    
    points_in_fit = empty_data;
    points_in_fit(indicies_inside_fit,:) = test_points_for_domain(indicies_inside_fit,:);
    set(h_plot,'XData',points_in_fit(:,1));
    set(h_plot,'YData',points_in_fit(:,2));

    % Update the percentage plot
    NpointsInCurrentFit = length(current_points_in_domain(:,1));
    percentage_of_fits(current_point_index,1) = min(sum(indicies_inside_fit)/NpointsInCurrentFit,1);
    set(h_percentage,'YData',percentage_of_fits);

    pause(0.01)
end

% % Perform the Hough fits of data using circle fitting to find the domains
% transverse_tolerance = 0.1;
% station_tolerance = 1;
% points_required_for_agreement = 20;
% flag_force_circle_fit = 0;
% expected_radii_range = [1 10];
% flag_find_only_best_agreement = [];
% flag_use_permutations = [];
% 
% inputPoints = onearc_test_points;
% domains_onearc_test_points  = ...
% fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
%         (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement), (flag_use_permutations), (fig_num));



%% Basic call with clean data
fig_num = 1;
figure(fig_num);
clf;
hold on;

% Filling perfect test data for arcs
figure(fig_num);
subplot(2,2,1);

arc_seed_points = [2 3; 4 5; 6 3];
[arc_true_circleCenter, arc_true_circleRadius] = fcn_geometry_circleCenterFrom3Points(arc_seed_points(1,:),arc_seed_points(2,:),arc_seed_points(3,:),fig_num);

figure(fig_num);
subplot(2,2,1);
axis([1 6 1 6]);
title('Generating circle');

% Convert arc to noisy points
figure(fig_num);
subplot(2,2,2);

M = 10; % Number of points per meter
sigma = 0.02;
onearc_test_points = fcn_geometry_fillArcTestPoints(arc_seed_points, M, sigma,fig_num);

figure(fig_num);
subplot(2,2,1);
axis([1 6 1 6]);
title('Noisy points');


% % Add outliers?
% % Corrupt the results
% probability_of_corruption = 0.3;
% magnitude_of_corruption = 1;
% 
% corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
%     (probability_of_corruption), (magnitude_of_corruption), (fig_num));


% Find domain via Hough fit
figure(fig_num);
subplot(2,2,3);

% Perform the Hough fits of data using circle fitting to find the domains
transverse_tolerance = 0.1;
station_tolerance = 1;
points_required_for_agreement = 20;
flag_force_circle_fit = 0;
expected_radii_range = [1 10];
flag_find_only_best_agreement = [];
flag_use_permutations = [];


inputPoints = onearc_test_points;
domains_onearc_test_points  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement), (flag_use_permutations), (fig_num));

figure(fig_num);
subplot(2,2,3);
axis([1 6 1 6]);
title('Hough domain');

% fig_num = 235;
% figure(fig_num); clf;
% inputPoints = corrupted_twoarc_test_points;
% domains_corrupted_twoarc_test_points  = ...
% fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
%         (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement), (flag_use_permutations), (fig_num));

% % BASIC call with circle data, fitting it with a circle by not specifying station tolerance
% fig_num = 1111;
% figure(fig_num); clf;
% 
% inputPoints = corrupted_circle_test_points;
% transverse_tolerance = 0.1;
% station_tolerance = [];
% points_required_for_agreement = [];
% flag_force_circle_fit = [];
% expected_radii_range = [];
% flag_find_only_best_agreement = [];
% flag_use_permutations = [];
% 
% 
% domains_corrupted_circle_test_points  = ...
% fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
%         (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement), (flag_use_permutations), (fig_num));


figure(fig_num);
subplot(2,2,4);
regression_domain  =  ...
    fcn_geometry_fitArcRegressionFromHoughFit(domains_onearc_test_points{1}, fig_num);
figure(fig_num);
subplot(2,2,4);
axis([1 6 1 6]);
title('Regression fit');


% Print results
true_params = [arc_true_circleCenter(1,1),arc_true_circleCenter(1,2), arc_true_circleRadius, arc_true_start_angle_in_radians, arc_true_end_angle_in_radians, 0.00 ];
fprintf(1,'\n\nResults of arc regression fitting:\n')
fprintf(1,'                        [centerX          centerY         radius         startAngle      endAngle        isCircle] (meters and degrees)\n');
params = true_params;
fcn_INTERNAL_printResults('True parameters', params)
params = domains_onearc_test_points{1}.best_fit_parameters;
fcn_INTERNAL_printResults('Hough parameters', params)
params = regression_domain.best_fit_parameters;
fcn_INTERNAL_printResults('Regression fit parameters', params)
fprintf(1,'ERRORS:\n');
params = abs( domains_onearc_test_points{1}.best_fit_parameters - true_params);
fcn_INTERNAL_printResults('Hough errors', params)
params = abs(regression_domain.best_fit_parameters - true_params);
fcn_INTERNAL_printResults('Regression errors', params)