%% script_test_fcn_geometry_fitArcConstrainedStart.m
% Exercises the function: fcn_geometry_fitArcConstrainedStart

% Revision history:
% 2024_03_30 - S. Brennan
% -- wrote the code

close all;



%% BASIC test - perfect fit, perfectly oriented
fig_num = 1;
figure(fig_num); clf;

arc_radius = 2;
arc_center = [2 0];
angles = [0; pi/4; pi/2]; % angles relative to vertical

% Calculate locations of test points
if arc_center(1,1)<0
    arc_seed_points = arc_radius*[cos(angles) sin(angles)] + arc_center;
else
    arc_seed_points = arc_radius*[cos(pi-angles) sin(pi-angles)] + arc_center;
end

[arc_true_circleCenter, arc_true_circleRadius] = fcn_geometry_circleCenterFrom3Points(arc_seed_points(1,:),arc_seed_points(2,:),arc_seed_points(3,:),-1);

M = 20; % Number of points per meter
sigma = 0;

onearc_test_points = fcn_geometry_fillArcTestPoints(arc_seed_points, M, sigma, -1);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Find the arc fit

inputPoints = corrupted_onearc_test_points;
[fitted_radius, fitted_arcCenter, arcLength, radial_fitting_error] = fcn_geometry_fitArcConstrainedStart(inputPoints, [], [], (fig_num));

% Check sizes
assert(length(fitted_radius(:,1))==1);
assert(length(fitted_radius(1,:))==1);
assert(length(fitted_arcCenter(:,1))==1);
assert(length(fitted_arcCenter(1,:))==2);
assert(length(arcLength(:,1))==1);
assert(length(arcLength(1,:))==1);
assert(length(radial_fitting_error(:,1))==length(inputPoints(:,1)));
assert(length(radial_fitting_error(1,:))==1);

% Check values
assert(isequal(round(arc_true_circleRadius,4),round(fitted_radius,4)));
assert(isequal(round(arc_true_circleCenter,4),round(fitted_arcCenter,4)));
assert(isequal(round(angles(end)-angles(1),4),round(arcLength,4)));
assert(max(radial_fitting_error)<0.001);


%% BASIC test - perfect fit, perfectly oriented, not at origin
fig_num = 2;
figure(fig_num); clf;


arc_radius = 2;
arc_center = [2 0];
angles = [0.1; pi/4; .7*pi];

% Calculate locations of test points
if arc_center(1,1)<0
    arc_seed_points = arc_radius*[cos(angles) sin(angles)] + arc_center;
else
    arc_seed_points = arc_radius*[cos(pi-angles) sin(pi-angles)] + arc_center;
end

[arc_true_circleCenter, arc_true_circleRadius] = fcn_geometry_circleCenterFrom3Points(arc_seed_points(1,:),arc_seed_points(2,:),arc_seed_points(3,:),-1);

M = 20; % Number of points per meter
sigma = 0;

onearc_test_points = fcn_geometry_fillArcTestPoints(arc_seed_points, M, sigma, -1);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Find the arc fit

inputPoints = corrupted_onearc_test_points;
[fitted_radius, fitted_arcCenter, arcLength, radial_fitting_error] = fcn_geometry_fitArcConstrainedStart(inputPoints, [], [], (fig_num));

% Check sizes
assert(length(fitted_radius(:,1))==1);
assert(length(fitted_radius(1,:))==1);
assert(length(fitted_arcCenter(:,1))==1);
assert(length(fitted_arcCenter(1,:))==2);
assert(length(arcLength(:,1))==1);
assert(length(arcLength(1,:))==1);
assert(length(radial_fitting_error(:,1))==length(inputPoints(:,1)));
assert(length(radial_fitting_error(1,:))==1);

% Check values
assert(isequal(round(arc_true_circleRadius,4),round(fitted_radius,4)));
assert(isequal(round(arc_true_circleCenter,4),round(fitted_arcCenter,4)));
assert(isequal(round(angles(end)-angles(1),4),round(arcLength,4)));
assert(max(radial_fitting_error)<0.001);

%% BASIC test - perfect fit, perfectly oriented, negative
fig_num = 3;
figure(fig_num); clf;

arc_radius = 2;
arc_center = -[2 0];
angles = [0; pi/6; pi/4];

% Calculate locations of test points
if arc_center(1,1)<0
    arc_seed_points = arc_radius*[cos(angles) sin(angles)] + arc_center;
else
    arc_seed_points = arc_radius*[cos(pi-angles) sin(pi-angles)] + arc_center;
end

[arc_true_circleCenter, arc_true_circleRadius] = fcn_geometry_circleCenterFrom3Points(arc_seed_points(1,:),arc_seed_points(2,:),arc_seed_points(3,:),-1);

M = 20; % Number of points per meter
sigma = 0;

onearc_test_points = fcn_geometry_fillArcTestPoints(arc_seed_points, M, sigma, -1);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Find the arc fit

inputPoints = corrupted_onearc_test_points;
[fitted_radius, fitted_arcCenter, arcLength, radial_fitting_error] = fcn_geometry_fitArcConstrainedStart(inputPoints, [], [], (fig_num));

% Check sizes
assert(length(fitted_radius(:,1))==1);
assert(length(fitted_radius(1,:))==1);
assert(length(fitted_arcCenter(:,1))==1);
assert(length(fitted_arcCenter(1,:))==2);
assert(length(radial_fitting_error(:,1))==length(inputPoints(:,1)));
assert(length(radial_fitting_error(1,:))==1);

% Check values
assert(isequal(round(arc_true_circleRadius,4),round(fitted_radius,4)));
assert(isequal(round(arc_true_circleCenter,4),round(fitted_arcCenter,4)));
assert(max(radial_fitting_error)<0.001);

%% BASIC test - perfect fit, perfectly oriented, negative and noisy
fig_num = 4;
figure(fig_num); clf;
rng(4747);

arc_radius = 2;
arc_center = -[2 0];
angles = [0; pi/6; pi/4];

% Calculate locations of test points
if arc_center(1,1)<0
    arc_seed_points = arc_radius*[cos(angles) sin(angles)] + arc_center;
else
    arc_seed_points = arc_radius*[cos(pi-angles) sin(pi-angles)] + arc_center;
end

[arc_true_circleCenter, arc_true_circleRadius] = fcn_geometry_circleCenterFrom3Points(arc_seed_points(1,:),arc_seed_points(2,:),arc_seed_points(3,:),-1);

M = 20; % Number of points per meter
sigma = 0.05;

onearc_test_points = fcn_geometry_fillArcTestPoints(arc_seed_points, M, sigma, -1);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Find the arc fit

inputPoints = corrupted_onearc_test_points;
[fitted_radius, fitted_arcCenter, arcLength, radial_fitting_error] = fcn_geometry_fitArcConstrainedStart(inputPoints, [], [], (fig_num));

% Check sizes
assert(length(fitted_radius(:,1))==1);
assert(length(fitted_radius(1,:))==1);
assert(length(fitted_arcCenter(:,1))==1);
assert(length(fitted_arcCenter(1,:))==2);
assert(length(radial_fitting_error(:,1))==length(inputPoints(:,1)));
assert(length(radial_fitting_error(1,:))==1);

% Check values
assert(abs(arc_true_circleRadius - fitted_radius)<(0.001+6*sigma));
center_error = sum((arc_true_circleCenter - fitted_arcCenter).^2,2).^0.5;
assert(center_error<(0.001+6*sigma));
assert(max(radial_fitting_error)<(0.001+6*sigma));

%% BASIC test - perfect fit, perfectly oriented, negative and noisy with outliers
fig_num = 5;
figure(fig_num); clf;
%rng(4747);

arc_radius = 2;
arc_center = -[2 0];
angles = [0; pi/6; pi/2];

% Calculate locations of test points
if arc_center(1,1)<0
    arc_seed_points = arc_radius*[cos(angles) sin(angles)] + arc_center;
else
    arc_seed_points = arc_radius*[cos(pi-angles) sin(pi-angles)] + arc_center;
end

[~, ~] = fcn_geometry_circleCenterFrom3Points(arc_seed_points(1,:),arc_seed_points(2,:),arc_seed_points(3,:),-1);

M = 20; % Number of points per meter
sigma = 0.05;

onearc_test_points = fcn_geometry_fillArcTestPoints(arc_seed_points, M, sigma, -1);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Find the arc fit

inputPoints = corrupted_onearc_test_points;
[fitted_radius, fitted_arcCenter, arcLength, radial_fitting_error] = fcn_geometry_fitArcConstrainedStart(inputPoints, [], [], (fig_num));

% Check sizes
assert(length(fitted_radius(:,1))==1);
assert(length(fitted_radius(1,:))==1);
assert(length(fitted_arcCenter(:,1))==1);
assert(length(fitted_arcCenter(1,:))==2);
assert(length(radial_fitting_error(:,1))==length(inputPoints(:,1)));
assert(length(radial_fitting_error(1,:))==1);

% Check values
% assert(abs(arc_true_circleRadius - fitted_radius)<(0.001+6*sigma));
% center_error = sum((arc_true_circleCenter - fitted_arcCenter).^2,2).^0.5;
% assert(center_error<(0.001+6*sigma));
% assert(max(radial_fitting_error)<(0.001+6*sigma));

%% BASIC test - perfect fit, misaligned
fig_num = 6;
figure(fig_num); clf;

arc_radius = 2;
arc_center = [2 0];
angles = [0; pi/4; pi/3]; % angles relative to vertical

initial_rotation = pi/4;
initial_offset = [2 3];

% Calculate locations of test points
if arc_center(1,1)<0
    arc_seed_points = arc_radius*[cos(angles) sin(angles)] + arc_center;
else
    arc_seed_points = arc_radius*[cos(pi-angles) sin(pi-angles)] + arc_center;
end

% Move the seed points
arc_seed_points = arc_seed_points*[cos(initial_rotation) sin(initial_rotation); -sin(initial_rotation) cos(initial_rotation)] + initial_offset;

[arc_true_circleCenter, arc_true_circleRadius] = fcn_geometry_circleCenterFrom3Points(arc_seed_points(1,:),arc_seed_points(2,:),arc_seed_points(3,:),-1);

M = 20; % Number of points per meter
sigma = 0;

onearc_test_points = fcn_geometry_fillArcTestPoints(arc_seed_points, M, sigma, -1);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Find the arc fit

inputPoints = corrupted_onearc_test_points;
[fitted_radius, fitted_arcCenter, arcLength, radial_fitting_error] = fcn_geometry_fitArcConstrainedStart(inputPoints, initial_rotation, initial_offset, (fig_num));

% Check sizes
assert(length(fitted_radius(:,1))==1);
assert(length(fitted_radius(1,:))==1);
assert(length(fitted_arcCenter(:,1))==1);
assert(length(fitted_arcCenter(1,:))==2);
assert(length(radial_fitting_error(:,1))==length(inputPoints(:,1)));
assert(length(radial_fitting_error(1,:))==1);

% Check values
assert(isequal(round(arc_true_circleRadius,4),round(fitted_radius,4)));
assert(isequal(round(arc_true_circleCenter,4),round(fitted_arcCenter,4)));
assert(max(radial_fitting_error)<0.001);

%% BASIC test - perfect fit, perfectly oriented, negative and noisy 
fig_num = 7;
figure(fig_num); clf;
%rng(4747);

arc_radius = 2;
arc_center = -[2 0];
angles = [0; pi/6; pi/2];

initial_rotation = pi/4;
initial_offset = [2 3];

% Calculate locations of test points
if arc_center(1,1)<0
    arc_seed_points = arc_radius*[cos(angles) sin(angles)] + arc_center;
else
    arc_seed_points = arc_radius*[cos(pi-angles) sin(pi-angles)] + arc_center;
end

% Move the seed points
arc_seed_points = arc_seed_points*[cos(initial_rotation) sin(initial_rotation); -sin(initial_rotation) cos(initial_rotation)] + initial_offset;

[~, ~] = fcn_geometry_circleCenterFrom3Points(arc_seed_points(1,:),arc_seed_points(2,:),arc_seed_points(3,:),-1);

M = 20; % Number of points per meter
sigma = 0.05;

onearc_test_points = fcn_geometry_fillArcTestPoints(arc_seed_points, M, sigma, -1);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.0;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Find the arc fit

inputPoints = corrupted_onearc_test_points;
[fitted_radius, fitted_arcCenter, arcLength, radial_fitting_error] = fcn_geometry_fitArcConstrainedStart(inputPoints, initial_rotation, initial_offset, (fig_num));

% Check sizes
assert(length(fitted_radius(:,1))==1);
assert(length(fitted_radius(1,:))==1);
assert(length(fitted_arcCenter(:,1))==1);
assert(length(fitted_arcCenter(1,:))==2);
assert(length(radial_fitting_error(:,1))==length(inputPoints(:,1)));
assert(length(radial_fitting_error(1,:))==1);

%% Test of fast mode
% 
% % Fill circle data
% % circle
% circle_center = [4 3];
% circle_radius = 2;
% M = 3; % 5 points per meter
% sigma = 0.02;
% fig_num = -1;
% 
% circle_test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma); % (fig_num));
% circle_true_parameters = [circle_center, circle_radius, 0, 2*pi, 1];
% 
% % Add outliers?
% % Corrupt the results
% probability_of_corruption = 0.3;
% magnitude_of_corruption = 1;
% 
% corrupted_circle_test_points = fcn_geometry_corruptPointsWithOutliers(circle_test_points,...
%     (probability_of_corruption), (magnitude_of_corruption), (fig_num));
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
% 
% % Perform the calculation in slow mode
% fig_num = [];
% REPS = 100; minTimeSlow = Inf;
% tic;
% for i=1:REPS
%     tstart = tic;
% 
%    regression_domain  =  ...
%     fcn_geometry_fitArcRegressionFromHoughFit(domains_corrupted_circle_test_points{1}, fig_num); 
% 
% 
%     telapsed = toc(tstart);
%     minTimeSlow = min(telapsed,minTimeSlow);
% end
% averageTimeSlow = toc/REPS;
% 
% % Perform the operation in fast mode
% fig_num = -1;
% minTimeFast = Inf;
% tic;
% for i=1:REPS
%     tstart = tic;
% 
%     regression_domain  =  ...
%     fcn_geometry_fitArcRegressionFromHoughFit(domains_corrupted_circle_test_points{1}, fig_num); 
% 
%     telapsed = toc(tstart);
%     minTimeFast = min(telapsed,minTimeFast);
% end
% averageTimeFast = toc/REPS;
% 
% fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_fitArcRegressionFromHoughFit:\n');
% fprintf(1,'N repetitions: %.0d\n',REPS);
% fprintf(1,'Slow mode average speed per call (seconds): %.8f\n',averageTimeSlow);
% fprintf(1,'Slow mode fastest speed over all calls (seconds): %.8f\n',minTimeSlow);
% fprintf(1,'Fast mode average speed per call (seconds): %.8f\n',averageTimeFast);
% fprintf(1,'Fast mode fastest speed over all calls (seconds): %.8f\n',minTimeFast);
% fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
% fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);


%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
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

% function fcn_INTERNAL_printResults(param_string, params)
% header_string = sprintf('%s',param_string);
% fixed_header_string = fcn_DebugTools_debugPrintStringToNCharacters(header_string,25);
% fprintf(1,'%s ',fixed_header_string)
% for ith_value = 1:length(params)
%     param_string = sprintf('%.4f',params(ith_value));
%     fixed_param_string = fcn_DebugTools_debugPrintStringToNCharacters(param_string,15);
%     fprintf(1,'%s ',fixed_param_string)   
% end
% fprintf(1,'\n');
% end