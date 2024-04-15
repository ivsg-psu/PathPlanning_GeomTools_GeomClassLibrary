% script_test_fcn_geometry_findArcAgreementIndicies
% Exercises the function: fcn_geometry_findArcAgreementIndicies

% Revision history:
% 2024_01_09
% -- wrote the code
% 2024_01_12
% -- fixed typo in one of the printouts to use correct function name.


%% Set up the workspace

close all

%% Examples for basic path operations and function testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ______                           _
% |  ____|                         | |
% | |__  __  ____ _ _ __ ___  _ __ | | ___  ___
% |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|
% | |____ >  < (_| | | | | | | |_) | |  __/\__ \
% |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/
%                            | |
%                            |_|
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% %% Fill test data 
% % fig_num = 21;
% % figure(fig_num);
% % clf;
% % hold on;
% % axis equal
% % grid on;
% 
% % circle
% circle_center = [4 3];
% circle_radius = 2;
% M = 5; % 5 points per meter
% sigma = 0.02;
% 
% circle_test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma); % (fig_num));
% 
% 
% % Add outliers?
% % Corrupt the results
% probability_of_corruption = 0.3;
% magnitude_of_corruption = 1;
% 
% corrupted_circle_test_points = fcn_geometry_corruptPointsWithOutliers(circle_test_points,...
%     (probability_of_corruption), (magnitude_of_corruption)); %, (fig_num));
% 
% 
% % 1 arc
% % fig_num = 23;
% % figure(fig_num);
% % clf;
% % hold on;
% % axis equal
% % grid on;
% 
% seed_points = [2 3; 4 5; 6 3];
% [true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
% trueParameters = [true_circleCenter true_circleRadius];
% 
% M = 10; % Number of points per meter
% sigma = 0.02;
% 
% onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);
% 
% % Add outliers?
% % Corrupt the results
% probability_of_corruption = 0.3;
% magnitude_of_corruption = 1;
% 
% corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
%     (probability_of_corruption), (magnitude_of_corruption)); %, (fig_num));
% 
% % Fill test data - 2 arcs
% twoarc_test_points = [onearc_test_points(1:30,:); onearc_test_points(50:60,:)];
% corrupted_twoarc_test_points = [corrupted_onearc_test_points(1:30,:); corrupted_onearc_test_points(50:60,:)];

%% Example - 1 - BASIC call with arc data, fitting it to source point 10
fig_num = 111;
figure(fig_num); clf;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
% trueParameters = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);


% Fill test data - 2 arcs
twoarc_test_points = [onearc_test_points(1:30,:); onearc_test_points(50:60,:)];

points = twoarc_test_points;
circleCenter = true_circleCenter;
circleRadius = true_circleRadius;
index_source_point = 10;
station_tolerance = 0.5;

[indicies_in_station_agreement, flag_is_a_circle, start_angle_in_radians, end_angle_in_radians] = ...
    fcn_geometry_findArcAgreementIndicies(points, circleCenter, circleRadius, index_source_point, station_tolerance, fig_num);

%% Example - 2 - BASIC call with arc data, fitting it to source point 55
fig_num = 222;
figure(fig_num); clf;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);


% Fill test data - 2 arcs
twoarc_test_points = [onearc_test_points(1:30,:); onearc_test_points(50:60,:)];

points = twoarc_test_points;
circleCenter = true_circleCenter;
circleRadius = true_circleRadius;
index_source_point = 40;
station_tolerance = 0.5;

[indicies_in_station_agreement, flag_is_a_circle, start_angle_in_radians, end_angle_in_radians] = ...
    fcn_geometry_findArcAgreementIndicies(points, circleCenter, circleRadius, index_source_point, station_tolerance, fig_num);

%% Example - 3 - BASIC call with arc data, fitting it to source point 10
% But with a larger tolerance that will bridge the arcs
fig_num = 333;
figure(fig_num); clf;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);


% Fill test data - 2 arcs
twoarc_test_points = [onearc_test_points(1:30,:); onearc_test_points(50:60,:)];

points = twoarc_test_points;
circleCenter = true_circleCenter;
circleRadius = true_circleRadius;
index_source_point = 10;
station_tolerance = 2;

[indicies_in_station_agreement, flag_is_a_circle, start_angle_in_radians, end_angle_in_radians] = ...
    fcn_geometry_findArcAgreementIndicies(points, circleCenter, circleRadius, index_source_point, station_tolerance, fig_num); %#ok<*ASGLU>

%% Example - 1 - BASIC call with circle data, fitting it with a circle
fig_num = 444;
figure(fig_num); clf;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Fill test data - 2 arcs
twoarc_test_points = [onearc_test_points(1:30,:); onearc_test_points(50:60,:)];

points = twoarc_test_points;
circleCenter = true_circleCenter;
circleRadius = true_circleRadius;
index_source_point = 10;
station_tolerance = 10;

[indicies_in_station_agreement, flag_is_a_circle, start_angle_in_radians, end_angle_in_radians] = ...
    fcn_geometry_findArcAgreementIndicies(points, circleCenter, circleRadius, index_source_point, station_tolerance, fig_num);



%% Test of fast implementation mode 

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Fill test data - 2 arcs
twoarc_test_points = [onearc_test_points(1:30,:); onearc_test_points(50:60,:)];

points = twoarc_test_points;
circleCenter = true_circleCenter;
circleRadius = true_circleRadius;
index_source_point = 10;
station_tolerance = 2;

% Perform the calculation in slow mode
REPS = 100; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    
    [indicies_in_station_agreement, flag_is_a_circle, start_angle_in_radians, end_angle_in_radians] = fcn_geometry_findArcAgreementIndicies(points, circleCenter, circleRadius, index_source_point, station_tolerance, []);

    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;

    [indicies_in_station_agreement, flag_is_a_circle, start_angle_in_radians, end_angle_in_radians] = fcn_geometry_findArcAgreementIndicies(points, circleCenter, circleRadius, index_source_point, station_tolerance, -1);
    
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_findArcAgreementIndicies:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

%% ERROR Example - 1 - Show how it works with corrupted data (it still works)
fig_num = 999;
figure(fig_num); clf;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption)); %, (fig_num));

% Fill test data - 2 arcs
twoarc_test_points = [onearc_test_points(1:30,:); onearc_test_points(50:60,:)];
corrupted_twoarc_test_points = [corrupted_onearc_test_points(1:30,:); corrupted_onearc_test_points(50:60,:)];



points = corrupted_twoarc_test_points;
circleCenter = true_circleCenter;
circleRadius = true_circleRadius;
index_source_point = 10;
station_tolerance = 0.5;

[indicies_in_station_agreement, flag_is_a_circle, start_angle_in_radians, end_angle_in_radians] = ...
    fcn_geometry_findArcAgreementIndicies(points, circleCenter, circleRadius, index_source_point, station_tolerance, fig_num);