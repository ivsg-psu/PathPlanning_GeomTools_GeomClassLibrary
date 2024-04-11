% script_test_fcn_geometry_HoughSegmentation
% Exercises the function: fcn_geometry_HoughSegmentation

% Revision history:
% 2023_12_15
% -- wrote the code

close all;

%% Fill in test points
% Single segment
seed_points = [5 0; 15 10];
M = 5; % 10 points per meter
sigma = 0.;
rng(3423)

probability_of_corruption = 0.1;
magnitude_of_corruption = 3;


single_segment_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma,-1);
corrupted_single_segment_test_points = fcn_geometry_corruptPointsWithOutliers(single_segment_test_points,...
    (probability_of_corruption), (magnitude_of_corruption),-1);


% Multiple segments
seed_points = [2 3; 4 5; 7 0; 9 5; 10 20; 13 14];
M = 5; % 10 points per meter
sigma = 0.;
rng(3423)

multi_segment_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma,-1);
corrupted_multi_segment_test_points = fcn_geometry_corruptPointsWithOutliers(multi_segment_test_points,...
    (probability_of_corruption), (magnitude_of_corruption),-1);

% Create circle data
circle_center = [3 5];
circle_radius = 2;
M = 5; % 5 points per meter
sigma = 0.02;

circle_test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma,-1); % (fig_num));
corrupted_circle_test_points = fcn_geometry_corruptPointsWithOutliers(circle_test_points,...
    (probability_of_corruption), (magnitude_of_corruption),-1);

% Seed test data for arcs
arc_seed_points = [2 3; 4 5; 6 3];
[arc_true_circleCenter, arc_true_circleRadius] = fcn_geometry_circleCenterFrom3Points(arc_seed_points(1,:),arc_seed_points(2,:),arc_seed_points(3,:),-1);

M = 10; % Points per meter
sigma = 0.02;

% Fill test data for 1 arc
onearc_test_points = fcn_geometry_fillArcTestPoints(arc_seed_points, M, sigma); %, fig_num);
corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

% Fill test data for 2 arcs
first_fraction = [0 0.5]; % data from 0 to 50 percent
second_fraction = [0.80 1]; % data from 80 percent to end
N_points = length(onearc_test_points(:,1));

first_fraction_indicies = round(first_fraction*N_points); % find closest indicies
first_fraction_indicies = max([first_fraction_indicies; 1 1],[],1); % Make sure none are below 1
first_fraction_indicies = min([first_fraction_indicies; N_points N_points],[],1); % Make sure none are above N_points

second_fraction_indicies = round(second_fraction*N_points); % find closest indicies
second_fraction_indicies = max([second_fraction_indicies; 1 1],[],1); % Make sure none are below 1
second_fraction_indicies = min([second_fraction_indicies; N_points N_points],[],1); % Make sure none are above N_points

twoarc_test_points = ...
    [onearc_test_points(first_fraction_indicies(1):first_fraction_indicies(2),:); ...
    onearc_test_points(second_fraction_indicies(1):second_fraction_indicies(2),:)];

corrupted_twoarc_test_points = ...
    [corrupted_onearc_test_points(first_fraction_indicies(1):first_fraction_indicies(2),:); ...
    corrupted_onearc_test_points(second_fraction_indicies(1):second_fraction_indicies(2),:)];


% % For debugging
% figure(33838);
% plot(corrupted_twoarc_test_points(:,1),corrupted_twoarc_test_points(:,2),'k.');

%% Basic example: find one line segment


seed_points = [5 0; 15 10];
M = 5; % 10 points per meter
sigma = 0.;
rng(3423)

probability_of_corruption = 0.1;
magnitude_of_corruption = 3;


single_segment_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma,-1);
corrupted_single_segment_test_points = fcn_geometry_corruptPointsWithOutliers(single_segment_test_points,...
    (probability_of_corruption), (magnitude_of_corruption),-1);

fig_num = 10; 
figure(fig_num); clf;

transverse_tolerance = 0.05; % Units are meters
station_tolerance = 1; % Units are meters. Usually station tolerance needs to be larger than transverse tolerance, and it needs to be large enough that it can span gaps in corrupted data
threshold_max_points = 10;
input_points = corrupted_single_segment_test_points;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

%% Basic example: find 5 lines within noisy data

% Multiple segments
seed_points = [2 3; 4 5; 7 0; 9 5; 10 20; 13 14];
M = 5; % 10 points per meter
sigma = 0.;
rng(3423)

probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

multi_segment_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma,-1);
corrupted_multi_segment_test_points = fcn_geometry_corruptPointsWithOutliers(multi_segment_test_points,...
    (probability_of_corruption), (magnitude_of_corruption),-1);

fig_num = 1;
transverse_tolerance = 0.1; % Units are meters
station_tolerance = inf; % Units are meters
threshold_max_points = 20;
input_points = corrupted_multi_segment_test_points;
domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

%% Basic example 2: find 5 line segments within same data
% segments are created by imposing constraints on separation

% Multiple segments
seed_points = [2 3; 4 5; 7 0; 9 5; 10 20; 13 14];
M = 5; % 10 points per meter
sigma = 0.;
rng(3423)

probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

multi_segment_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma,-1);
corrupted_multi_segment_test_points = fcn_geometry_corruptPointsWithOutliers(multi_segment_test_points,...
    (probability_of_corruption), (magnitude_of_corruption),-1);


% Call the segmentation function
fig_num = 2;
transverse_tolerance = 0.1; % Units are meters
station_tolerance = 1; % Units are meters. Usually station tolerance needs to be larger than transverse tolerance, and it needs to be large enough that it can span gaps in corrupted data
threshold_max_points = 10;
input_points = corrupted_multi_segment_test_points;
domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

%% Basic example 3: find circle data

% Create circle data
circle_center = [3 5];
circle_radius = 2;
M = 5; % 5 points per meter
sigma = 0.02;
rng(3423)

probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

circle_test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma,-1); % (fig_num));
corrupted_circle_test_points = fcn_geometry_corruptPointsWithOutliers(circle_test_points,...
    (probability_of_corruption), (magnitude_of_corruption),-1);

fig_num = 3;
transverse_tolerance = 0.1; % Units are meters
station_tolerance = 2; % Units are meters. Usually station tolerance needs to be larger than transverse tolerance, and it needs to be large enough that it can span gaps in corrupted data
threshold_max_points = 20;
input_points = corrupted_circle_test_points;
domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

%% Basic example 3: find arc data

% Seed test data for arcs
arc_seed_points = [2 3; 4 5; 6 3];
[arc_true_circleCenter, arc_true_circleRadius] = fcn_geometry_circleCenterFrom3Points(arc_seed_points(1,:),arc_seed_points(2,:),arc_seed_points(3,:),-1);

M = 10; % Points per meter
sigma = 0.02;

rng(3423)

probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

% Fill test data for 1 arc
onearc_test_points = fcn_geometry_fillArcTestPoints(arc_seed_points, M, sigma); %, fig_num);
corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));

fig_num = 4;
transverse_tolerance = 0.03; % Units are meters
station_tolerance = 1; % Units are meters. Usually station tolerance needs to be larger than transverse tolerance, and it needs to be large enough that it can span gaps in corrupted data
threshold_max_points = 20;
input_points = corrupted_onearc_test_points;
domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

%% Advanced example: find line segments and circles in same data set

fig_num = 11; 
figure(fig_num); clf;

% Create circle data
circle_center = [3 5];
circle_radius = 2;
M = 5; % 5 points per meter
sigma = 0.02;
rng(3423)

probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

circle_test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma,-1); % (fig_num));
corrupted_circle_test_points = fcn_geometry_corruptPointsWithOutliers(circle_test_points,...
    (probability_of_corruption), (magnitude_of_corruption),-1);

seed_points = [5 0; 15 10];
M = 5; % 10 points per meter
sigma = 0.;
rng(3423)

single_segment_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma,-1);
corrupted_single_segment_test_points = fcn_geometry_corruptPointsWithOutliers(single_segment_test_points,...
    (probability_of_corruption), (magnitude_of_corruption),-1);


transverse_tolerance = 0.05; % Units are meters
station_tolerance = 3; % Units are meters. Usually station tolerance needs to be larger than transverse tolerance, and it needs to be large enough that it can span gaps in corrupted data
threshold_max_points = 10;
input_points = [corrupted_circle_test_points; corrupted_single_segment_test_points];
domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);


%% Advanced example: find line segments and circles and arcs in same data set
fig_num = 22; 
figure(fig_num); clf;

fig_num = 11; 
figure(fig_num); clf;

% Create circle data
circle_center = [3 5];
circle_radius = 2;
M = 5; % 5 points per meter
sigma = 0.02;
rng(3423)

probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

circle_test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma,-1); % (fig_num));
corrupted_circle_test_points = fcn_geometry_corruptPointsWithOutliers(circle_test_points,...
    (probability_of_corruption), (magnitude_of_corruption),-1);

seed_points = [5 0; 15 10];
M = 5; % 10 points per meter
sigma = 0.;
rng(3423)

single_segment_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma,-1);
corrupted_single_segment_test_points = fcn_geometry_corruptPointsWithOutliers(single_segment_test_points,...
    (probability_of_corruption), (magnitude_of_corruption),-1);

% Seed test data for arcs
arc_seed_points = [2 3; 4 5; 6 3];
[arc_true_circleCenter, arc_true_circleRadius] = fcn_geometry_circleCenterFrom3Points(arc_seed_points(1,:),arc_seed_points(2,:),arc_seed_points(3,:),-1);

M = 10; % Points per meter
sigma = 0.02;

rng(3423)

probability_of_corruption = 0.1;
magnitude_of_corruption = 3;

% Fill test data for 1 arc
onearc_test_points = fcn_geometry_fillArcTestPoints(arc_seed_points, M, sigma); %, fig_num);
corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (-1));



transverse_tolerance = 0.05; % Units are meters
station_tolerance = 1.5; % Units are meters. Usually station tolerance needs to be larger than transverse tolerance, and it needs to be large enough that it can span gaps in corrupted data
threshold_max_points = 20;
input_points = [corrupted_circle_test_points; corrupted_single_segment_test_points; corrupted_onearc_test_points+[0 8]];
domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);


%% Advanced example 3: find segments within a chevron
M = 10; % points per meter

rng(234)
sigma = 0.02;

multi_segment_test_points = [];

seed_points = [0 0; 10 0];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
multi_segment_test_points = [multi_segment_test_points; subtest_points]; %#ok<*NASGU>

seed_points = [0 0; 10 5];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
multi_segment_test_points = [multi_segment_test_points; subtest_points];

seed_points = [2 0; 3 1.5];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
multi_segment_test_points = [multi_segment_test_points; subtest_points];

seed_points = [4 0; 5 2.5];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
multi_segment_test_points = [multi_segment_test_points; subtest_points];

seed_points = [6 0; 7 3.5];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
multi_segment_test_points = [multi_segment_test_points; subtest_points];

seed_points = [8 0; 9 4.5];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
multi_segment_test_points = [multi_segment_test_points; subtest_points];

seed_points = [10 0; 10 5];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
multi_segment_test_points = [multi_segment_test_points; subtest_points];

% Add outliers to corrupt the results
outliers = [10*rand(100,1) 5*rand(100,1)];
multi_segment_test_points = [multi_segment_test_points; outliers];


% Call the segmentation function
fig_num = 3;
transverse_tolerance = 0.1; % Units are meters
station_tolerance = 0.6; % Units are meters
threshold_max_points = 30;
input_points = multi_segment_test_points;
domains = fcn_geometry_HoughSegmentation(multi_segment_test_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

fig_num = 338;
figure(fig_num); clf;
fcn_geometry_plotFitDomains(domains, fig_num);

%% Advanced example 3: find segments within a hashtag
M = 10; % 40 points per meter

rng(234)
sigma = 0.02;

multi_segment_test_points = [];

seed_points = [2 0; 2 6];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
multi_segment_test_points = [multi_segment_test_points; subtest_points]; %#ok<*NASGU>

seed_points = [4 0; 4 6];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
multi_segment_test_points = [multi_segment_test_points; subtest_points];

seed_points = [0 2; 6 2];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
multi_segment_test_points = [multi_segment_test_points; subtest_points];

seed_points = [0 4; 6 4];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
multi_segment_test_points = [multi_segment_test_points; subtest_points];

% Add outliers to corrupt the results
N_outliers = 30;
outliers = [6*rand(N_outliers,1) 6*rand(N_outliers,1)];
multi_segment_test_points = [multi_segment_test_points; outliers];


% Call the segmentation function
fig_num = 4;
transverse_tolerance = 0.1; % Units are meters
station_tolerance = 0.6; % Units are meters
threshold_max_points = 30;
input_points = multi_segment_test_points;
domains = fcn_geometry_HoughSegmentation(multi_segment_test_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

fig_num = 383;
figure(fig_num); clf;
fcn_geometry_plotFitDomains(domains, fig_num);

%% Test of fast mode
% Perform the calculation in slow mode

M = 10; % 40 points per meter

rng(234)
sigma = 0.02;

multi_segment_test_points = [];

seed_points = [2 0; 2 6];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
multi_segment_test_points = [multi_segment_test_points; subtest_points]; %#ok<*NASGU>

seed_points = [4 0; 4 6];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
multi_segment_test_points = [multi_segment_test_points; subtest_points];

seed_points = [0 2; 6 2];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
multi_segment_test_points = [multi_segment_test_points; subtest_points];

seed_points = [0 4; 6 4];
subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
multi_segment_test_points = [multi_segment_test_points; subtest_points];

% Add outliers to corrupt the results
N_outliers = 30;
outliers = [6*rand(N_outliers,1) 6*rand(N_outliers,1)];
multi_segment_test_points = [multi_segment_test_points; outliers];

% Call the segmentation function
fig_num = 4;
transverse_tolerance = 0.1; % Units are meters
station_tolerance = 0.6; % Units are meters
threshold_max_points = 30;
input_points = multi_segment_test_points;



fig_num = [];
REPS = 1; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    domains = fcn_geometry_HoughSegmentation(multi_segment_test_points, threshold_max_points, transverse_tolerance, station_tolerance, (fig_num));
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
fig_num = -1;
minTimeFast = Inf; 
tic;
for i=1:REPS
    tstart = tic;

    domains = fcn_geometry_HoughSegmentation(multi_segment_test_points, threshold_max_points, transverse_tolerance, station_tolerance, (fig_num));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'Comparison of fast and slow modes of fcn_geometry_HoughSegmentation:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.8f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.8f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.8f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.8f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);



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

