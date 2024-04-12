% script_test_fcn_geometry_fitHoughCircle
% Exercises the function: fcn_geometry_fitHoughCircle

% Revision history:
% 2023_12_15
% -- wrote the code
% 2024_01_04
% -- fixed the argument inputs
% 2024_01_17
% -- changed output types to domains

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

%% Fill test data 
fig_num = 21;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

% circle
circle_center = [4 3];
circle_radius = 2;
M = 3; % 5 points per meter
sigma = 0.02;

circle_test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma); % (fig_num));


% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_circle_test_points = fcn_geometry_corruptPointsWithOutliers(circle_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));


% 1 arc
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters_onearc_test_points = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

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


% For debugging
figure(33838);
clf;
hold on;
grid on;
grid minor;
axis equal;
plot(corrupted_twoarc_test_points(:,1),corrupted_twoarc_test_points(:,2),'k.');

% 1 outlier arc
seed_points = [6 6; 9 3; 6 0];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters = [true_circleCenter true_circleRadius];

M = 8; % Number of points per meter
sigma = 0.02;

outlieronearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
% For debugging
figure(234);
clf;
hold on;
grid on;
grid minor;
axis equal;

probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_outlieronearc_test_points= fcn_geometry_corruptPointsWithOutliers(outlieronearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (234));


%% BASIC call with arc data, fitting it with a circle by not specifying station tolerance

% 1 arc
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters_onearc_test_points = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

fig_num = 111;
figure(fig_num); clf;

inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.1;
station_tolerance = [];
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [];
flag_find_only_best_agreement = []; flag_use_permutations = [];


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement), (flag_use_permutations), (fig_num));

%% BASIC call with arc data, fitting it with an arc by specifying low station tolerance

% 1 arc
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters_onearc_test_points = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));


fig_num = 222;
figure(fig_num); clf;

inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.1;
station_tolerance = 1;
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [];
flag_find_only_best_agreement = []; flag_use_permutations = [];


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range),(flag_find_only_best_agreement),(flag_use_permutations), (fig_num));

%% BASIC call with arc data, fitting it with a circle and some arcs by specifying medium station tolerance

% 1 arc
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters_onearc_test_points = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));


fig_num = 333;
figure(fig_num); clf;

inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.1;
station_tolerance = 4;
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [];
flag_find_only_best_agreement = []; flag_use_permutations = [];


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),(flag_use_permutations), (fig_num));

%% BASIC call with arc data, fitting it with all circles by specifying large station tolerance

% 1 arc
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters_onearc_test_points = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));


fig_num = 334;
figure(fig_num); clf;

inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.1;
station_tolerance = 10;
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [];
flag_find_only_best_agreement = []; flag_use_permutations = [];


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),(flag_use_permutations), (fig_num));

%% BASIC call with arc data, forcing circle fit to fail, by using flag_force_circle_fit to give ONLY circles that meet all criteria

% 1 arc
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters_onearc_test_points = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));


fig_num = 444;
figure(fig_num); clf;

inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.1;
station_tolerance = 1;
points_required_for_agreement = [];
flag_force_circle_fit = 1;
expected_radii_range = [];
flag_find_only_best_agreement = []; flag_use_permutations = [];


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),(flag_use_permutations), (fig_num));

%% BASIC call with arc data, forcing circle fit by using flag_force_circle_fit to give ONLY circles that meet all criteria

% 1 arc
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters_onearc_test_points = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));


% Passes because station tolerance is larger
fig_num = 445;
figure(fig_num); clf;

inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.1;
station_tolerance = 4;
points_required_for_agreement = [];
flag_force_circle_fit = 1;
expected_radii_range = [];
flag_find_only_best_agreement = []; flag_use_permutations = [];


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));

%% BASIC call with circle data, fitting it with a circle by not specifying station tolerance

fig_num = 21;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

% circle
circle_center = [4 3];
circle_radius = 2;
M = 3; % 5 points per meter
sigma = 0.02;

circle_test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma); % (fig_num));


% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_circle_test_points = fcn_geometry_corruptPointsWithOutliers(circle_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

fig_num = 1111;
figure(fig_num); clf;

inputPoints = corrupted_circle_test_points;
transverse_tolerance = 0.1;
station_tolerance = [];
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [];
flag_find_only_best_agreement = []; flag_use_permutations = [];


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));

%% BASIC call with circle data, fitting it with a circle by specifying station tolerance that winds all the way around

fig_num = 21;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

% circle
circle_center = [4 3];
circle_radius = 2;
M = 3; % 5 points per meter
sigma = 0.02;

circle_test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma); % (fig_num));


% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_circle_test_points = fcn_geometry_corruptPointsWithOutliers(circle_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));


% NOTE: notice how this is much slower than the previous call, as it takes
% siginificant computation to check arcs, which is required when
% station_tolerance is given
fig_num = 2222;
figure(fig_num); clf;

inputPoints = corrupted_circle_test_points;
transverse_tolerance = 0.1;
station_tolerance = 2;
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [];
flag_find_only_best_agreement = []; flag_use_permutations = [];


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));

%% BASIC call with circle data, fitting it with an arc by specifying station tolerance is too small

fig_num = 21;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

% circle
circle_center = [4 3];
circle_radius = 2;
M = 3; % 5 points per meter
sigma = 0.02;

circle_test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma); % (fig_num));


% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_circle_test_points = fcn_geometry_corruptPointsWithOutliers(circle_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));


fig_num = 3333;
figure(fig_num); clf;

inputPoints = corrupted_circle_test_points;
transverse_tolerance = 0.1;
station_tolerance = 0.8;
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [];
flag_find_only_best_agreement = []; flag_use_permutations = [];


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));

%% BASIC call with circle AND arc data, fitting it with an arc because this has the most points

fig_num = 21;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

% circle
circle_center = [4 3];
circle_radius = 2;
M = 3; % 5 points per meter
sigma = 0.02;

circle_test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma); % (fig_num));


% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_circle_test_points = fcn_geometry_corruptPointsWithOutliers(circle_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));


% 1 outlier arc
seed_points = [6 6; 9 3; 6 0];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters = [true_circleCenter true_circleRadius];

M = 8; % Number of points per meter
sigma = 0.02;

outlieronearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
% For debugging
figure(234);
clf;
hold on;
grid on;
grid minor;
axis equal;

probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_outlieronearc_test_points= fcn_geometry_corruptPointsWithOutliers(outlieronearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (234));


fig_num = 4444;
figure(fig_num); clf;

inputPoints = [corrupted_circle_test_points; corrupted_outlieronearc_test_points];
% inputPoints = [circle_test_points; outlieronearc_test_points];
% inputPoints = [outlieronearc_test_points];


transverse_tolerance = 0.1;
station_tolerance = 2;
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [];
flag_find_only_best_agreement = []; flag_use_permutations = [];


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));

%% BASIC call with circle AND arc data, fitting it with a circle because of flag setting

fig_num = 21;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

% circle
circle_center = [4 3];
circle_radius = 2;
M = 3; % 5 points per meter
sigma = 0.02;

circle_test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma); % (fig_num));


% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_circle_test_points = fcn_geometry_corruptPointsWithOutliers(circle_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));


% 1 outlier arc
seed_points = [6 6; 9 3; 6 0];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters = [true_circleCenter true_circleRadius];

M = 8; % Number of points per meter
sigma = 0.02;

outlieronearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
% For debugging
figure(234);
clf;
hold on;
grid on;
grid minor;
axis equal;

probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_outlieronearc_test_points= fcn_geometry_corruptPointsWithOutliers(outlieronearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (234));


fig_num = 5555;
figure(fig_num); clf;

inputPoints = [corrupted_circle_test_points; corrupted_outlieronearc_test_points];
% inputPoints = [circle_test_points; outlieronearc_test_points];
% inputPoints = [outlieronearc_test_points];


transverse_tolerance = 0.1;
station_tolerance = 1;
points_required_for_agreement = [];
flag_force_circle_fit = 1;
expected_radii_range = [];
flag_find_only_best_agreement = []; flag_use_permutations = [];


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));

%% Test using expected radii range

% 1 arc
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters_onearc_test_points = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));


fig_num = 2;
figure(fig_num); clf;

inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.1;
station_tolerance = [];
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [1 3];
flag_find_only_best_agreement = []; flag_use_permutations = [];


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));

%% Test using expected radii range that is bad

% 1 arc
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters_onearc_test_points = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));


fig_num = 2;
figure(fig_num); clf;

inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.1;
station_tolerance = [];
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [10 30];
flag_find_only_best_agreement = []; flag_use_permutations = [];


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));

%% Speed test effect of adding radii range to show this speeds up calculations

% 1 arc
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters_onearc_test_points = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));


% Note, there is more speed-up the more corrupted and larger the data is
% used
inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.1;
station_tolerance = [];
points_required_for_agreement = [];
flag_force_circle_fit = [];
flag_find_only_best_agreement = []; flag_use_permutations = [];

% Perform the calculation in slow mode
expected_radii_range = [];
fig_num = -1;
REPS = 1; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;

    domains  = ...
        fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));

    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
expected_radii_range = [1 3];
fig_num = -1;
minTimeFast = Inf;
tic;
for i=1:REPS
    tstart = tic;

    domains  = ...
        fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));

    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fcn_geometry_fitHoughCircle without radii range (slow) and with radii range (fast):\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

%% Test arc fitting using station_tolerance


% 1 arc
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters_onearc_test_points = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

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

corrupted_twoarc_test_points = ...
    [corrupted_onearc_test_points(first_fraction_indicies(1):first_fraction_indicies(2),:); ...
    corrupted_onearc_test_points(second_fraction_indicies(1):second_fraction_indicies(2),:)];


inputPoints = corrupted_twoarc_test_points;
transverse_tolerance = 0.1;
station_tolerance = [];
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [1 5];
flag_find_only_best_agreement = []; 
flag_use_permutations = [];


% Use station tolerance low to find only largest arc
station_tolerance = 0.3;
fig_num = 7777;
figure(fig_num); clf;


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));

% Make station tolerance larger so it finds entire arc, connecting together
% but not finding a circle
station_tolerance = 3;
fig_num = 7788;
figure(fig_num); clf;


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));

% Fit a circle by shutting station tolerance off
station_tolerance = [];
fig_num = 7799;
figure(fig_num); clf;


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));

% Force a circle that forces station tolerance to be met by using flag
station_tolerance = 10;
flag_force_circle_fit = 1;
fig_num = 7766;
figure(fig_num); clf;


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));



%% Test using flag_use_permutations to search only points in sequence
% This forces the search to assume the points are ordered. This speeds up
% the search process quite a bit, since much of the search is simply
% sorting points. This will not work well if the points are not well
% sorted, for example if there are a large number of outliers.

% 1 arc
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters_onearc_test_points = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));


fig_num = 44;
figure(fig_num); clf;

inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.3;
station_tolerance = 0.5;
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [1 3];
flag_find_only_best_agreement = []; 

flag_use_permutations = 0;


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));

%% Test using flag_use_permutations for fractional setting (e.g. only 50%)
% This assumes that the points are over-fitted, e.g. that there are way
% more points than needed to fit. By setting a fraction, we can specify to
% only use a fraction of the input data, which speeds things up.


% 1 arc
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters_onearc_test_points = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));


fig_num = 444;
figure(fig_num); clf;

inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.3;
station_tolerance = 0.5;
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [1 3];
flag_find_only_best_agreement = []; 

flag_use_permutations = 0.5;


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));

fprintf(1,'\n\nTrue parameters (X Y radius, all in meters): %.3f %.3f %.3f\n',trueParameters_onearc_test_points(1),trueParameters_onearc_test_points(2),trueParameters_onearc_test_points(3));
fprintf(1,'Results of flag_use_permutations set to: %.5f\n',flag_use_permutations);
% fprintf(1,'Fit parameters (X Y radius, all in meters): %.3f %.3f %.3f\n',domains{1}.best_fit_parameters(1),domains{1}.best_fit_parameters(2),domains{1}.best_fit_parameters(3));

%% Test using flag_use_permutations for N-point setting
% This specifies the number of points to use as N points. This should only
% be used if the data is VERY clean, wherein all the data is quite good.
% However, if this is the case, this is quite fast.


% 1 arc
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters_onearc_test_points = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));


fig_num = 55;
figure(fig_num); clf;

inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.3;
station_tolerance = 0.5;
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [1 3];
flag_find_only_best_agreement = []; 

flag_use_permutations = 30;


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));

fprintf(1,'\n\nTrue parameters (X Y radius, all in meters): %.3f %.3f %.3f\n',trueParameters_onearc_test_points(1),trueParameters_onearc_test_points(2),trueParameters_onearc_test_points(3));
fprintf(1,'Results of flag_use_permutations set to: %.5f\n',flag_use_permutations);
% fprintf(1,'Fit parameters (X Y radius, all in meters): %.3f %.3f %.3f\n',domains{1}.best_fit_parameters(1),domains{1}.best_fit_parameters(2),domains{1}.best_fit_parameters(3));

% Now show how it works with "clean" points
fig_num = 66;
figure(fig_num); clf;

inputPoints = onearc_test_points;
transverse_tolerance = 0.3;
station_tolerance = 0.5;
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [1 3];
flag_find_only_best_agreement = []; 

flag_use_permutations = 30;


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));

fprintf(1,'\n\nTrue parameters (X Y radius, all in meters): %.3f %.3f %.3f\n',trueParameters_onearc_test_points(1),trueParameters_onearc_test_points(2),trueParameters_onearc_test_points(3));
fprintf(1,'Results of flag_use_permutations set to: %.5f\n',flag_use_permutations);
fprintf(1,'Fit parameters (X Y radius, all in meters): %.3f %.3f %.3f\n',domains{1}.best_fit_parameters(1),domains{1}.best_fit_parameters(2),domains{1}.best_fit_parameters(3));





%% Speed test effect of flag_use_permutations
% 1 arc
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters_onearc_test_points = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.1;
station_tolerance = 1;
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [1 3];
flag_find_only_best_agreement = []; 

slow_flag_use_permutations = 1;
fast_flag_use_permutations = 0;



% Perform the calculation in slow mode
fig_num = [];
REPS = 3; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;

    slowdomains  = ...
        fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement), (slow_flag_use_permutations), (fig_num));

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

    fastdomains  = ...
        fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (fast_flag_use_permutations), (fig_num));
    
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fcn_geometry_fitHoughCircle with flag_use_permutations = %.2d (slow) and with flag_use_permutations=%.2f (fast):\n',slow_flag_use_permutations,fast_flag_use_permutations);
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'True parameters (X Y radius, all in meters): %.3f %.3f %.3f\n',trueParameters_onearc_test_points(1),trueParameters_onearc_test_points(2),trueParameters_onearc_test_points(3));
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
% fprintf(1,'Slow parameters (X Y radius, all in meters): %.3f %.3f %.3f\n',slowdomains{1}.best_fit_parameters(1),slowdomains{1}.best_fit_parameters(2),slowdomains{1}.best_fit_parameters(3));
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
% fprintf(1,'Fast parameters (X Y radius, all in meters): %.3f %.3f %.3f\n',fastdomains{1}.best_fit_parameters(1),fastdomains{1}.best_fit_parameters(2),fastdomains{1}.best_fit_parameters(3));
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

%% Test of fast implementation mode 


% 1 arc
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters_onearc_test_points = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));


inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.1;
station_tolerance = 1;
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [1 3];
flag_find_only_best_agreement = []; flag_use_permutations = [];

% Perform the calculation in slow mode
fig_num = [];
REPS = 3; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;

    domains  = ...
        fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));
    
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

    domains  = ...
        fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));

    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_fitHoughCircle:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

%% Basic test illustrating flag_find_only_best_agreement

inputPoints = [
    2.1337    1.6760
    2.1132    1.6480
    2.0894    1.6178
    2.0873    1.6045
    2.0835    1.5901
    2.0843    1.5792
    2.0773    1.5625
    2.0702    1.5459
    2.0721    1.5366
    2.0759    1.5282
    2.0783    1.5183
    2.0725    1.5030
    2.0758    1.4949
    2.0731    1.4820
    2.0782    1.4747
    2.0818    1.4663
    2.0901    1.4613
    2.0790    1.4428
    2.0873    1.4378
    2.0891    1.4283
    2.0973    1.4232
    2.1022    1.4158
    2.0990    1.4030
    2.1069    1.3982
    2.1099    1.3902
    2.1213    1.3871
    2.1245    1.3786
    2.1177    1.3637
    2.1291    1.3606
    2.1372    1.3552
    2.1552    1.3561
    2.1549    1.3455
    2.1529    1.3338
    2.1609    1.3283
    2.1722    1.3249
    2.1768    1.3173
    2.1881    1.3137
    2.1893    1.3041
    2.1990    1.2994
    2.2103    1.2957
    2.2249    1.2939
    2.2277    1.2852
    2.2424    1.2832
    2.2554    1.2802
    2.2667    1.2762
    2.2832    1.2750
    2.2945    1.2708
    2.3092    1.2685
    2.3378    1.2735
    2.3699    1.2803];

station_tolerance = 1;
transverse_tolerance = 0.02;
points_required_for_agreement = 10;
flag_force_circle_fit = 1;
expected_radii_range = [0.2 0.4];
flag_find_only_best_agreement = 1; 
flag_use_permutations = [];

Hough_domains  = fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), -1);

fcn_geometry_plotFitDomains(Hough_domains, 10);



%% Speed test effect of flag_find_only_best_agreement


inputPoints = [
    2.1337    1.6760
    2.1132    1.6480
    2.0894    1.6178
    2.0873    1.6045
    2.0835    1.5901
    2.0843    1.5792
    2.0773    1.5625
    2.0702    1.5459
    2.0721    1.5366
    2.0759    1.5282
    2.0783    1.5183
    2.0725    1.5030
    2.0758    1.4949
    2.0731    1.4820
    2.0782    1.4747
    2.0818    1.4663
    2.0901    1.4613
    2.0790    1.4428
    2.0873    1.4378
    2.0891    1.4283
    2.0973    1.4232
    2.1022    1.4158
    2.0990    1.4030
    2.1069    1.3982
    2.1099    1.3902
    2.1213    1.3871
    2.1245    1.3786
    2.1177    1.3637
    2.1291    1.3606
    2.1372    1.3552
    2.1552    1.3561
    2.1549    1.3455
    2.1529    1.3338
    2.1609    1.3283
    2.1722    1.3249
    2.1768    1.3173
    2.1881    1.3137
    2.1893    1.3041
    2.1990    1.2994
    2.2103    1.2957
    2.2249    1.2939
    2.2277    1.2852
    2.2424    1.2832
    2.2554    1.2802
    2.2667    1.2762
    2.2832    1.2750
    2.2945    1.2708
    2.3092    1.2685
    2.3378    1.2735
    2.3699    1.2803];

transverse_tolerance = 0.1;
station_tolerance = 1;
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [1 3];
flag_use_permutations = [];



% Perform the calculation in slow mode
fig_num = -1;
flag_find_only_best_agreement = 0; 
REPS = 3; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;

    slowdomains  = ...
        fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement), (flag_use_permutations), (fig_num));

    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
fig_num = -1;
flag_find_only_best_agreement = 1; 
minTimeFast = Inf;
tic;
for i=1:REPS
    tstart = tic;

    fastdomains  = ...
        fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement),  (flag_use_permutations), (fig_num));
    
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fcn_geometry_fitHoughCircle with flag_find_only_best_agreement = 0 (slow) and with flag_use_permutations= 1 (fast):\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Slow parameters (X Y radius, all in meters): %.3f %.3f %.3f\n',slowdomains{1}.best_fit_parameters(1),slowdomains{1}.best_fit_parameters(2),slowdomains{1}.best_fit_parameters(3));
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Fast parameters (X Y radius, all in meters): %.3f %.3f %.3f\n',fastdomains{1}.best_fit_parameters(1),fastdomains{1}.best_fit_parameters(2),fastdomains{1}.best_fit_parameters(3));
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);


