%% script_test_fcn_geometry_comparePointsToPoints
% Exercises the function: fcn_geometry_comparePointsToPoints
% Revision history:
% 2024_04_14
% -- wrote the code
% -- revised from script_test_fcn_geometry_fitVectorToNPoints
%
close all;


%% Test 1: origin as reference point, 1 point as test points.
fig_num = 1;
figure(fig_num); clf;

reference_points_XY = [0 0];
test_points_XY      = [2 0];
threshold           = [];

[flag_is_similar, minimum_distance_to_each_point, indicies_of_nearest_reference_points, mean_error, max_error, std_dev_error] = fcn_geometry_comparePointsToPoints(reference_points_XY, test_points_XY, (threshold), (fig_num));


% Test sizes of variables
assert(islogical(flag_is_similar));
assert(length(minimum_distance_to_each_point(1,:))==1);
assert(length(minimum_distance_to_each_point(:,1))==length(test_points_XY(:,1)));
assert(length(indicies_of_nearest_reference_points(:,1))==length(test_points_XY(:,1)));
assert(isequal(size(mean_error),[1 1]));
assert(isequal(size(max_error),[1 1]));
assert(isequal(size(std_dev_error),[1 1]));

% Test contents of variables
assert(flag_is_similar);
assert(isequal(minimum_distance_to_each_point,2));
assert(isequal(indicies_of_nearest_reference_points,1));
assert(isequal(mean_error,2));
assert(isequal(max_error,2));
assert(isequal(std_dev_error,0));

%% Test 2: 2 points as reference point, 3 point as test points
fig_num = 2;
figure(fig_num); clf;

reference_points_XY = [0 0; 1 0; 2 0];
test_points_XY      = [0 1; 1 1; 2 2];
threshold           = [];

[flag_is_similar, minimum_distance_to_each_point, indicies_of_nearest_reference_points, mean_error, max_error, std_dev_error] = fcn_geometry_comparePointsToPoints(reference_points_XY, test_points_XY, (threshold), (fig_num));


% Test sizes of variables
assert(islogical(flag_is_similar));
assert(length(minimum_distance_to_each_point(1,:))==1);
assert(length(minimum_distance_to_each_point(:,1))==length(test_points_XY(:,1)));
assert(length(indicies_of_nearest_reference_points(:,1))==length(test_points_XY(:,1)));
assert(isequal(size(mean_error),[1 1]));
assert(isequal(size(max_error),[1 1]));
assert(isequal(size(std_dev_error),[1 1]));

% Test contents of variables
assert(flag_is_similar);
assert(isequal(minimum_distance_to_each_point,[1; 1; 2]));
assert(isequal(indicies_of_nearest_reference_points,[1; 2; 3]));
assert(isequal(round(mean_error,4),1.3333));
assert(isequal(max_error,2));
assert(isequal(round(std_dev_error,4),0.5774));

%% Test 3: origin as reference point, unit circle test points
fig_num = 3;
figure(fig_num); clf;

reference_points_XY = [0 0];
angles  = (0:10:350)'*pi/180;
test_points_XY      =[cos(angles) sin(angles)];
threshold           = [];

[flag_is_similar, minimum_distance_to_each_point, indicies_of_nearest_reference_points, mean_error, max_error, std_dev_error] = fcn_geometry_comparePointsToPoints(reference_points_XY, test_points_XY, (threshold), (fig_num));


% Test sizes of variables
assert(islogical(flag_is_similar));
assert(length(minimum_distance_to_each_point(1,:))==1);
assert(length(minimum_distance_to_each_point(:,1))==length(test_points_XY(:,1)));
assert(length(indicies_of_nearest_reference_points(:,1))==length(test_points_XY(:,1)));
assert(isequal(size(mean_error),[1 1]));
assert(isequal(size(max_error),[1 1]));
assert(isequal(size(std_dev_error),[1 1]));

% Test contents of variables
assert(flag_is_similar);
assert(all(round(minimum_distance_to_each_point,4)==1));
assert(all(indicies_of_nearest_reference_points==1));
assert(isequal(round(mean_error,4),1));
assert(isequal(max_error,1));
assert(isequal(round(std_dev_error,4),0));

%% Test 4: origin as reference point, unit circle test points, add high tolerance
fig_num = 4;
figure(fig_num); clf;

reference_points_XY = [0 0];
angles  = (0:10:350)'*pi/180;
test_points_XY      =[cos(angles) sin(angles)];
threshold           = 2;

[flag_is_similar, minimum_distance_to_each_point, indicies_of_nearest_reference_points, mean_error, max_error, std_dev_error] = fcn_geometry_comparePointsToPoints(reference_points_XY, test_points_XY, (threshold), (fig_num));


% Test sizes of variables
assert(islogical(flag_is_similar));
assert(length(minimum_distance_to_each_point(1,:))==1);
assert(length(minimum_distance_to_each_point(:,1))==length(test_points_XY(:,1)));
assert(length(indicies_of_nearest_reference_points(:,1))==length(test_points_XY(:,1)));
assert(isequal(size(mean_error),[1 1]));
assert(isequal(size(max_error),[1 1]));
assert(isequal(size(std_dev_error),[1 1]));

% Test contents of variables
assert(flag_is_similar);
assert(all(round(minimum_distance_to_each_point,4)==1));
assert(all(indicies_of_nearest_reference_points==1));
assert(isequal(round(mean_error,4),1));
assert(isequal(max_error,1));
assert(isequal(round(std_dev_error,4),0));


%% Test 5: origin as reference point, unit circle test points, add low tolerance
fig_num = 5;
figure(fig_num); clf;

reference_points_XY = [0 0];
angles  = (0:10:350)'*pi/180;
test_points_XY      =[cos(angles) sin(angles)];
threshold           = 0.5;

[flag_is_similar, minimum_distance_to_each_point, indicies_of_nearest_reference_points, mean_error, max_error, std_dev_error] = fcn_geometry_comparePointsToPoints(reference_points_XY, test_points_XY, (threshold), (fig_num));


% Test sizes of variables
assert(islogical(flag_is_similar));
assert(length(minimum_distance_to_each_point(1,:))==1);
assert(length(minimum_distance_to_each_point(:,1))==length(test_points_XY(:,1)));
assert(length(indicies_of_nearest_reference_points(:,1))==length(test_points_XY(:,1)));
assert(isequal(size(mean_error),[1 1]));
assert(isequal(size(max_error),[1 1]));
assert(isequal(size(std_dev_error),[1 1]));

% Test contents of variables
assert(~flag_is_similar);
assert(all(round(minimum_distance_to_each_point,4)==1));
assert(all(indicies_of_nearest_reference_points==1));
assert(isequal(round(mean_error,4),1));
assert(isequal(max_error,1));
assert(isequal(round(std_dev_error,4),0));

%% Test 100: many points along a line versus many points along same line
fig_num = 100;
figure(fig_num); clf;

seed_points = [2 3; 4 5];
M = 10;
sigma = 0.02;

reference_points_XY = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);
test_points_XY      = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);
threshold           = 5*sigma;

[flag_is_similar, minimum_distance_to_each_point, indicies_of_nearest_reference_points, mean_error, max_error, std_dev_error] = fcn_geometry_comparePointsToPoints(reference_points_XY, test_points_XY, (threshold), (fig_num));


% Test sizes of variables
assert(islogical(flag_is_similar));
assert(length(minimum_distance_to_each_point(1,:))==1);
assert(length(minimum_distance_to_each_point(:,1))==length(test_points_XY(:,1)));
assert(length(indicies_of_nearest_reference_points(:,1))==length(test_points_XY(:,1)));
assert(isequal(size(mean_error),[1 1]));
assert(isequal(size(max_error),[1 1]));
assert(isequal(size(std_dev_error),[1 1]));

% Test contents of variables
assert(flag_is_similar);

%% Test 101: many points along a line versus many points along same line, low tolerance
fig_num = 101;
figure(fig_num); clf;

seed_points = [2 3; 4 5];
M = 10;
sigma = 0.2;

reference_points_XY = fcn_geometry_fillLineTestPoints(seed_points, M, 0.02, -1);
test_points_XY      = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);
threshold           = 0.3;

[flag_is_similar, minimum_distance_to_each_point, indicies_of_nearest_reference_points, mean_error, max_error, std_dev_error] = fcn_geometry_comparePointsToPoints(reference_points_XY, test_points_XY, (threshold), (fig_num));


% Test sizes of variables
assert(islogical(flag_is_similar));
assert(length(minimum_distance_to_each_point(1,:))==1);
assert(length(minimum_distance_to_each_point(:,1))==length(test_points_XY(:,1)));
assert(length(indicies_of_nearest_reference_points(:,1))==length(test_points_XY(:,1)));
assert(isequal(size(mean_error),[1 1]));
assert(isequal(size(max_error),[1 1]));
assert(isequal(size(std_dev_error),[1 1]));

% Test contents of variables
assert(~flag_is_similar);

%% Test 102: many points along a line versus many points along different line, high tolerance
fig_num = 102;
figure(fig_num); clf;

seed_points = [2 3; 4 5];
M = 10;
sigma = 0.02;

reference_points_XY = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);

seed_points = [2 3; 4 4];
M = 10;
sigma = 0.02;

test_points_XY      = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);
threshold           = 1;

[flag_is_similar, minimum_distance_to_each_point, indicies_of_nearest_reference_points, mean_error, max_error, std_dev_error] = fcn_geometry_comparePointsToPoints(reference_points_XY, test_points_XY, (threshold), (fig_num));


% Test sizes of variables
assert(islogical(flag_is_similar));
assert(length(minimum_distance_to_each_point(1,:))==1);
assert(length(minimum_distance_to_each_point(:,1))==length(test_points_XY(:,1)));
assert(length(indicies_of_nearest_reference_points(:,1))==length(test_points_XY(:,1)));
assert(isequal(size(mean_error),[1 1]));
assert(isequal(size(max_error),[1 1]));
assert(isequal(size(std_dev_error),[1 1]));

% Test contents of variables
assert(flag_is_similar);

%% Test 103: many points along a line versus many points along different line, low tolerance
fig_num = 103;
figure(fig_num); clf;

seed_points = [2 3; 4 5];
M = 10;
sigma = 0.02;

reference_points_XY = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);

seed_points = [2 3; 4 4];
M = 10;
sigma = 0.02;

test_points_XY      = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);
threshold           = 0.5;

[flag_is_similar, minimum_distance_to_each_point, indicies_of_nearest_reference_points, mean_error, max_error, std_dev_error] = fcn_geometry_comparePointsToPoints(reference_points_XY, test_points_XY, (threshold), (fig_num));


% Test sizes of variables
assert(islogical(flag_is_similar));
assert(length(minimum_distance_to_each_point(1,:))==1);
assert(length(minimum_distance_to_each_point(:,1))==length(test_points_XY(:,1)));
assert(length(indicies_of_nearest_reference_points(:,1))==length(test_points_XY(:,1)));
assert(isequal(size(mean_error),[1 1]));
assert(isequal(size(max_error),[1 1]));
assert(isequal(size(std_dev_error),[1 1]));

% Test contents of variables
assert(~flag_is_similar);

%% Test 201: many points along a line versus many points along an arc
fig_num = 103;
figure(fig_num); clf;

seed_points = [2 3; 4 5];
M = 10;
sigma = 0.02;

reference_points_XY = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);

seed_points = [2 3; 2.5 4; 4 5];
M = 10;
sigma = 0.02;

test_points_XY = fcn_geometry_fillArcTestPoints(seed_points, M, sigma, -1);

threshold           = 0.5;

[flag_is_similar, minimum_distance_to_each_point, indicies_of_nearest_reference_points, mean_error, max_error, std_dev_error] = fcn_geometry_comparePointsToPoints(reference_points_XY, test_points_XY, (threshold), (fig_num));


% Test sizes of variables
assert(islogical(flag_is_similar));
assert(length(minimum_distance_to_each_point(1,:))==1);
assert(length(minimum_distance_to_each_point(:,1))==length(test_points_XY(:,1)));
assert(length(indicies_of_nearest_reference_points(:,1))==length(test_points_XY(:,1)));
assert(isequal(size(mean_error),[1 1]));
assert(isequal(size(max_error),[1 1]));
assert(isequal(size(std_dev_error),[1 1]));

% Test contents of variables
assert(flag_is_similar);

%% Test of fast implementation mode 
seed_points = [2 3; 4 5];
M = 10;
sigma = 0.02;

reference_points_XY = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, -1);

seed_points = [2 3; 2.5 4; 4 5];
M = 10;
sigma = 0.02;

test_points_XY = fcn_geometry_fillArcTestPoints(seed_points, M, sigma, -1);

threshold           = 0.5;



% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    [flag_is_similar, minimum_distance_to_each_point, indicies_of_nearest_reference_points, mean_error, max_error, std_dev_error] = fcn_geometry_comparePointsToPoints(reference_points_XY, test_points_XY, (threshold), (fig_num));
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
fig_num = -1;
minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    [flag_is_similar, minimum_distance_to_each_point, indicies_of_nearest_reference_points, mean_error, max_error, std_dev_error] = fcn_geometry_comparePointsToPoints(reference_points_XY, test_points_XY, (threshold), (fig_num));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_comparePointsToPoints:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);
%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [root_point, unit_vector] = fcn_geometry_comparePointsToPoints(points,fig_num);
    fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));
end
