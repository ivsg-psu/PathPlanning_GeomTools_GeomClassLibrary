% script_test_fcn_geometry_findAgreementsOfPointsToLineVector
% Exercises the function: fcn_geometry_findAgreementsOfPointsToLineVector
% Revision history:
% 2024_01_14
% -- wrote the code

close all;

%% Fill in test data
rng(383);

% Fill test data - 3 segments
seed_points = [2 3; 4 5; 7 0; 9 5]; 
M = 10;
sigma = 0.2;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 4;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption));


% Shuffle points?
test_points = fcn_geometry_shufflePointOrdering(test_points);

assert(length(test_points)>1);
assert(length(test_points(1,:)) == 2);

%% Test 1: a basic test of line segment fitting, specifying index-type base_point_index
fig_num = 1;
figure(fig_num); clf;

rng(383);

% Fill test data - 3 segments
seed_points = [2 3; 4 5; 7 0; 9 5]; 
M = 10;
sigma = 0.2;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 4;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption));


% Shuffle points?
test_points = fcn_geometry_shufflePointOrdering(test_points);

base_point_index = 1;
base_point = test_points(base_point_index,:);
end_point = [4 5];
unit_projection_vector = fcn_geometry_calcUnitVector(end_point-base_point,-1);
transverse_tolerance = 0.2;
station_tolerance = 2;

[agreement_indicies,station_distances] = fcn_geometry_findAgreementsOfPointsToLineVector( test_points, unit_projection_vector, base_point_index, transverse_tolerance, station_tolerance, (fig_num));

assert(length(agreement_indicies)>1);
assert(length(station_distances(1,:)) == 2);


%% Test 2: a basic test of line segment fitting, specifying point-type base_point_index
fig_num = 2;
figure(fig_num); clf;

rng(383);

% Fill test data - 3 segments
seed_points = [2 3; 4 5; 7 0; 9 5]; 
M = 10;
sigma = 0.2;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 4;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption));


% Shuffle points?
test_points = fcn_geometry_shufflePointOrdering(test_points);

base_point_index = [2 3];
end_point = [4 5];
unit_projection_vector = fcn_geometry_calcUnitVector(end_point-base_point_index,-1);
transverse_tolerance = 0.2;
station_tolerance = 2;

agreement_indicies = fcn_geometry_findAgreementsOfPointsToLineVector( test_points, unit_projection_vector, base_point_index, transverse_tolerance, station_tolerance, (fig_num));

assert(length(agreement_indicies)>1);

%% Test 3: a basic test of line segment fitting, showing unit vector effect
fig_num = 3;
figure(fig_num); clf;

rng(383);

% Fill test data - 3 segments
seed_points = [2 3; 4 5; 7 0; 9 5]; 
M = 10;
sigma = 0.2;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 4;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption));


% Shuffle points?
test_points = fcn_geometry_shufflePointOrdering(test_points);

base_point_index = [2 3];
end_point = [3 3];
unit_projection_vector = fcn_geometry_calcUnitVector(end_point-base_point_index,-1);
transverse_tolerance = 0.2;
station_tolerance = 2;

agreement_indicies = fcn_geometry_findAgreementsOfPointsToLineVector( test_points, unit_projection_vector, base_point_index, transverse_tolerance, station_tolerance, (fig_num));

assert(length(agreement_indicies)>1);

%% Test 4: a basic test of line segment fitting, showing unit vector effect
fig_num = 4;
figure(fig_num); clf;

rng(383);

% Fill test data - 3 segments
seed_points = [2 3; 4 5; 7 0; 9 5]; 
M = 10;
sigma = 0.2;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 4;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption));


% Shuffle points?
test_points = fcn_geometry_shufflePointOrdering(test_points);

base_point_index = [2 3];
end_point = [9 4];
unit_projection_vector = fcn_geometry_calcUnitVector(end_point-base_point_index,-1);
transverse_tolerance = 0.2;
station_tolerance = 2;

agreement_indicies = fcn_geometry_findAgreementsOfPointsToLineVector( test_points, unit_projection_vector, base_point_index, transverse_tolerance, station_tolerance, (fig_num));

assert(length(agreement_indicies)>1);

%% Test 5: a basic test of line segment fitting, showing transverse_tolerance effect
fig_num = 5;
figure(fig_num); clf;

rng(383);

% Fill test data - 3 segments
seed_points = [2 3; 4 5; 7 0; 9 5]; 
M = 10;
sigma = 0.2;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 4;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption));


% Shuffle points?
test_points = fcn_geometry_shufflePointOrdering(test_points);

base_point_index = [2 3];
end_point = [9 4];
unit_projection_vector = fcn_geometry_calcUnitVector(end_point-base_point_index,-1);
transverse_tolerance = 1;
station_tolerance = 2;

agreement_indicies = fcn_geometry_findAgreementsOfPointsToLineVector( test_points, unit_projection_vector, base_point_index, transverse_tolerance, station_tolerance, (fig_num));

assert(length(agreement_indicies)>1);

%% Test 6: a basic test of line segment fitting, showing station_tolerance effect
% NOTE: the lines indicate that the stations are NOT re-sorted if
% station_tolerance is not specified
fig_num = 6;
figure(fig_num); clf;

rng(383);

% Fill test data - 3 segments
seed_points = [2 3; 4 5; 7 0; 9 5]; 
M = 10;
sigma = 0.2;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 4;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption));


% Shuffle points?
test_points = fcn_geometry_shufflePointOrdering(test_points);

base_point_index = [2 3];
end_point = [9 4];
unit_projection_vector = fcn_geometry_calcUnitVector(end_point-base_point_index,-1);
transverse_tolerance = 1;
station_tolerance = [];

agreement_indicies = fcn_geometry_findAgreementsOfPointsToLineVector( test_points, unit_projection_vector, base_point_index, transverse_tolerance, station_tolerance, (fig_num));

assert(length(agreement_indicies)>1);

%% Test of fast mode

rng(383);

% Fill test data - 3 segments
seed_points = [2 3; 4 5; 7 0; 9 5]; 
M = 10;
sigma = 0.2;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 4;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption));


% Shuffle points?
test_points = fcn_geometry_shufflePointOrdering(test_points);

% Set inputs to a "normal" mode
base_point_index = 1;
end_point = [9 4];
unit_projection_vector = fcn_geometry_calcUnitVector(end_point-base_point_index,-1);
transverse_tolerance = 1;
station_tolerance = 2;

% Perform the calculation in slow mode
fig_num = [];
REPS = 10; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    agreement_indicies = fcn_geometry_findAgreementsOfPointsToLineVector( test_points, unit_projection_vector, base_point_index, transverse_tolerance, station_tolerance, (fig_num));
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
    agreement_indicies = fcn_geometry_findAgreementsOfPointsToLineVector( test_points, unit_projection_vector, base_point_index, transverse_tolerance, station_tolerance, (fig_num));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_findAgreementsOfPointsToLineVector:\n');
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
    [slope,intercept] = fcn_geometry_findAgreementsOfPointsToLineVector(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end