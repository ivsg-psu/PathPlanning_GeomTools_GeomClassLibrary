% script_test_fcn_geometry_findAngleBetweenAngles
% Exercises the function: fcn_geometry_findAngleBetweenAngles
% Revision history:
% 2024_01_07
% -- wrote the code

close all;


%% Test 1: a basic test 
fig_num = 1;


start_angle_in_radians = 0*pi/180;
end_angle_in_radians = 90*pi/180;
angles_to_test_in_radians = 45*pi/180;
direction = 1;

[isAngleBetween]  = fcn_geometry_findAngleBetweenAngles(start_angle_in_radians, end_angle_in_radians, direction, angles_to_test_in_radians, (fig_num));

assert(isequal(isAngleBetween,1));

%% Test 2: a basic test - vectorized
fig_num = 2;


start_angle_in_radians = 0*pi/180;
change_in_angle = 90*pi/180;
end_angle_in_radians = start_angle_in_radians+change_in_angle;
angles_to_test_in_radians = (start_angle_in_radians:(start_angle_in_radians+2*pi))';
direction = 1;

[isAngleBetween]  = fcn_geometry_findAngleBetweenAngles(start_angle_in_radians, end_angle_in_radians, direction, angles_to_test_in_radians, (fig_num));
correct_results = (angles_to_test_in_radians-start_angle_in_radians)<=end_angle_in_radians;
assert(isequal(isAngleBetween,correct_results));

%% Test 3: systematic testing
fig_num = 3;
figure(fig_num); clf;
set(fig_num,'UserData',[]);

increment_angle = 15;
increment_angle_in_radians = (increment_angle*pi/180);
changes_in_angles = (0:increment_angle:(360))'*pi/180;
start_angles_in_radians = (0:increment_angle:360)'*pi/180;
direction = 1;

for ith_start = 1:length(start_angles_in_radians)
    start_angle_in_radians = start_angles_in_radians(ith_start);
    for ith_test = 1:length(changes_in_angles)
        change_in_angle = changes_in_angles(ith_test);
        end_angle_in_radians = start_angle_in_radians+change_in_angle;
        angles_to_test_in_radians = (start_angle_in_radians:increment_angle_in_radians:(start_angle_in_radians + 2*pi - increment_angle_in_radians))';

        [isAngleBetween]  = fcn_geometry_findAngleBetweenAngles(start_angle_in_radians, end_angle_in_radians, direction, angles_to_test_in_radians, (fig_num));
        correct_results = round((angles_to_test_in_radians-start_angle_in_radians),8)<=round(mod(change_in_angle,2*pi),8);
        drawnow;
        assert(isequal(isAngleBetween,correct_results));
    end
end

direction = -1;
for ith_start = 1:length(start_angles_in_radians)
    start_angle_in_radians = start_angles_in_radians(ith_start);
    for ith_test = 1:length(changes_in_angles)
        change_in_angle = changes_in_angles(ith_test);
        end_angle_in_radians = start_angle_in_radians+change_in_angle;
        angles_to_test_in_radians = (start_angle_in_radians:increment_angle_in_radians:(start_angle_in_radians + 2*pi - increment_angle_in_radians))';

        [isAngleBetween]  = fcn_geometry_findAngleBetweenAngles(start_angle_in_radians, end_angle_in_radians, direction, angles_to_test_in_radians, (fig_num));
        correct_results = round((angles_to_test_in_radians-start_angle_in_radians),8)>=round(mod(change_in_angle,2*pi),8);
        drawnow;
        assert(isequal(isAngleBetween,correct_results));
    end
end

%% Test of fast mode

start_angle_in_radians = 0*pi/180;
change_in_angle = 90*pi/180;
end_angle_in_radians = start_angle_in_radians+change_in_angle;
angles_to_test_in_radians = (start_angle_in_radians:(start_angle_in_radians+2*pi))';
direction = 1;


% Perform the calculation in slow mode
REPS = 10; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    [isAngleBetween]  = fcn_geometry_findAngleBetweenAngles(start_angle_in_radians, end_angle_in_radians, direction, angles_to_test_in_radians, ([]));
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    [isAngleBetween]  = fcn_geometry_findAngleBetweenAngles(start_angle_in_radians, end_angle_in_radians, direction, angles_to_test_in_radians, (-1));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_findAngleBetweenAngles:\n');
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
