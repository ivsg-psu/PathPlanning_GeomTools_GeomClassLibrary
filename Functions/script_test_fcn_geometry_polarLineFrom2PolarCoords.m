% script_test_fcn_geometry_polarLineFrom2PolarCoords
% Exercises the function: fcn_geometry_polarLineFrom2PolarCoords

% Revision history:
% 2021_05_27
% -- wrote the code

close all;



%% script_test_fcn_geometry_polarLineFrom2PolarCoords
fig_num = 1;

start_point_cart = [1 1];
end_point_cart   = [4 3];

[theta1,r1] = cart2pol(start_point_cart(:,1),start_point_cart(:,2));
[theta2,r2]   = cart2pol(end_point_cart(:,1),end_point_cart(:,2));
points = [theta1,r1;theta2,r2];

% Calculate the line
[phi,rho] = fcn_geometry_polarLineFrom2PolarCoords(points,fig_num);
axis([0 5 0 5]);
assert(isequal(round(phi,4),-0.9828));
assert(isequal(round(rho,4),-0.2774));

%% Test of fast implementation mode 


start_point_cart = [1 1];
end_point_cart   = [4 3];

[theta1,r1] = cart2pol(start_point_cart(:,1),start_point_cart(:,2));
[theta2,r2]   = cart2pol(end_point_cart(:,1),end_point_cart(:,2));
points = [theta1,r1;theta2,r2];

% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    [phi,rho] = fcn_geometry_polarLineFrom2PolarCoords(points, (fig_num));
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
    [phi,rho] = fcn_geometry_polarLineFrom2PolarCoords(points, (fig_num));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_fitVectorToNPoints:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

