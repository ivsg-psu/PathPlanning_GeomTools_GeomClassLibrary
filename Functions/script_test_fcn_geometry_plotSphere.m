% script_test_fcn_geometry_plotSphere
% Tests fcn_geometry_plotSphere

% Revision history:
%      2024_01_23 by S. Brennan
%      -- Edited for new function

close all

%% BASIC example for one circle
fig_num = 1;
figure(fig_num); clf;

center = [1 3 5];
radius = [2]; %#ok<*NBRAK>
color_vector = [];
fcn_geometry_plotSphere(center, radius, color_vector, fig_num);

%% BASIC example for one circle, with color
fig_num = 2;
figure(fig_num); clf;

center = [1 3 5];
radius = [2]; %#ok<*NBRAK>
color_vector = [0 0 1];
fcn_geometry_plotSphere(center, radius, color_vector, fig_num);

%% BASIC example for one circle, with color, saving data
fig_num = 21;
figure(fig_num); clf;

center = [1 3 5];
radius = [2]; %#ok<*NBRAK>
color_vector = [0 0 1];
XYZ_data = fcn_geometry_plotSphere(center, radius, color_vector, fig_num);

assert(length(XYZ_data(1,:))==3);

%% BASIC example for multiple spheres
fig_num = 3;
figure(fig_num); clf;

centers = [1 3 0; 2 4 5];
radii = [0.4; 1];
color_vector = [0 1 0];
XYZ_data = fcn_geometry_plotSphere(centers,radii,color_vector,fig_num);
assert(length(XYZ_data)==2);
assert(length(XYZ_data{1}(1,:))==3);
assert(length(XYZ_data{2}(1,:))==3);


%% BASIC example for multiple spheres with colors
% % Does not work - does not change color
% fig_num = 4;
% figure(fig_num); clf;
% 
% centers = [1 3 0; 2 4 0];
% radii = [0.4; 1];
% color_vector = [0 0 1; 1 0 0];
% fcn_geometry_plotSphere(centers,radii,color_vector,fig_num);


%% Test of fast implementation mode 

centers = [1 3 0; 2 4 5];
radii = [0.4; 1];
color_vector = [0 1 0];

% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    fcn_geometry_plotSphere(centers,radii,color_vector,fig_num);
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
    fcn_geometry_plotSphere(centers,radii,color_vector,fig_num);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_plotSphere:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

%% Break cases follow
% - these are ones that intentionally crash the code by passing invalid
% arguments
if 1==0
%% BREAK CASES 1 - break on centers
fig_num = 999;

centers  = [1 2; 2 4; 3 5];
radii = [3; 4];
fcn_geometry_plotSphere(centers,radii,[0.1 0.1 1])

end
