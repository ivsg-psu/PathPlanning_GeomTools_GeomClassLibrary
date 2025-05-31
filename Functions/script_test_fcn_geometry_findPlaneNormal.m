% script_test_fcn_geometry_findPlaneNormal
% Exercises the function: fcn_geometry_findPlaneNormal
% Revision history:
% 2025_05_30
% -- wrote the code

close all;


%% Demonstration Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  _____                                 _             _   _               ______                           _
% |  __ \                               | |           | | (_)             |  ____|                         | |
% | |  | | ___ _ __ ___   ___  _ __  ___| |_ _ __ __ _| |_ _  ___  _ __   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___
% | |  | |/ _ \ '_ ` _ \ / _ \| '_ \/ __| __| '__/ _` | __| |/ _ \| '_ \  |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|
% | |__| |  __/ | | | | | (_) | | | \__ \ |_| | | (_| | |_| | (_) | | | | | |____ >  < (_| | | | | | | |_) | |  __/\__ \
% |_____/ \___|_| |_| |_|\___/|_| |_|___/\__|_|  \__,_|\__|_|\___/|_| |_| |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/
%                                                                                                    | |
%                                                                                                    |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Demonstration%20Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Demonstration case 1: many points in XY plane
fig_num = 0001;
figure(fig_num);
clf;

points = [
    0 0 0;
    1 0 0; 
    2 3 0; 
    4 5 0; 
    2 7 1; 
    0 8 1; 
    -1 7.5 1; 
    -2 6 0; 
    -4 5 0; 
    -2 3 0; 
    -1 -1 0]*2;

[unit_normal_vector, base_point, flags_in_directional_agreement, flags_in_magnitude_agreement] = fcn_geometry_findPlaneNormal(points,(fig_num));

view(70, 20);

% Check variable types
assert(isnumeric(unit_normal_vector));
assert(isnumeric(base_point));
assert(islogical(flags_in_directional_agreement));
assert(islogical(flags_in_magnitude_agreement));

% Check variable lengths
Npoints = length(points(:,1));
assert(isequal(size(unit_normal_vector),[1 3]));
assert(isequal(size(base_point),[1 3]));
assert(isequal(size(flags_in_directional_agreement),[Npoints 1]));
assert(isequal(size(flags_in_magnitude_agreement),[Npoints 1]));

% Check variable values
vectorLength = sum(unit_normal_vector.^2,2).^0.5;
assert(isequal(round(vectorLength,4),1));
assert(isequal(round(base_point,4),round(mean(points,1),4)));
assert(~all(flags_in_directional_agreement==1));
assert(~all(flags_in_magnitude_agreement==1));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic testing examples in 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____            _        _______        _   _                ______                           _                       ____  _____
% |  _ \          (_)      |__   __|      | | (_)              |  ____|                         | |                     |___ \|  __ \
% | |_) | __ _ ___ _  ___     | | ___  ___| |_ _ _ __   __ _   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___    ______    __) | |  | |
% |  _ < / _` / __| |/ __|    | |/ _ \/ __| __| | '_ \ / _` |  |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|  |______|  |__ <| |  | |
% | |_) | (_| \__ \ | (__     | |  __/\__ \ |_| | | | | (_| |  | |____ >  < (_| | | | | | | |_) | |  __/\__ \            ___) | |__| |
% |____/ \__,_|___/_|\___|    |_|\___||___/\__|_|_| |_|\__, |  |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/           |____/|_____/
%                                                       __/ |                             | |
%                                                      |___/                              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20Testing%20%20Examples%20%20-%203D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic test case 1: basic plane with 3 points in XY plane
fig_num = 1001;
figure(fig_num);
clf;

points = [0 0 0; 1 0 0; 1 1 0]*3;

[unit_normal_vector, base_point, flags_in_directional_agreement, flags_in_magnitude_agreement] = fcn_geometry_findPlaneNormal(points,(fig_num));

% Check variable types
assert(isnumeric(unit_normal_vector));
assert(isnumeric(base_point));
assert(islogical(flags_in_directional_agreement));
assert(islogical(flags_in_magnitude_agreement));

% Check variable lengths
Npoints = length(points(:,1));
assert(isequal(size(unit_normal_vector),[1 3]));
assert(isequal(size(base_point),[1 3]));
assert(isequal(size(flags_in_directional_agreement),[Npoints 1]));
assert(isequal(size(flags_in_magnitude_agreement),[Npoints 1]));

% Check variable values
vectorLength = sum(unit_normal_vector.^2,2).^0.5;
assert(isequal(round(vectorLength,4),1));
assert(isequal(round(base_point,4),round(mean(points,1),4)));
assert(all(flags_in_directional_agreement==1));
assert(all(flags_in_magnitude_agreement==1));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic test case 2: basic plane with 4 points in XY plane
fig_num = 1002;
figure(fig_num);
clf;

points = [0 0 0; 0 1 0; 1 0 0; 0 0 0];

[unit_normal_vector, base_point, flags_in_directional_agreement, flags_in_magnitude_agreement] = fcn_geometry_findPlaneNormal(points,(fig_num));

% Check variable types
assert(isnumeric(unit_normal_vector));
assert(isnumeric(base_point));
assert(islogical(flags_in_directional_agreement));
assert(islogical(flags_in_magnitude_agreement));

% Check variable lengths
Npoints = length(points(:,1));
assert(isequal(size(unit_normal_vector),[1 3]));
assert(isequal(size(base_point),[1 3]));
assert(isequal(size(flags_in_directional_agreement),[Npoints 1]));
assert(isequal(size(flags_in_magnitude_agreement),[Npoints 1]));

% Check variable values
vectorLength = sum(unit_normal_vector.^2,2).^0.5;
assert(isequal(round(vectorLength,4),1));
assert(isequal(round(base_point,4),round(mean(points,1),4)));
assert(all(flags_in_directional_agreement==1));
assert(all(flags_in_magnitude_agreement==1));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Basic test case 3: basic plane with many points in XY plane
fig_num = 1003;
figure(fig_num);
clf;

points = [
    0 0 0;
    1 0 0; 
    2 3 0; 
    4 5 0; 
    2 7 1; 
    0 8 1; 
    -1 7.5 1; 
    -2 6 0; 
    -4 5 0; 
    -2 3 0; 
    -1 -1 0]*2;

[unit_normal_vector, base_point, flags_in_directional_agreement, flags_in_magnitude_agreement] = fcn_geometry_findPlaneNormal(points,(fig_num));

view(70, 20);

% Check variable types
assert(isnumeric(unit_normal_vector));
assert(isnumeric(base_point));
assert(islogical(flags_in_directional_agreement));
assert(islogical(flags_in_magnitude_agreement));

% Check variable lengths
Npoints = length(points(:,1));
assert(isequal(size(unit_normal_vector),[1 3]));
assert(isequal(size(base_point),[1 3]));
assert(isequal(size(flags_in_directional_agreement),[Npoints 1]));
assert(isequal(size(flags_in_magnitude_agreement),[Npoints 1]));

% Check variable values
vectorLength = sum(unit_normal_vector.^2,2).^0.5;
assert(isequal(round(vectorLength,4),1));
assert(isequal(round(base_point,4),round(mean(points,1),4)));
assert(~all(flags_in_directional_agreement==1));
assert(~all(flags_in_magnitude_agreement==1));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic test case 4: basic Vertial plane with many points in XZ plane
fig_num = 1004;
figure(fig_num);
clf;

old_points = [
    0 0 0;
    1 0 0; 
    2 3 0; 
    4 5 0; 
    2 7 1; 
    0 8 1; 
    -1 7.5 1; 
    -2 6 0; 
    -4 5 0; 
    -2 3 0; 
    -1 -1 0]*2;
points = [old_points(:,1) old_points(:,3) old_points(:,2)];

[unit_normal_vector, base_point, flags_in_directional_agreement, flags_in_magnitude_agreement] = fcn_geometry_findPlaneNormal(points,(fig_num));

view(70, 20);

% Check variable types
assert(isnumeric(unit_normal_vector));
assert(isnumeric(base_point));
assert(islogical(flags_in_directional_agreement));
assert(islogical(flags_in_magnitude_agreement));

% Check variable lengths
Npoints = length(points(:,1));
assert(isequal(size(unit_normal_vector),[1 3]));
assert(isequal(size(base_point),[1 3]));
assert(isequal(size(flags_in_directional_agreement),[Npoints 1]));
assert(isequal(size(flags_in_magnitude_agreement),[Npoints 1]));

% Check variable values
vectorLength = sum(unit_normal_vector.^2,2).^0.5;
assert(isequal(round(vectorLength,4),1));
assert(isequal(round(base_point,4),round(mean(points,1),4)));
assert(~all(flags_in_directional_agreement==1));
assert(~all(flags_in_magnitude_agreement==1));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Fast Mode Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ______        _     __  __           _        _______        _
% |  ____|      | |   |  \/  |         | |      |__   __|      | |
% | |__ __ _ ___| |_  | \  / | ___   __| | ___     | | ___  ___| |_ ___
% |  __/ _` / __| __| | |\/| |/ _ \ / _` |/ _ \    | |/ _ \/ __| __/ __|
% | | | (_| \__ \ |_  | |  | | (_) | (_| |  __/    | |  __/\__ \ |_\__ \
% |_|  \__,_|___/\__| |_|  |_|\___/ \__,_|\___|    |_|\___||___/\__|___/
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Fast%20Mode%20Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic example - NO FIGURE
fig_num = 9901;
figure(fig_num);
close(fig_num);

points = [
    0 0 0;
    1 0 0; 
    2 3 0; 
    4 5 0; 
    2 7 1; 
    0 8 1; 
    -1 7.5 1; 
    -2 6 0; 
    -4 5 0; 
    -2 3 0; 
    -1 -1 0]*2;

[unit_normal_vector, base_point, flags_in_directional_agreement, flags_in_magnitude_agreement] = fcn_geometry_findPlaneNormal(points,([]));

view(70, 20);

% Check variable types
assert(isnumeric(unit_normal_vector));
assert(isnumeric(base_point));
assert(islogical(flags_in_directional_agreement));
assert(islogical(flags_in_magnitude_agreement));

% Check variable lengths
Npoints = length(points(:,1));
assert(isequal(size(unit_normal_vector),[1 3]));
assert(isequal(size(base_point),[1 3]));
assert(isequal(size(flags_in_directional_agreement),[Npoints 1]));
assert(isequal(size(flags_in_magnitude_agreement),[Npoints 1]));

% Check variable values
vectorLength = sum(unit_normal_vector.^2,2).^0.5;
assert(isequal(round(vectorLength,4),1));
assert(isequal(round(base_point,4),round(mean(points,1),4)));
assert(~all(flags_in_directional_agreement==1));
assert(~all(flags_in_magnitude_agreement==1));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Basic example of vertex calculation - non-normal wall shrinking, NO FIGURE, FAST MODE
fig_num = 9902;
figure(fig_num);
close(fig_num);

points = [
    0 0 0;
    1 0 0; 
    2 3 0; 
    4 5 0; 
    2 7 1; 
    0 8 1; 
    -1 7.5 1; 
    -2 6 0; 
    -4 5 0; 
    -2 3 0; 
    -1 -1 0]*2;

[unit_normal_vector, base_point, flags_in_directional_agreement, flags_in_magnitude_agreement] = fcn_geometry_findPlaneNormal(points,(-1));

view(70, 20);

% Check variable types
assert(isnumeric(unit_normal_vector));
assert(isnumeric(base_point));
assert(islogical(flags_in_directional_agreement));
assert(islogical(flags_in_magnitude_agreement));

% Check variable lengths
Npoints = length(points(:,1));
assert(isequal(size(unit_normal_vector),[1 3]));
assert(isequal(size(base_point),[1 3]));
assert(isequal(size(flags_in_directional_agreement),[Npoints 1]));
assert(isequal(size(flags_in_magnitude_agreement),[Npoints 1]));

% Check variable values
vectorLength = sum(unit_normal_vector.^2,2).^0.5;
assert(isequal(round(vectorLength,4),1));
assert(isequal(round(base_point,4),round(mean(points,1),4)));
assert(~all(flags_in_directional_agreement==1));
assert(~all(flags_in_magnitude_agreement==1));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 9903;
figure(fig_num);
close(fig_num);
rng(1823);

points = [
    0 0 0;
    1 0 0; 
    2 3 0; 
    4 5 0; 
    2 7 1; 
    0 8 1; 
    -1 7.5 1; 
    -2 6 0; 
    -4 5 0; 
    -2 3 0; 
    -1 -1 0]*2;


% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    [unit_normal_vector, base_point, flags_in_directional_agreement, flags_in_magnitude_agreement] = fcn_geometry_findPlaneNormal(points,fig_num);

    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
slow_method = toc;

% Perform the operation in fast mode
fig_num = -1;
minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    [unit_normal_vector, base_point, flags_in_directional_agreement, flags_in_magnitude_agreement] = fcn_geometry_findPlaneNormal(points,fig_num);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
fast_method = toc;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_findPlaneNormal:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',slow_method/REPS);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',fast_method/REPS);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',slow_method/fast_method);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

% Plot results as bar chart
figure(373737);
clf;
X = categorical({'Normal mode','Fast mode'});
X = reordercats(X,{'Normal mode','Fast mode'}); % Forces bars to appear in this exact order, not alphabetized
Y = [slow_method fast_method ]*1000/REPS;
bar(X,Y)
ylabel('Execution time (Milliseconds)')


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [root_point, unit_vector] = fcn_geometry_findPlaneNormal(points,fig_num);
    fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));
end