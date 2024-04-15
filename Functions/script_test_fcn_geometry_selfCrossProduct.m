% script_test_fcn_geometry_selfCrossProduct.m
% Tests fcn_geometry_selfCrossProduct.m

% Revision history:
%      2021_04_25:
%      -- first write of the code copying functionality from
%      fcn_FastestTraversal_checkInputsToFunctions

close all;

%% BASIC example - find the  cross-products for positive bend
fig_num = 1;
path = [0 0; 1 1; 0 2];
[cross_products] = ...
    fcn_geometry_selfCrossProduct(...
    path, fig_num);

assert(isequal(cross_products,2));

%% BASIC example - find all the cross-products for negative bend
fig_num = 2;
path = [0 0; 1 1; 2 0];
[cross_products] = ...
    fcn_geometry_selfCrossProduct(...
    path, fig_num);

assert(isequal(cross_products,-2));

%% BASIC example -
% find all the cross-products for many bends
fig_num = 3;
path = [0 0; 1 1; 0 2; 2 4; 4 2; 6 2; 2 7];
[cross_products] = ...
    fcn_geometry_selfCrossProduct(...
    path, fig_num);
assert(isequal(cross_products,[ 2    -4    -8     4    10]'));

%% Test of fast implementation mode 

path = [0 0; 1 1; 0 2; 2 4; 4 2; 6 2; 2 7];

% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    [cross_products] = ...
    fcn_geometry_selfCrossProduct(...
    path, fig_num);
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
    [cross_products] = ...
    fcn_geometry_selfCrossProduct(...
    path, fig_num);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_selfCrossProduct:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);


assert(averageTimeFast<averageTimeSlow);


if 1==0

    %% ERROR example - find cross product for straight
    fig_num = 4;
    path = [0 0; 1 1; 2 2];
    [~,~] = ...
        fcn_geometry_selfCrossProduct(...
        path, fig_num);

    %% ERROR example - find cross product for bent back on self
    fig_num = 5;
    path = [0 0; 1 1; 0 0];
    [~,~] = ...
        fcn_geometry_selfCrossProduct(...
        path, fig_num);

    %% ERROR example - find cross product for zero-length segment
    fig_num = 6;
    path = [1 1; 1 1; 0 0];
    [~,~] = ...
        fcn_geometry_selfCrossProduct(...
        path, fig_num);

end

