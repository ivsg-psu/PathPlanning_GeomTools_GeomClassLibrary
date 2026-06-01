% script_test_fcn_geometry_selfCrossProduct.m
% Tests fcn_geometry_selfCrossProduct.m

% Revision history:
%      2021_04_25:
%      -- first write of the code copying functionality from
%      fcn_FastestTraversal_checkInputsToFunctions

close all;

%% BASIC example - find the  cross-products for positive bend
figNum = 1;
path = [0 0; 1 1; 0 2];
[cross_products] = ...
    fcn_geometry_selfCrossProduct(...
    path, figNum);

assert(isequal(cross_products,2));

%% BASIC example - find all the cross-products for negative bend
figNum = 2;
path = [0 0; 1 1; 2 0];
[cross_products] = ...
    fcn_geometry_selfCrossProduct(...
    path, figNum);

assert(isequal(cross_products,-2));

%% BASIC example -
% find all the cross-products for many bends
figNum = 3;
path = [0 0; 1 1; 0 2; 2 4; 4 2; 6 2; 2 7];
[cross_products] = ...
    fcn_geometry_selfCrossProduct(...
    path, figNum);
assert(isequal(cross_products,[ 2    -4    -8     4    10]'));

%% Test of fast implementation mode 

path = [0 0; 1 1; 0 2; 2 4; 4 2; 6 2; 2 7];

% Perform the calculation in slow mode
figNum = [];
REPS = 100; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    [cross_products] = ...
    fcn_geometry_selfCrossProduct(...
    path, figNum);
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
figNum = -1;
minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    [cross_products] = ...
    fcn_geometry_selfCrossProduct(...
    path, figNum);
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
    figNum = 4;
    path = [0 0; 1 1; 2 2];
    [~,~] = ...
        fcn_geometry_selfCrossProduct(...
        path, figNum);

    %% ERROR example - find cross product for bent back on self
    figNum = 5;
    path = [0 0; 1 1; 0 0];
    [~,~] = ...
        fcn_geometry_selfCrossProduct(...
        path, figNum);

    %% ERROR example - find cross product for zero-length segment
    figNum = 6;
    path = [1 1; 1 1; 0 0];
    [~,~] = ...
        fcn_geometry_selfCrossProduct(...
        path, figNum);

end

